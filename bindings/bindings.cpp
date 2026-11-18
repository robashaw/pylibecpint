#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "api.hpp"
#include "config.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"

#ifndef LIBECPINT_VERSION_STR
#define LIBECPINT_VERSION_STR "unknown"
#endif

namespace py = pybind11;
using namespace libecpint;

namespace {

template <typename T>
using carray = py::array_t<T, py::array::c_style | py::array::forcecast>;

// Copies a TwoIndex into a new 2D numpy array
py::array_t<double> to_numpy(const TwoIndex<double>& m) {
  return py::array_t<double>({m.dims[0], m.dims[1]}, m.data.data());
}

py::list to_numpy_list(const std::vector<TwoIndex<double>>& ms) {
  py::list result;
  for (const auto& m : ms) result.append(to_numpy(m));
  return result;
}

void check_size(py::ssize_t actual, py::ssize_t expected, const char* name) {
  if (actual != expected)
    throw py::value_error(std::string(name) + " has " + std::to_string(actual) +
                          " elements, expected " + std::to_string(expected));
}

py::ssize_t total_length(const carray<int>& lengths) {
  const int* p = lengths.data();
  py::ssize_t total = 0;
  for (py::ssize_t i = 0; i < lengths.size(); i++) {
    if (p[i] < 0) throw py::value_error("shell_lengths must be non-negative");
    total += p[i];
  }
  return total;
}

// The C++ library only asserts these preconditions, and asserts are compiled out of release
// builds, so they are checked here to raise instead of reading out of bounds
void require_init(const ECPIntegrator& self) {
  if (!self.ecpint)
    throw std::runtime_error("ECPIntegrator.init must be called before computing anything");
}

void require_deriv(const ECPIntegrator& self, int order, const char* what) {
  require_init(self);
  if (self.deriv < order)
    throw std::runtime_error(std::string(what) + " requires init(deriv=" + std::to_string(order) +
                             ") or higher, but the integrator was initialised with deriv=" +
                             std::to_string(self.deriv));
}

void check_density(const ECPIntegrator& self, const carray<double>& density) {
  if (density.ndim() != 2 || density.shape(0) != self.ncart || density.shape(1) != self.ncart)
    throw py::value_error("density must have shape (ncart, ncart) = (" +
                          std::to_string(self.ncart) + ", " + std::to_string(self.ncart) + ")");
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Python bindings for libecpint";
  m.attr("MAX_L") = LIBECPINT_MAX_L;
  m.attr("__libecpint_version__") = LIBECPINT_VERSION_STR;

  // --- GaussianECP ---
  py::class_<GaussianECP>(m, "GaussianECP")
      .def(py::init<int, int, double, double>(), py::arg("n"), py::arg("l"), py::arg("a"),
           py::arg("d"))
      .def(py::init<>())
      .def_readonly("n", &GaussianECP::n)
      .def_readonly("l", &GaussianECP::l)
      .def_readonly("a", &GaussianECP::a)
      .def_readonly("d", &GaussianECP::d);

  // --- ECP ---
  py::class_<ECP>(m, "ECP")
      .def(py::init<>())
      .def(py::init([](carray<double> center) {
             check_size(center.size(), 3, "center");
             return ECP(center.data());
           }),
           py::arg("center"))
      .def("addPrimitive", &ECP::addPrimitive, py::arg("n"), py::arg("l"), py::arg("a"),
           py::arg("d"), py::arg("needSort") = true,
           "Add a primitive d r^(n-2) exp(-a r^2) to channel l. Raises ValueError if l is outside "
           "0..MAX_L or n outside 0..4.")
      .def("center",
           [](const ECP& self) {
             const double* c = self.center();
             return py::array_t<double>(3, c);
           })
      .def("setPos", &ECP::setPos, py::arg("x"), py::arg("y"), py::arg("z"))
      .def("sort", &ECP::sort)
      .def(
          "getGaussian",
          [](const ECP& self, int i) -> const GaussianECP& {
            if (i < 0 || i >= self.getN()) throw py::index_error("primitive index out of range");
            return self.getGaussian(i);
          },
          py::arg("i"), py::return_value_policy::reference_internal)
      .def("getN", &ECP::getN)
      .def("noType1", &ECP::noType1)
      .def("evaluate", &ECP::evaluate, py::arg("r"), py::arg("l"))
      .def("getL", &ECP::getL)
      .def_readonly("N", &ECP::N)
      .def_readonly("L", &ECP::L);

  // --- ECPBasis ---
  py::class_<ECPBasis>(m, "ECPBasis")
      .def(py::init<>())
      .def("addECP", &ECPBasis::addECP, py::arg("U"), py::arg("atom"))
      .def(
          "getECP",
          [](ECPBasis& self, int i) -> ECP& {
            if (i < 0 || i >= self.getN()) throw py::index_error("ECP index out of range");
            return self.getECP(i);
          },
          py::arg("i"), py::return_value_policy::reference_internal)
      .def("getECPCore", &ECPBasis::getECPCore, py::arg("q"))
      .def(
          "getAtom",
          [](const ECPBasis& self, int i) {
            if (i < 0 || i >= self.getN()) throw py::index_error("ECP index out of range");
            return self.getAtom(i);
          },
          py::arg("i"))
      .def("getMaxL", &ECPBasis::getMaxL)
      .def("getN", &ECPBasis::getN)
      .def(
          "addECP_from_file",
          [](ECPBasis& self, int q, carray<double> coords, const std::string& filename) {
            check_size(coords.size(), 3, "coords");
            std::array<double, 3> arr = {coords.data()[0], coords.data()[1], coords.data()[2]};
            self.addECP_from_file(q, arr, filename);
          },
          py::arg("q"), py::arg("coords"), py::arg("filename"),
          "Read the ECP for atomic number q from a MOLPRO-format file and add it at coords. "
          "Raises RuntimeError if the file cannot be read or has no ECP for q.")
      .def_readwrite("core_electrons", &ECPBasis::core_electrons);

  // --- GaussianShell ---
  py::class_<GaussianShell>(m, "GaussianShell")
      .def(py::init([](carray<double> center, int l) {
             check_size(center.size(), 3, "center");
             std::array<double, 3> arr = {center.data()[0], center.data()[1], center.data()[2]};
             return GaussianShell(arr, l);
           }),
           py::arg("center"), py::arg("l"))
      .def("addPrim", &GaussianShell::addPrim, py::arg("exp"), py::arg("c"))
      .def("nprimitive", &GaussianShell::nprimitive)
      .def("ncartesian", &GaussianShell::ncartesian)
      .def("center",
           [](const GaussianShell& self) {
             double* c = self.center();
             return py::array_t<double>(3, c);
           })
      .def("exp", &GaussianShell::exp, py::arg("i"))
      .def("coef", &GaussianShell::coef, py::arg("i"))
      .def("am", &GaussianShell::am);

  // --- ShellPairData ---
  py::class_<ShellPairData>(m, "ShellPairData")
      .def(py::init<>())
      .def_readonly("LA", &ShellPairData::LA)
      .def_readonly("LB", &ShellPairData::LB)
      .def_readonly("maxLBasis", &ShellPairData::maxLBasis)
      .def_readonly("ncartA", &ShellPairData::ncartA)
      .def_readonly("ncartB", &ShellPairData::ncartB)
      .def_readonly("A2", &ShellPairData::A2)
      .def_readonly("Am", &ShellPairData::Am)
      .def_readonly("B2", &ShellPairData::B2)
      .def_readonly("Bm", &ShellPairData::Bm)
      .def_readonly("RAB2", &ShellPairData::RAB2)
      .def_readonly("RABm", &ShellPairData::RABm)
      .def_readonly("A_on_ecp", &ShellPairData::A_on_ecp)
      .def_readonly("B_on_ecp", &ShellPairData::B_on_ecp)
      .def_property_readonly("A",
                             [](const ShellPairData& self) {
                               return py::array_t<double>(3, self.A);
                             })
      .def_property_readonly("B", [](const ShellPairData& self) {
        return py::array_t<double>(3, self.B);
      });

  // --- TwoIndex<double> ---
  py::class_<TwoIndex<double>>(m, "TwoIndex")
      .def(py::init<>())
      .def(py::init<int, int>(), py::arg("dim1"), py::arg("dim2"))
      .def(py::init<int, int, double>(), py::arg("dim1"), py::arg("dim2"), py::arg("value"))
      .def("__getitem__",
           [](const TwoIndex<double>& self, std::pair<int, int> idx) {
             if (idx.first < 0 || idx.first >= self.dims[0] || idx.second < 0 ||
                 idx.second >= self.dims[1])
               throw py::index_error("TwoIndex index out of range");
             return self(idx.first, idx.second);
           })
      .def("__setitem__",
           [](TwoIndex<double>& self, std::pair<int, int> idx, double val) {
             if (idx.first < 0 || idx.first >= self.dims[0] || idx.second < 0 ||
                 idx.second >= self.dims[1])
               throw py::index_error("TwoIndex index out of range");
             self(idx.first, idx.second) = val;
           })
      .def_property_readonly("shape",
                             [](const TwoIndex<double>& self) {
                               return py::make_tuple(self.dims[0], self.dims[1]);
                             })
      .def("to_numpy", &to_numpy);

  // --- ECPIntegral ---
  py::class_<ECPIntegral>(m, "ECPIntegral")
      .def(py::init<int, int, int, double>(), py::arg("maxLB"), py::arg("maxLU"),
           py::arg("deriv") = 0, py::arg("thresh") = 1e-15,
           "Create a shell-pair integral engine for orbital angular momenta up to maxLB and ECP "
           "angular momenta up to maxLU, able to compute derivatives up to order deriv. thresh is "
           "the screening tolerance of the radial integrals (0 skips nothing). Raises ValueError "
           "if maxLB + deriv or maxLU exceeds MAX_L.")
      .def(
          "compute_shell_pair",
          [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
             const GaussianShell& shellB, int shiftA, int shiftB) {
            TwoIndex<double> values;
            self.compute_shell_pair(U, shellA, shellB, values, shiftA, shiftB);
            return values;
          },
          py::arg("U"), py::arg("shellA"), py::arg("shellB"), py::arg("shiftA") = 0,
          py::arg("shiftB") = 0)
      .def(
          "compute_shell_pair_derivative",
          [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
             const GaussianShell& shellB) {
            std::array<TwoIndex<double>, 9> results;
            self.compute_shell_pair_derivative(U, shellA, shellB, results);
            return std::vector<TwoIndex<double>>(results.begin(), results.end());
          },
          py::arg("U"), py::arg("shellA"), py::arg("shellB"),
          "First derivatives over the shell pair (A|U|B) on ECP centre C, in order "
          "[Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz]. Requires an ECPIntegral constructed with "
          "deriv >= 1.")
      .def(
          "compute_shell_pair_second_derivative",
          [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
             const GaussianShell& shellB) {
            std::array<TwoIndex<double>, 45> results;
            self.compute_shell_pair_second_derivative(U, shellA, shellB, results);
            return std::vector<TwoIndex<double>>(results.begin(), results.end());
          },
          py::arg("U"), py::arg("shellA"), py::arg("shellB"),
          "Second derivatives over the shell pair (A|U|B) on ECP centre C, in order "
          "[AA, AB, AC, BB, BC, CC], with components [xx, xy, xz, yy, yz, zz] for AA, BB, CC and "
          "[xx, xy, xz, yx, yy, yz, zx, zy, zz] for AB, AC, BC. Requires an ECPIntegral "
          "constructed with deriv >= 2.")
      .def(
          "compute_shell_pair_derivative_Bfield",
          [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
             const GaussianShell& shellB) {
            std::array<TwoIndex<double>, 3> results;
            self.compute_shell_pair_derivative_Bfield(U, shellA, shellB, results);
            return std::vector<TwoIndex<double>>(results.begin(), results.end());
          },
          py::arg("U"), py::arg("shellA"), py::arg("shellB"),
          "Derivatives [x, y, z] of the shell pair (A|U|B) of gauge-including atomic orbitals "
          "with respect to a uniform magnetic field. Each matrix is real, and the derivative is "
          "i/2 times it; the (B|A) matrices are minus the transposes. Requires an ECPIntegral "
          "constructed with deriv >= 1.")
      .def_readonly("skipped", &ECPIntegral::skipped)
      .def_readonly("zero", &ECPIntegral::zero)
      .def_readonly("nonzero", &ECPIntegral::nonzero);

  // --- ECPIntegrator ---
  py::class_<ECPIntegrator>(m, "ECPIntegrator")
      .def(py::init<>())
      .def(
          "set_gaussian_basis",
          [](ECPIntegrator& self, carray<double> coords, carray<double> exponents,
             carray<double> coefs, carray<int> ams, carray<int> shell_lengths) {
            const auto nshells = ams.size();
            check_size(shell_lengths.size(), nshells, "shell_lengths");
            check_size(coords.size(), 3 * nshells, "coords");
            const auto nprim = total_length(shell_lengths);
            check_size(exponents.size(), nprim, "exponents");
            check_size(coefs.size(), nprim, "coefs");
            self.set_gaussian_basis(static_cast<int>(nshells), coords.data(), exponents.data(),
                                    coefs.data(), ams.data(), shell_lengths.data());
          },
          py::arg("coords"), py::arg("exponents"), py::arg("coefs"), py::arg("ams"),
          py::arg("shell_lengths"))
      .def(
          "set_ecp_basis",
          [](ECPIntegrator& self, carray<double> coords, carray<double> exponents,
             carray<double> coefs, carray<int> ams, carray<int> ns, carray<int> shell_lengths) {
            const auto necps = shell_lengths.size();
            check_size(coords.size(), 3 * necps, "coords");
            const auto nprim = total_length(shell_lengths);
            check_size(exponents.size(), nprim, "exponents");
            check_size(coefs.size(), nprim, "coefs");
            check_size(ams.size(), nprim, "ams");
            check_size(ns.size(), nprim, "ns");
            self.set_ecp_basis(static_cast<int>(necps), coords.data(), exponents.data(),
                               coefs.data(), ams.data(), ns.data(), shell_lengths.data());
          },
          py::arg("coords"), py::arg("exponents"), py::arg("coefs"), py::arg("ams"), py::arg("ns"),
          py::arg("shell_lengths"))
      .def(
          "set_ecp_basis_from_library",
          [](ECPIntegrator& self, carray<double> coords, carray<int> charges,
             const std::vector<std::string>& names, const std::string& share_dir) {
            const auto necps = charges.size();
            check_size(coords.size(), 3 * necps, "coords");
            check_size(static_cast<py::ssize_t>(names.size()), necps, "names");
            self.set_ecp_basis_from_library(static_cast<int>(necps), coords.data(),
                                            charges.data(), names, share_dir);
          },
          py::arg("coords"), py::arg("charges"), py::arg("names"), py::arg("share_dir"))
      .def(
          "update_gaussian_basis_coords",
          [](ECPIntegrator& self, carray<double> coords) {
            const auto nshells = static_cast<py::ssize_t>(self.shells.size());
            check_size(coords.size(), 3 * nshells, "coords");
            self.update_gaussian_basis_coords(static_cast<int>(nshells), coords.data());
          },
          py::arg("coords"))
      .def(
          "update_ecp_basis_coords",
          [](ECPIntegrator& self, carray<double> coords) {
            const py::ssize_t necps = self.ecps.getN();
            check_size(coords.size(), 3 * necps, "coords");
            self.update_ecp_basis_coords(static_cast<int>(necps), coords.data());
          },
          py::arg("coords"))
      .def(
          "init",
          [](ECPIntegrator& self, int deriv, double tol, double radial_tol) {
            if (!self.basis_is_set || !self.ecp_is_set)
              throw std::runtime_error(
                  "the Gaussian and ECP bases must both be set before calling init");
            self.init(deriv, tol, radial_tol);
          },
          py::arg("deriv") = 0, py::arg("tol") = 1e-12, py::arg("radial_tol") = 1e-15)
      .def("compute_integrals",
           [](ECPIntegrator& self) {
             require_init(self);
             py::gil_scoped_release release;
             self.compute_integrals();
           })
      .def("compute_first_derivs",
           [](ECPIntegrator& self) {
             require_deriv(self, 1, "compute_first_derivs");
             py::gil_scoped_release release;
             self.first_derivs.clear();
             self.compute_first_derivs();
           })
      .def("compute_second_derivs",
           [](ECPIntegrator& self) {
             require_deriv(self, 2, "compute_second_derivs");
             py::gil_scoped_release release;
             self.second_derivs.clear();
             self.compute_second_derivs();
           })
      .def("compute_bfield_derivs",
           [](ECPIntegrator& self) {
             require_deriv(self, 1, "compute_bfield_derivs");
             py::gil_scoped_release release;
             self.compute_bfield_derivs();
           })
      .def(
          "compute_gradient",
          [](ECPIntegrator& self, carray<double> density) {
            require_deriv(self, 1, "compute_gradient");
            check_density(self, density);
            std::vector<double> g;
            {
              py::gil_scoped_release release;
              g = self.compute_gradient(density.data());
            }
            return py::array_t<double>(static_cast<py::ssize_t>(g.size()), g.data());
          },
          py::arg("density"))
      .def(
          "compute_hessian",
          [](ECPIntegrator& self, carray<double> density) {
            require_deriv(self, 2, "compute_hessian");
            check_density(self, density);
            std::vector<double> h;
            {
              py::gil_scoped_release release;
              h = self.compute_hessian(density.data());
            }
            const py::ssize_t n = 3 * self.natoms;
            return py::array_t<double>({n, n}, h.data());
          },
          py::arg("density"))
      .def("get_integrals", [](const ECPIntegrator& self) { return to_numpy(self.integrals); })
      .def("get_first_derivs",
           [](const ECPIntegrator& self) { return to_numpy_list(self.first_derivs); })
      .def("get_second_derivs",
           [](const ECPIntegrator& self) { return to_numpy_list(self.second_derivs); })
      .def("get_bfield_derivs",
           [](const ECPIntegrator& self) { return to_numpy_list(self.bfield_derivs); })
      .def_readwrite("tolerance", &ECPIntegrator::tolerance)
      .def_readonly("ncart", &ECPIntegrator::ncart)
      .def_readonly("maxLB", &ECPIntegrator::maxLB)
      // natoms and deriv are only set by init
      .def_property_readonly("natoms",
                             [](const ECPIntegrator& self) { return self.ecpint ? self.natoms : 0; })
      .def_property_readonly("deriv",
                             [](const ECPIntegrator& self) { return self.ecpint ? self.deriv : 0; });
}
