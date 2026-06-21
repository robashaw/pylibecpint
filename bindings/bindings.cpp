#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <vector>

#include "api.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"
#include "multiarr.hpp"

namespace py = pybind11;
using namespace libecpint;

PYBIND11_MODULE(_core, m) {
  m.doc() = "Python bindings for libecpint";

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
      .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> center) {
        if (center.size() != 3) throw std::runtime_error("center must have 3 elements");
        return ECP(center.data());
      }))
      .def("addPrimitive", &ECP::addPrimitive, py::arg("n"), py::arg("l"), py::arg("a"),
           py::arg("d"), py::arg("needSort") = true)
      .def("center",
           [](const ECP& self) {
             const double* c = self.center();
             return py::array_t<double>(3, c);
           })
      .def("setPos", &ECP::setPos)
      .def("sort", &ECP::sort)
      .def("getGaussian",
           [](const ECP& self, int i) -> const GaussianECP& { return self.getGaussian(i); },
           py::return_value_policy::reference_internal)
      .def("getN", &ECP::getN)
      .def("noType1", &ECP::noType1)
      .def("evaluate", &ECP::evaluate)
      .def("getL", &ECP::getL)
      .def_readonly("N", &ECP::N)
      .def_readonly("L", &ECP::L);

  // --- ECPBasis ---
  py::class_<ECPBasis>(m, "ECPBasis")
      .def(py::init<>())
      .def("addECP", &ECPBasis::addECP)
      .def(
          "getECP", [](ECPBasis& self, int i) -> ECP& { return self.getECP(i); },
          py::return_value_policy::reference_internal)
      .def("getECPCore", &ECPBasis::getECPCore)
      .def("getAtom", &ECPBasis::getAtom)
      .def("getMaxL", &ECPBasis::getMaxL)
      .def("getN", &ECPBasis::getN)
#ifdef HAS_PUGIXML
      .def("addECP_from_file",
           [](ECPBasis& self, int q,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords,
              const std::string& filename) {
             if (coords.size() != 3)
               throw std::runtime_error("coords must have 3 elements");
             std::array<double, 3> arr = {coords.data()[0], coords.data()[1], coords.data()[2]};
             self.addECP_from_file(q, arr, filename);
           })
#endif
      .def_readwrite("core_electrons", &ECPBasis::core_electrons);

  // --- GaussianShell ---
  py::class_<GaussianShell>(m, "GaussianShell")
      .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> center,
                       int l) {
             if (center.size() != 3)
               throw std::runtime_error("center must have 3 elements");
             std::array<double, 3> arr = {center.data()[0], center.data()[1], center.data()[2]};
             return GaussianShell(arr, l);
           }),
           py::arg("center"), py::arg("l"))
      .def("addPrim", &GaussianShell::addPrim)
      .def("nprimitive", &GaussianShell::nprimitive)
      .def("ncartesian", &GaussianShell::ncartesian)
      .def("center",
           [](const GaussianShell& self) {
             double* c = self.center();
             return py::array_t<double>(3, c);
           })
      .def("exp", &GaussianShell::exp)
      .def("coef", &GaussianShell::coef)
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
      .def(py::init<int, int>())
      .def(py::init<int, int, double>())
      .def("__getitem__",
           [](const TwoIndex<double>& self, std::pair<int, int> idx) {
             return self(idx.first, idx.second);
           })
      .def("__setitem__",
           [](TwoIndex<double>& self, std::pair<int, int> idx, double val) {
             self(idx.first, idx.second) = val;
           })
      .def_property_readonly("shape",
                             [](const TwoIndex<double>& self) {
                               return py::make_tuple(self.dims[0], self.dims[1]);
                             })
      .def("to_numpy", [](const TwoIndex<double>& self) {
        return py::array_t<double>({self.dims[0], self.dims[1]}, self.data.data());
      });

  // --- ECPIntegral ---
  py::class_<ECPIntegral>(m, "ECPIntegral")
      .def(py::init<int, int, int, double, unsigned, unsigned>(), py::arg("maxLB"),
           py::arg("maxLU"), py::arg("deriv") = 0, py::arg("thresh") = 1e-15,
           py::arg("smallGrid") = 256, py::arg("bigGrid") = 1024)
      .def("compute_shell_pair",
           [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
              const GaussianShell& shellB, int shiftA, int shiftB) {
             TwoIndex<double> values;
             self.compute_shell_pair(U, shellA, shellB, values, shiftA, shiftB);
             return values;
           },
           py::arg("U"), py::arg("shellA"), py::arg("shellB"), py::arg("shiftA") = 0,
           py::arg("shiftB") = 0)
      .def("compute_shell_pair_derivative",
           [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
              const GaussianShell& shellB) {
             std::array<TwoIndex<double>, 9> results;
             self.compute_shell_pair_derivative(U, shellA, shellB, results);
             std::vector<TwoIndex<double>> out(results.begin(), results.end());
             return out;
           })
      .def("compute_shell_pair_second_derivative",
           [](ECPIntegral& self, const ECP& U, const GaussianShell& shellA,
              const GaussianShell& shellB) {
             std::array<TwoIndex<double>, 45> results;
             self.compute_shell_pair_second_derivative(U, shellA, shellB, results);
             std::vector<TwoIndex<double>> out(results.begin(), results.end());
             return out;
           })
      .def_readonly("skipped", &ECPIntegral::skipped)
      .def_readonly("zero", &ECPIntegral::zero)
      .def_readonly("nonzero", &ECPIntegral::nonzero);

  // --- ECPIntegrator ---
  py::class_<ECPIntegrator>(m, "ECPIntegrator")
      .def(py::init<>())
      .def("set_gaussian_basis",
           [](ECPIntegrator& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords,
              py::array_t<double, py::array::c_style | py::array::forcecast> exponents,
              py::array_t<double, py::array::c_style | py::array::forcecast> coefs,
              py::array_t<int, py::array::c_style | py::array::forcecast> ams,
              py::array_t<int, py::array::c_style | py::array::forcecast> shell_lengths) {
             int nshells = static_cast<int>(ams.size());
             self.set_gaussian_basis(nshells, coords.data(), exponents.data(), coefs.data(),
                                    ams.data(), shell_lengths.data());
           })
      .def("set_ecp_basis",
           [](ECPIntegrator& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords,
              py::array_t<double, py::array::c_style | py::array::forcecast> exponents,
              py::array_t<double, py::array::c_style | py::array::forcecast> coefs,
              py::array_t<int, py::array::c_style | py::array::forcecast> ams,
              py::array_t<int, py::array::c_style | py::array::forcecast> ns,
              py::array_t<int, py::array::c_style | py::array::forcecast> shell_lengths) {
             int necps = static_cast<int>(shell_lengths.size());
             self.set_ecp_basis(necps, coords.data(), exponents.data(), coefs.data(), ams.data(),
                                ns.data(), shell_lengths.data());
           })
#ifdef HAS_PUGIXML
      .def("set_ecp_basis_from_library",
           [](ECPIntegrator& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords,
              py::array_t<int, py::array::c_style | py::array::forcecast> charges,
              std::vector<std::string> names, const std::string& share_dir) {
             int necps = static_cast<int>(charges.size());
             self.set_ecp_basis_from_library(necps, coords.data(), charges.data(), names,
                                            share_dir);
           })
#endif
      .def("update_gaussian_basis_coords",
           [](ECPIntegrator& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords) {
             int nshells = static_cast<int>(coords.size()) / 3;
             self.update_gaussian_basis_coords(nshells, coords.data());
           })
      .def("update_ecp_basis_coords",
           [](ECPIntegrator& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> coords) {
             int necps = static_cast<int>(coords.size()) / 3;
             self.update_ecp_basis_coords(necps, coords.data());
           })
      .def("init", &ECPIntegrator::init, py::arg("deriv") = 0, py::arg("tol") = 1e-12,
           py::arg("radial_tol") = 1e-15, py::arg("smallGrid") = 256, py::arg("bigGrid") = 1024)
      .def("compute_integrals", &ECPIntegrator::compute_integrals)
      .def("compute_first_derivs", &ECPIntegrator::compute_first_derivs)
      .def("compute_second_derivs", &ECPIntegrator::compute_second_derivs)
      .def("get_integrals",
           [](ECPIntegrator& self) {
             auto ptr = self.get_integrals();
             return py::array_t<double>({self.ncart, self.ncart}, ptr->data());
           })
      .def("get_first_derivs",
           [](ECPIntegrator& self) {
             auto ptrs = self.get_first_derivs();
             py::list result;
             for (auto& ptr : ptrs)
               result.append(py::array_t<double>({self.ncart, self.ncart}, ptr->data()));
             return result;
           })
      .def("get_second_derivs",
           [](ECPIntegrator& self) {
             auto ptrs = self.get_second_derivs();
             py::list result;
             for (auto& ptr : ptrs)
               result.append(py::array_t<double>({self.ncart, self.ncart}, ptr->data()));
             return result;
           })
      .def_readwrite("tolerance", &ECPIntegrator::tolerance)
      .def_readwrite("shell_tolerance", &ECPIntegrator::shell_tolerance)
      .def_readonly("ncart", &ECPIntegrator::ncart)
      .def_readonly("natoms", &ECPIntegrator::natoms)
      .def_readonly("maxLB", &ECPIntegrator::maxLB)
      .def_readonly("deriv", &ECPIntegrator::deriv);
}
