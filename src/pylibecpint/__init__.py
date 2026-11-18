"""Python bindings for libecpint — ECP integral evaluation."""

from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files
from pathlib import Path

import numpy as np

from pylibecpint._core import (
    ECP,
    MAX_L,
    ECPBasis,
    ECPIntegral,
    ECPIntegrator as _ECPIntegrator,
    GaussianECP,
    GaussianShell,
    ShellPairData,
    TwoIndex,
    __libecpint_version__,
)

try:
    __version__ = version("pylibecpint")
except PackageNotFoundError:  # pragma: no cover - running from an uninstalled source tree
    __version__ = "unknown"

__all__ = [
    "ECPIntegrator",
    "ECPIntegral",
    "GaussianShell",
    "ECP",
    "GaussianECP",
    "ECPBasis",
    "ShellPairData",
    "TwoIndex",
    "get_data_dir",
    "available_ecps",
    "MAX_L",
    "__version__",
    "__libecpint_version__",
]


def get_data_dir() -> str:
    """Return the path to the bundled ECP library.

    This is the libecpint ``share`` directory: the ECP definitions are the MOLPRO-format files
    ``raw/<name>.ecp`` inside it. It is the default ``share_dir`` of
    :meth:`ECPIntegrator.set_ecp_basis_from_library`.
    """
    installed = Path(str(files("pylibecpint"))) / "data"
    if (installed / "raw").is_dir():
        return str(installed)
    # Editable install: resolve from the submodule
    repo_root = Path(__file__).resolve().parent.parent.parent
    submodule = repo_root / "libecpint" / "share" / "libecpint"
    if (submodule / "raw").is_dir():
        return str(submodule)
    raise FileNotFoundError(
        f"Cannot find the bundled ECP library (looked for {installed / 'raw'} "
        f"and {submodule / 'raw'})"
    )


def available_ecps() -> list[str]:
    """Return the names of the ECP families in the bundled library, e.g. ``"ecp10mdf"``.

    These are the names accepted by :meth:`ECPIntegrator.set_ecp_basis_from_library`.
    """
    return sorted(p.stem for p in (Path(get_data_dir()) / "raw").glob("*.ecp"))


def _f64(a) -> np.ndarray:
    return np.ascontiguousarray(a, dtype=np.float64)


def _i32(a) -> np.ndarray:
    return np.ascontiguousarray(a, dtype=np.int32)


class ECPIntegrator:
    """High-level interface for computing ECP integrals and their derivatives.

    Wraps the C++ ECPIntegrator with numpy array inputs/outputs
    and automatic resolution of the bundled ECP data directory.

    Usage::

        integrator = ECPIntegrator()
        integrator.set_gaussian_basis(coords, exponents, coefs, ams, shell_lengths)
        integrator.set_ecp_basis(coords, exponents, coefs, ams, ns, shell_lengths)
        integrator.init()
        integrals = integrator.compute_integrals()

    Results are in canonical Cartesian order, with shells in the order given to
    :meth:`set_gaussian_basis`. Derivatives are with respect to the coordinates of each distinct
    atom, in the order in which the atoms first appear in :meth:`set_gaussian_basis`.
    """

    def __init__(self) -> None:
        self._impl = _ECPIntegrator()

    @property
    def tolerance(self) -> float:
        """Screening threshold for shells and shell pairs (set by :meth:`init`)."""
        return self._impl.tolerance

    @tolerance.setter
    def tolerance(self, value: float) -> None:
        self._impl.tolerance = value

    @property
    def ncart(self) -> int:
        """Number of Cartesian basis functions."""
        return self._impl.ncart

    @property
    def natoms(self) -> int:
        """Number of distinct atoms (0 before :meth:`init`)."""
        return self._impl.natoms

    @property
    def maxLB(self) -> int:
        """Maximum angular momentum in the Gaussian basis."""
        return self._impl.maxLB

    @property
    def deriv(self) -> int:
        """Maximum derivative order the integrator was initialised for (0 before :meth:`init`)."""
        return self._impl.deriv

    def set_gaussian_basis(
        self,
        coords: np.ndarray,
        exponents: np.ndarray,
        coefs: np.ndarray,
        ams: np.ndarray,
        shell_lengths: np.ndarray,
    ) -> None:
        """Set the Gaussian orbital basis from arrays.

        Args:
            coords: Cartesian coordinates in bohr, shape (3*nshells,).
            exponents: Primitive exponents, shape (sum(shell_lengths),).
            coefs: Contraction coefficients, same shape as exponents.
            ams: Angular momentum per shell, shape (nshells,).
            shell_lengths: Number of primitives per shell, shape (nshells,).
        """
        self._impl.set_gaussian_basis(
            _f64(coords), _f64(exponents), _f64(coefs), _i32(ams), _i32(shell_lengths)
        )

    def set_ecp_basis(
        self,
        coords: np.ndarray,
        exponents: np.ndarray,
        coefs: np.ndarray,
        ams: np.ndarray,
        ns: np.ndarray,
        shell_lengths: np.ndarray,
    ) -> None:
        """Set the ECP basis from raw arrays.

        Args:
            coords: ECP coordinates in bohr, shape (3*necps,).
            exponents: Primitive exponents.
            coefs: Primitive coefficients.
            ams: Angular momentum per primitive (up to :data:`MAX_L`).
            ns: Power of r per primitive, with n=2 meaning r^0.
            shell_lengths: Number of primitives per ECP, shape (necps,).

        Raises:
            ValueError: if an angular momentum or power is out of range.
        """
        self._impl.set_ecp_basis(
            _f64(coords), _f64(exponents), _f64(coefs), _i32(ams), _i32(ns), _i32(shell_lengths)
        )

    def set_ecp_basis_from_library(
        self,
        coords: np.ndarray,
        charges: np.ndarray,
        names: list[str],
        share_dir: str | None = None,
    ) -> None:
        """Load ECPs from the bundled library of MOLPRO-format ``.ecp`` files.

        Args:
            coords: ECP coordinates in bohr, shape (3*necps,).
            charges: Atomic numbers, shape (necps,).
            names: ECP family names (e.g. "ecp10mdf"), one per ECP; see :func:`available_ecps`.
            share_dir: Directory containing ``raw/<name>.ecp``. Defaults to the bundled library.

        Raises:
            RuntimeError: if a file cannot be read or has no ECP for the given atom.
        """
        if share_dir is None:
            share_dir = get_data_dir()
        self._impl.set_ecp_basis_from_library(_f64(coords), _i32(charges), list(names), share_dir)

    def update_gaussian_basis_coords(self, coords: np.ndarray) -> None:
        """Update Gaussian shell positions. Order must match set_gaussian_basis."""
        self._impl.update_gaussian_basis_coords(_f64(coords))

    def update_ecp_basis_coords(self, coords: np.ndarray) -> None:
        """Update ECP positions. Order must match set_ecp_basis."""
        self._impl.update_ecp_basis_coords(_f64(coords))

    def init(self, deriv: int = 0, tol: float = 1e-12, radial_tol: float = 1e-15) -> None:
        """Initialise the integrator. Must be called after setting both bases.

        Args:
            deriv: Maximum derivative order to be computed (0, 1 or 2). Magnetic field
                derivatives need 1.
            tol: Screening tolerance for shells and shell pairs.
            radial_tol: Screening tolerance of the radial integrals (0 skips nothing).

        Raises:
            ValueError: if the basis angular momentum plus ``deriv`` exceeds :data:`MAX_L`.
        """
        self._impl.init(deriv, tol, radial_tol)

    def compute_integrals(self) -> np.ndarray:
        """Compute and return the ECP integrals as an (ncart, ncart) array."""
        self._impl.compute_integrals()
        return self._impl.get_integrals()

    def compute_first_derivs(self) -> list[np.ndarray]:
        """Compute and return the first derivatives as a list of (ncart, ncart) arrays.

        Order: [Ax, Ay, Az, Bx, By, Bz, ...] matching atom order from set_gaussian_basis.
        All 3*natoms matrices are held in memory, which grows as natoms * ncart^2 (e.g. ~7 GB
        for 77 atoms with 1900 basis functions); use :meth:`compute_gradient` when only the
        contraction with a density matrix is needed. Requires ``init(deriv >= 1)``.
        """
        self._impl.compute_first_derivs()
        return self._impl.get_first_derivs()

    def compute_second_derivs(self) -> list[np.ndarray]:
        """Compute and return the second derivatives as a list of (ncart, ncart) arrays.

        The order is {AA, AB, AC, ..., BB, BC, ..., CC, ...} over atom pairs, with components
        {xx, xy, xz, yy, yz, zz} for each diagonal pair (AA, BB, ...) and
        {xx, xy, xz, yx, yy, yz, zx, zy, zz} for each off-diagonal pair (AB, AC, ...), for a total
        of 3*natoms*(3*natoms+1)/2 matrices. Memory grows as natoms^2 * ncart^2, so use
        :meth:`compute_hessian` beyond small systems. Requires ``init(deriv >= 2)``.
        """
        self._impl.compute_second_derivs()
        return self._impl.get_second_derivs()

    def compute_gradient(self, density: np.ndarray) -> np.ndarray:
        """Contract the first derivatives of the integrals V with a density matrix.

        ``gradient[3*X + i] = sum_kl density[k, l] * dV[k, l]/dX_i`` for each atom X and
        i = x, y, z. The derivative matrices are never stored, so this is the recommended way to
        get derivatives for large systems. Requires ``init(deriv >= 1)``.

        Args:
            density: (ncart, ncart) matrix in the order of the integrals (need not be symmetric).

        Returns:
            Array of shape (3*natoms,).
        """
        return self._impl.compute_gradient(_f64(density))

    def compute_hessian(self, density: np.ndarray) -> np.ndarray:
        """Contract the second derivatives of the integrals V with a density matrix.

        ``hessian[3*X + i, 3*Y + j] = sum_kl density[k, l] * d2V[k, l]/dX_i dY_j``, as a full
        symmetric matrix. The derivative matrices are never stored. Requires
        ``init(deriv >= 2)``.

        Args:
            density: (ncart, ncart) matrix in the order of the integrals (need not be symmetric).

        Returns:
            Array of shape (3*natoms, 3*natoms).
        """
        return self._impl.compute_hessian(_f64(density))

    def compute_bfield_derivs(self) -> list[np.ndarray]:
        """Compute the derivatives with respect to a uniform magnetic field, for GIAOs.

        Returns three (ncart, ncart) matrices, for the field components x, y, z. Each is real
        and antisymmetric, and the derivative of the integrals is i/2 times it.
        Requires ``init(deriv >= 1)``.
        """
        self._impl.compute_bfield_derivs()
        return self._impl.get_bfield_derivs()
