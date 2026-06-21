"""Python bindings for libecpint — ECP integral evaluation."""

from pathlib import Path

from importlib.resources import files

import numpy as np

from pylibecpint._core import (
    ECP,
    ECPBasis,
    ECPIntegral,
    ECPIntegrator as _ECPIntegrator,
    GaussianECP,
    GaussianShell,
    ShellPairData,
    TwoIndex,
)

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
]


def get_data_dir() -> str:
    """Return the path to the bundled ECP XML data directory."""
    installed = Path(str(files("pylibecpint"))) / "data" / "xml"
    if installed.is_dir():
        return str(installed)
    # Editable install: resolve from submodule
    repo_root = Path(__file__).resolve().parent.parent.parent
    submodule = repo_root / "libecpint" / "share" / "libecpint" / "xml"
    if submodule.is_dir():
        return str(submodule)
    raise FileNotFoundError("Cannot find ECP XML data directory")


class ECPIntegrator:
    """High-level interface for computing ECP integrals.

    Wraps the C++ ECPIntegrator with numpy array inputs/outputs
    and automatic resolution of the bundled ECP data directory.

    Usage::

        integrator = ECPIntegrator()
        integrator.set_gaussian_basis(coords, exponents, coefs, ams, shell_lengths)
        integrator.set_ecp_basis(coords, exponents, coefs, ams, ns, shell_lengths)
        integrator.init()
        integrals = integrator.compute_integrals()
    """

    def __init__(self) -> None:
        self._impl = _ECPIntegrator()

    @property
    def tolerance(self) -> float:
        return self._impl.tolerance

    @tolerance.setter
    def tolerance(self, value: float) -> None:
        self._impl.tolerance = value

    @property
    def shell_tolerance(self) -> float:
        return self._impl.shell_tolerance

    @shell_tolerance.setter
    def shell_tolerance(self, value: float) -> None:
        self._impl.shell_tolerance = value

    @property
    def ncart(self) -> int:
        return self._impl.ncart

    @property
    def natoms(self) -> int:
        return self._impl.natoms

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
            np.ascontiguousarray(coords, dtype=np.float64),
            np.ascontiguousarray(exponents, dtype=np.float64),
            np.ascontiguousarray(coefs, dtype=np.float64),
            np.ascontiguousarray(ams, dtype=np.int32),
            np.ascontiguousarray(shell_lengths, dtype=np.int32),
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
            ams: Angular momentum per primitive.
            ns: Power of r per primitive.
            shell_lengths: Number of primitives per ECP, shape (necps,).
        """
        self._impl.set_ecp_basis(
            np.ascontiguousarray(coords, dtype=np.float64),
            np.ascontiguousarray(exponents, dtype=np.float64),
            np.ascontiguousarray(coefs, dtype=np.float64),
            np.ascontiguousarray(ams, dtype=np.int32),
            np.ascontiguousarray(ns, dtype=np.int32),
            np.ascontiguousarray(shell_lengths, dtype=np.int32),
        )

    def set_ecp_basis_from_library(
        self,
        coords: np.ndarray,
        charges: np.ndarray,
        names: list[str],
        share_dir: str | None = None,
    ) -> None:
        """Load ECPs from the bundled XML library.

        Args:
            coords: ECP coordinates in bohr, shape (3*necps,).
            charges: Atomic numbers, shape (necps,).
            names: ECP names (e.g. "ecp10mdf"), one per ECP.
            share_dir: Path to XML data directory. Defaults to bundled data.
        """
        if share_dir is None:
            share_dir = get_data_dir()
        self._impl.set_ecp_basis_from_library(
            np.ascontiguousarray(coords, dtype=np.float64),
            np.ascontiguousarray(charges, dtype=np.int32),
            names,
            share_dir,
        )

    def update_gaussian_basis_coords(self, coords: np.ndarray) -> None:
        """Update Gaussian shell positions. Order must match set_gaussian_basis."""
        self._impl.update_gaussian_basis_coords(
            np.ascontiguousarray(coords, dtype=np.float64),
        )

    def update_ecp_basis_coords(self, coords: np.ndarray) -> None:
        """Update ECP positions. Order must match set_ecp_basis."""
        self._impl.update_ecp_basis_coords(
            np.ascontiguousarray(coords, dtype=np.float64),
        )

    def init(
        self,
        deriv: int = 0,
        tol: float = 1e-12,
        radial_tol: float = 1e-15,
        small_grid: int = 256,
        big_grid: int = 1024,
    ) -> None:
        """Initialize the integrator. Must be called after setting bases."""
        self._impl.init(deriv, tol, radial_tol, small_grid, big_grid)

    def compute_integrals(self) -> np.ndarray:
        """Compute and return ECP integrals as a 2D numpy array."""
        self._impl.compute_integrals()
        return self._impl.get_integrals()

    def compute_first_derivs(self) -> list[np.ndarray]:
        """Compute and return first derivatives as a list of 2D numpy arrays.

        Order: [Ax, Ay, Az, Bx, By, Bz, ...] matching atom order from set_gaussian_basis.
        """
        self._impl.compute_first_derivs()
        return self._impl.get_first_derivs()

    def compute_second_derivs(self) -> list[np.ndarray]:
        """Compute and return second derivatives as a list of 2D numpy arrays."""
        self._impl.compute_second_derivs()
        return self._impl.get_second_derivs()
