# pylibecpint

[![Documentation](https://readthedocs.org/projects/pylibecpint/badge/?version=latest)](https://pylibecpint.readthedocs.io)

Python bindings for [libecpint](https://github.com/robashaw/libecpint), a C++ library for the efficient evaluation of integrals over effective core potentials (ECPs).

pylibecpint 1.1.0 bundles libecpint 1.1.0; see libecpint's [What's new in v1.1](https://github.com/robashaw/libecpint#whats-new-in-v11) for the changes in the underlying library.

**[Full documentation](https://pylibecpint.readthedocs.io)** — installation, configuration, working examples, and API reference.

## Installation

```bash
pip install pylibecpint
```

### From source

Requires a C++17 compiler and CMake >= 3.15.

```bash
git clone --recurse-submodules https://github.com/robashaw/pylibecpint.git
cd pylibecpint
pip install .
```

Angular momentum up to 8 is always supported (for the basis, this includes the derivative order), so no build options are needed. Optionally, more of the integral code can be unrolled (default 1), or the per-thread radial cache enlarged (default 4 MB):

```bash
CMAKE_ARGS="-DLIBECPINT_MAX_UNROL=2 -DLIBECPINT_RADIAL_CACHE_MB=8" pip install .
```

## Quick start

```python
import numpy as np
from pylibecpint import ECPIntegrator

integrator = ECPIntegrator()

# Set up Gaussian basis (coords in bohr)
integrator.set_gaussian_basis(coords, exponents, coefs, ams, shell_lengths)

# Load ECPs from the built-in library
integrator.set_ecp_basis_from_library(ecp_coords, charges, ["ecp10mdf"])

# Compute integrals
integrator.init(deriv=1)
integrals = integrator.compute_integrals()  # returns 2D numpy array

# Gradient contribution for a density matrix, without storing the derivative matrices
gradient = integrator.compute_gradient(density)  # shape (3 * natoms,)
```

## Features

- Full access to the libecpint high-level API (`ECPIntegrator`)
- Lower-level shell-pair integrals via `ECPIntegral`
- First and second derivatives, either as full matrices or contracted directly with a density matrix (`compute_gradient`, `compute_hessian`) to save memory
- Derivatives with respect to a magnetic field, for gauge-including atomic orbitals
- Angular momentum up to 8 (l functions)
- 91 bundled ECP families in MOLPRO format (Stuttgart/Cologne, Karlsruhe def2/dhf, LANL2DZ, CRENBL/CRENBS, SBKJC, ...); list them with `pylibecpint.available_ecps()`
- All inputs/outputs as numpy arrays

## API

### High-level: `ECPIntegrator`

The main interface. Set up bases, initialize, compute:

- `set_gaussian_basis(coords, exponents, coefs, ams, shell_lengths)` — define the orbital basis
- `set_ecp_basis(coords, exponents, coefs, ams, ns, shell_lengths)` — define ECPs from raw data
- `set_ecp_basis_from_library(coords, charges, names)` — load ECPs from the bundled library
- `init(deriv=0)` — initialize the integrator
- `compute_integrals()` — returns `np.ndarray` of shape `(ncart, ncart)`
- `compute_gradient(density)` / `compute_hessian(density)` — derivatives contracted with a density matrix, shapes `(3*natoms,)` and `(3*natoms, 3*natoms)`
- `compute_first_derivs()` — returns list of `np.ndarray`, one per atomic coordinate
- `compute_second_derivs()` — returns list of `np.ndarray`
- `compute_bfield_derivs()` — returns three antisymmetric `np.ndarray`, for the x, y, z field components

### Low-level: `ECPIntegral`

For per-shell-pair calculations:

- `ECPIntegral(maxLB, maxLU, deriv=0)`
- `compute_shell_pair(U, shellA, shellB)` — returns `TwoIndex`
- `compute_shell_pair_derivative`, `compute_shell_pair_second_derivative`, `compute_shell_pair_derivative_Bfield` — return lists of `TwoIndex`

### Data structures

- `GaussianShell(center, l)` — contracted Gaussian shell
- `ECP(center)` — effective core potential
- `ECPBasis()` — container of ECPs

All coordinates must be in **bohr**. Results are in **canonical Cartesian order**.

## License

MIT — see [LICENSE](LICENSE).
