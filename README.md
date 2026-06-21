# pylibecpint

[![Documentation](https://readthedocs.org/projects/pylibecpint/badge/?version=latest)](https://pylibecpint.readthedocs.io)

Python bindings for [libecpint](https://github.com/robashaw/libecpint), a C++ library for the efficient evaluation of integrals over effective core potentials (ECPs).

**[Full documentation](https://pylibecpint.readthedocs.io)** — installation, configuration, working examples, and API reference.

## Installation

```bash
pip install pylibecpint
```

### From source

Requires a C++ compiler and CMake >= 3.15.

```bash
git clone --recurse-submodules https://github.com/robashaw/pylibecpint.git
cd pylibecpint
pip install .
```

To build with a higher maximum angular momentum (default is 5) or max unrol (default is 1):

```bash
CMAKE_ARGS="-DLIBECPINT_MAX_L=7 -DLIBECPINT_MAX_UNROL=2" pip install .
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
integrator.init()
integrals = integrator.compute_integrals()  # returns 2D numpy array
```

## Features

- Full access to the libecpint high-level API (`ECPIntegrator`)
- Lower-level shell-pair integrals via `ECPIntegral`
- First and second derivative support
- Bundled ECP definitions (Stuttgart-Dresden, LANL2DZ, etc.)
- All inputs/outputs as numpy arrays

## API

### High-level: `ECPIntegrator`

The main interface. Set up bases, initialize, compute:

- `set_gaussian_basis(coords, exponents, coefs, ams, shell_lengths)` — define the orbital basis
- `set_ecp_basis(coords, exponents, coefs, ams, ns, shell_lengths)` — define ECPs from raw data
- `set_ecp_basis_from_library(coords, charges, names)` — load ECPs from bundled XML files
- `init(deriv=0)` — initialize the integrator
- `compute_integrals()` — returns `np.ndarray` of shape `(ncart, ncart)`
- `compute_first_derivs()` — returns list of `np.ndarray`, one per atomic coordinate
- `compute_second_derivs()` — returns list of `np.ndarray`

### Low-level: `ECPIntegral`

For per-shell-pair calculations:

- `ECPIntegral(maxLB, maxLU, deriv=0)`
- `compute_shell_pair(U, shellA, shellB)` — returns `TwoIndex`

### Data structures

- `GaussianShell(center, l)` — contracted Gaussian shell
- `ECP(center)` — effective core potential
- `ECPBasis()` — container of ECPs

All coordinates must be in **bohr**. Results are in **canonical Cartesian order**.

## License

MIT — see [LICENSE](LICENSE).
