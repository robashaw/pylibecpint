.. _introduction:

Introduction
============

pylibecpint provides Python bindings for
`libecpint <https://github.com/robashaw/libecpint>`_, enabling the
evaluation of integrals over effective core potentials (ECPs) directly from
Python. All inputs and outputs use NumPy arrays.

pylibecpint 1.1.0 bundles libecpint 1.1.0, which is available at runtime as
``pylibecpint.__libecpint_version__``. See libecpint's
`What's new in v1.1 <https://github.com/robashaw/libecpint#whats-new-in-v11>`_ for the
changes in the underlying library.

Installation
------------

From PyPI
^^^^^^^^^

.. code-block:: bash

   pip install pylibecpint

From source
^^^^^^^^^^^

Building from source requires a C++17 compiler and CMake >= 3.15.

.. code-block:: bash

   git clone --recurse-submodules https://github.com/robashaw/pylibecpint.git
   cd pylibecpint
   pip install .

Configuration
-------------

Angular momentum
^^^^^^^^^^^^^^^^

Angular momentum up to 8 (l functions) is always supported, both for the ECPs and for the
orbital basis including the derivative order: integrals up to l functions, first derivatives
up to k, and second derivatives up to i. The limit is available as ``pylibecpint.MAX_L``, and
exceeding it raises :class:`ValueError`. No build options are needed for this.

Build options
^^^^^^^^^^^^^

Two optional settings can be passed when building from source:

.. code-block:: bash

   CMAKE_ARGS="-DLIBECPINT_MAX_UNROL=2 -DLIBECPINT_RADIAL_CACHE_MB=8" pip install .

``LIBECPINT_MAX_UNROL`` (default 1) is the angular momentum below which integrals are evaluated
with generated, unrolled code instead of loops, at the cost of extra compile time.
``LIBECPINT_RADIAL_CACHE_MB`` (default 4; one of 1, 2, 4, 8, 16, 32 or 64) is the size of the
per-thread cache for radial factors.

ECP library
^^^^^^^^^^^

pylibecpint bundles libecpint's library of 91 ECP families, as MOLPRO-format ``.ecp`` files.
These are used by :meth:`~pylibecpint.ECPIntegrator.set_ecp_basis_from_library`
and are resolved automatically, so you do not need to specify a path.

The library includes the Stuttgart/Cologne ECPs (e.g. ``ecp10mdf``, ``ecp28mwb``,
``ecp60mdf``), the Karlsruhe ``def2-ecp`` and ``dhf-ecp``, ``lanl2dz``, ``crenbl``,
``crenbs``, ``sbkjc`` and ``lcecp-*``. To list every name:

.. code-block:: python

   import pylibecpint
   print(pylibecpint.available_ecps())

The files are in the ``raw`` subdirectory of :func:`~pylibecpint.get_data_dir`. Your own ECPs in
the same format (for example, a MOLPRO-format export from the
`Basis Set Exchange <https://www.basissetexchange.org>`_) can be loaded with
``share_dir`` pointing to a directory that contains ``raw/<name>.ecp``, or with
:meth:`pylibecpint.ECPBasis.addECP_from_file`.

Example: HBr ECP integrals
---------------------------

This example computes ECP integrals for HBr using the AVDZ-PP basis set with
a Stuttgart-Dresden 10-electron ECP on bromine. It is a port of the HBr test from
libecpint v1.0, and is one of pylibecpint's tests.

.. code-block:: python

   import numpy as np
   from pylibecpint import ECPIntegrator

   # Bond length in bohr
   R = 2.66

   # Gaussian basis: 12 shells on Br (at origin) + 5 shells on H (at z=R)
   # Each shell is defined by coordinates, exponents, coefficients,
   # angular momentum, and number of primitives.
   g_coords = np.array([
       # 12 Br shells at origin (each needs xyz = 0, 0, 0)
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       # 5 H shells at (0, 0, R)
       0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0, R,
   ])

   g_exps = np.array([
       2.808600e+03, 4.211800e+02, 5.034570e+01, 1.791330e+01,
       3.805310e+00, 1.749680e+00, 4.485550e-01, 1.644980e-01,
       2.808600e+03, 4.211800e+02, 5.034570e+01, 1.791330e+01,
       3.805310e+00, 1.749680e+00, 4.485550e-01, 1.644980e-01,
       4.485550e-01, 1.644980e-01, 5.020000e-02,
       1.057520e+02, 2.763680e+01, 6.596560e+00, 2.785220e+00,
       1.078120e+00, 3.935370e-01, 1.274690e-01,
       1.057520e+02, 2.763680e+01, 6.596560e+00, 2.785220e+00,
       1.078120e+00, 3.935370e-01, 1.274690e-01,
       1.274690e-01, 3.940000e-02,
       1.438650e+02, 4.611630e+01, 1.736940e+01, 6.951070e+00,
       2.756070e+00, 1.011780e+00, 4.291000e-01,
       4.291000e-01, 1.548000e-01,
       1.301000e+01, 1.962000e+00, 4.446000e-01, 1.220000e-01,
       1.220000e-01, 2.974000e-02, 7.270000e-01, 1.410000e-01,
   ])

   g_coefs = np.array([
        1.606000e-03,  8.393000e-03,  6.957800e-02, -3.899080e-01,
        6.944970e-01,  4.913540e-01,  2.263700e-02, -3.723000e-03,
       -6.350000e-04, -3.492000e-03, -2.519500e-02,  1.501130e-01,
       -3.662260e-01, -3.834220e-01,  7.144680e-01,  5.352530e-01,
        1.0, 1.0, 1.0,
        5.341000e-03, -8.308400e-02,  4.477660e-01,  5.506170e-01,
        1.235000e-01, -3.771000e-03,  2.278000e-03,
       -1.308000e-03,  2.292100e-02, -1.450290e-01, -2.090370e-01,
        9.373000e-02,  6.050210e-01,  4.571230e-01,
        1.0, 1.0,
        1.023700e-02,  7.608300e-02,  2.298070e-01,  4.033470e-01,
        4.097280e-01,  1.627900e-01,  6.991000e-03,
        1.0, 1.0,
        1.968500e-02,  1.379770e-01,  4.781480e-01,  5.012400e-01,
        1.0, 1.0, 1.0, 1.0,
   ])

   g_ams = np.array(
       [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1],
       dtype=np.int32,
   )
   g_lengths = np.array(
       [8, 8, 1, 1, 1, 7, 7, 1, 1, 7, 1, 1, 4, 1, 1, 1, 1],
       dtype=np.int32,
   )

   # ECP basis: one center at the Br position (origin)
   u_coords = np.array([0.0, 0.0, 0.0])
   u_exps = np.array([
       1.0, 70.024257, 31.178412, 7.156593,
       46.773471, 46.184120, 21.713858, 20.941792,
       50.698839, 50.644764, 15.447509, 15.500259,
       2.800391, 1.077480, 14.465606, 21.234065,
   ])
   u_coefs = np.array([
       0.0, 49.962834, 370.014205, 10.241439,
       99.112244, 198.253046, 28.261740, 56.623366,
       -18.605853, -27.923280, -0.379693, -0.780583,
       0.035968, 0.094397, -1.091269, -2.887691,
   ])
   u_ams = np.array(
       [4, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3],
       dtype=np.int32,
   )
   u_ns = np.array([2] * 16, dtype=np.int32)
   u_lengths = np.array([16], dtype=np.int32)

   # Create integrator and compute
   integrator = ECPIntegrator()
   integrator.set_gaussian_basis(g_coords, g_exps, g_coefs, g_ams, g_lengths)
   integrator.set_ecp_basis(u_coords, u_exps, u_coefs, u_ams, u_ns, u_lengths)
   integrator.init()
   integrals = integrator.compute_integrals()

   print(f"Integral matrix shape: {integrals.shape}")
   print(f"Trace: {np.trace(integrals):.6f}")

Running this produces::

   Integral matrix shape: (44, 44)
   Trace: 92.205935

Using the ECP library
^^^^^^^^^^^^^^^^^^^^^

Instead of specifying the ECP primitives manually, you can load them from
the bundled library. For the same HBr system with ECP10MDF on bromine
(atomic number 35):

.. code-block:: python

   integrator = ECPIntegrator()
   integrator.set_gaussian_basis(g_coords, g_exps, g_coefs, g_ams, g_lengths)
   integrator.set_ecp_basis_from_library(
       coords=np.array([0.0, 0.0, 0.0]),
       charges=np.array([35], dtype=np.int32),
       names=["ecp10mdf"],
   )
   integrator.init()
   integrals = integrator.compute_integrals()

The ``share_dir`` argument defaults to the bundled library, so there is no need
to locate the files manually.

Computing derivatives
^^^^^^^^^^^^^^^^^^^^^

To compute first or second derivatives, pass the derivative order to
:meth:`~pylibecpint.ECPIntegrator.init`. Derivatives are with respect to the
coordinates of each distinct atom, in the order in which the atoms first appear in
the Gaussian basis.

Usually only the contraction of the derivatives with a density matrix is needed,
for example for a gradient or Hessian. :meth:`~pylibecpint.ECPIntegrator.compute_gradient`
and :meth:`~pylibecpint.ECPIntegrator.compute_hessian` compute this directly,
without ever storing the derivative matrices:

.. code-block:: python

   integrator.init(deriv=2)
   density = ...  # (ncart, ncart), in the order of the integrals

   gradient = integrator.compute_gradient(density)  # shape (3 * natoms,)
   hessian = integrator.compute_hessian(density)    # shape (3 * natoms, 3 * natoms)

   # gradient[3 * X + i] is the derivative with respect to coordinate i (x, y, z)
   # of atom X, and hessian[3 * X + i, 3 * Y + j] likewise

This is the recommended way to get derivatives for anything but small systems.
The full derivative matrices take memory proportional to natoms * ncart\ :sup:`2`
(first derivatives) and natoms\ :sup:`2` * ncart\ :sup:`2` (second derivatives),
e.g. about 7 GB and 800 GB respectively for 77 atoms with 1900 basis functions.
When you do need them:

.. code-block:: python

   integrator.init(deriv=1)
   first_derivs = integrator.compute_first_derivs()

   # first_derivs is a list of 3*natoms arrays, in order:
   # [Ax, Ay, Az, Bx, By, Bz, ...]
   print(f"Number of derivative components: {len(first_derivs)}")
   print(f"Each component shape: {first_derivs[0].shape}")

:meth:`~pylibecpint.ECPIntegrator.compute_second_derivs` returns the
3*natoms*(3*natoms+1)/2 unique second derivatives; its documentation gives the order.

Magnetic field derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^

For gauge-including atomic orbitals (London orbitals),
:meth:`~pylibecpint.ECPIntegrator.compute_bfield_derivs` gives the derivatives of the
ECP integrals with respect to a uniform magnetic field. It requires ``init(deriv=1)``
or higher:

.. code-block:: python

   integrator.init(deriv=1)
   bx, by, bz = integrator.compute_bfield_derivs()

Each matrix is real and antisymmetric, and the derivative of the integrals is i/2
times it. The conventions are described in the
`libecpint documentation <https://libecpint.readthedocs.io>`_.

Low-level API
^^^^^^^^^^^^^

For finer control, you can compute integrals over individual shell pairs
using :class:`~pylibecpint.ECPIntegral`:

.. code-block:: python

   from pylibecpint import ECP, ECPIntegral, GaussianShell

   # Create an ECP at position A
   A = np.array([0.5, 1.5, 0.9])
   ecp = ECP(A)
   ecp.addPrimitive(2, 0, 2.0, 0.44)
   ecp.addPrimitive(2, 1, 1.2, 0.42)

   # Create two Gaussian shells
   shell_a = GaussianShell(A, 1)  # p-type at A
   shell_a.addPrim(0.15, 0.3)
   shell_a.addPrim(2.49, 0.9)

   B = np.array([-0.5, -0.5, -0.5])
   shell_b = GaussianShell(B, 0)  # s-type at B
   shell_b.addPrim(1.0, 1.0)

   # Compute: maxLB=1 (max L in orbital basis), maxLU=1 (max L in ECP)
   ecpint = ECPIntegral(1, 1)
   result = ecpint.compute_shell_pair(ecp, shell_a, shell_b)

   # result is a TwoIndex; convert to numpy
   arr = result.to_numpy()
   print(f"Shell pair integral shape: {arr.shape}")

Notes
-----

- All coordinates must be in **bohr** (atomic units).
- Integrals are returned in **canonical Cartesian order**
  (e.g. for L=2: xx, xy, xz, yy, yz, zz).
- Calling a method that the integrator was not initialised for (for example
  :meth:`~pylibecpint.ECPIntegrator.compute_gradient` after ``init(deriv=0)``) raises
  :class:`RuntimeError`, and inputs of the wrong size raise :class:`ValueError`.
- ECP powers (``n``) follow the Stuttgart-Dresden convention where ``n=2``
  is the default (includes the r\ :sup:`2` Jacobian factor). Internally the
  library stores ``n - 2``.
- For further details on the underlying C++ library, see the
  `libecpint documentation <https://libecpint.readthedocs.io>`_.
