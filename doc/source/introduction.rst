.. _introduction:

Introduction
============

pylibecpint provides Python bindings for
`libecpint <https://github.com/robashaw/libecpint>`_, enabling the
evaluation of integrals over effective core potentials (ECPs) directly from
Python. All inputs and outputs use NumPy arrays.

Installation
------------

From PyPI
^^^^^^^^^

.. code-block:: bash

   pip install pylibecpint

From source
^^^^^^^^^^^

Building from source requires a C++ compiler and CMake >= 3.15.

.. code-block:: bash

   git clone --recurse-submodules https://github.com/robashaw/pylibecpint.git
   cd pylibecpint
   pip install .

Configuration
-------------

Maximum angular momentum
^^^^^^^^^^^^^^^^^^^^^^^^^

The library is compiled with a maximum angular momentum (``LIBECPINT_MAX_L``,
default 5, supporting up to h-type functions) and a maximum unrolling level
(``LIBECPINT_MAX_UNROL``, default 1). If your calculation requires higher
angular momentum, rebuild from source with:

.. code-block:: bash

   CMAKE_ARGS="-DLIBECPINT_MAX_L=7 -DLIBECPINT_MAX_UNROL=2" pip install .

Higher values of ``MAX_L`` increase compile time and binary size. Higher values
of ``MAX_UNROL`` (up to ``MAX_L``) replace loop-based integral evaluation with
generated unrolled code for better performance, at the cost of additional
compile time.

ECP library
^^^^^^^^^^^

pylibecpint bundles the ECP definitions from libecpint as XML data files.
These are used by :meth:`~pylibecpint.ECPIntegrator.set_ecp_basis_from_library`
and are resolved automatically — you do not need to specify a path.

The available ECPs are:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Name
     - Atoms
   * - ECP10MDF
     - K -- Kr (Z = 19 -- 36)
   * - ECP28MDF
     - Rb -- Xe (Z = 37 -- 54)
   * - ECP46MDF
     - Cs, Ba (Z = 55, 56)
   * - ECP60MDF
     - Hf -- Rn (Z = 72 -- 86), Ac -- U (Z = 89 -- 92)
   * - ECP78MDF
     - Fr, Ra (Z = 87, 88)
   * - LANL2DZ
     - Na -- La (Z = 11 -- 57), Hf -- Bi (Z = 72 -- 83), U -- Pu (Z = 92 -- 94)

Example: HBr ECP integrals
---------------------------

This example computes ECP integrals for HBr using the AVDZ-PP basis set with
a Stuttgart-Dresden 10-electron ECP on bromine. It is a direct port of
``tests/lib/api_test1`` from the C++ library.

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
       # 5 H shells at (R, 0, 0)
       R, 0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0, R, 0.0, 0.0,
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
   Trace: 132.810798

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

The ``share_dir`` argument defaults to the bundled data directory, so there is
no need to locate the XML files manually.

Computing derivatives
^^^^^^^^^^^^^^^^^^^^^

To compute first or second derivatives, pass the derivative order to
:meth:`~pylibecpint.ECPIntegrator.init`:

.. code-block:: python

   integrator.init(deriv=1)
   integrals = integrator.compute_integrals()
   first_derivs = integrator.compute_first_derivs()

   # first_derivs is a list of 3*natoms arrays, in order:
   # [Ax, Ay, Az, Bx, By, Bz, ...]
   # where atom order matches set_gaussian_basis.
   print(f"Number of derivative components: {len(first_derivs)}")
   print(f"Each component shape: {first_derivs[0].shape}")

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
- ECP powers (``n``) follow the Stuttgart-Dresden convention where ``n=2``
  is the default (includes the r\ :sup:`2` Jacobian factor). Internally the
  library stores ``n - 2``.
- For further details on the underlying C++ library, see the
  `libecpint documentation <https://libecpint.readthedocs.io>`_.
