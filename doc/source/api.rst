.. _api:

API Reference
=============

High-level API
--------------

.. autoclass:: pylibecpint.ECPIntegrator
   :members:
   :undoc-members:

.. autofunction:: pylibecpint.get_data_dir

Low-level API
-------------

These classes are re-exported from the C++ extension module
(``pylibecpint._core``) and provide direct access to the underlying
libecpint objects.

ECPIntegral
^^^^^^^^^^^

.. autoclass:: pylibecpint.ECPIntegral
   :members:
   :undoc-members:

Data structures
---------------

GaussianShell
^^^^^^^^^^^^^

.. autoclass:: pylibecpint.GaussianShell
   :members:
   :undoc-members:

ECP
^^^

.. autoclass:: pylibecpint.ECP
   :members:
   :undoc-members:

GaussianECP
^^^^^^^^^^^

.. autoclass:: pylibecpint.GaussianECP
   :members:
   :undoc-members:

ECPBasis
^^^^^^^^

.. autoclass:: pylibecpint.ECPBasis
   :members:
   :undoc-members:

TwoIndex
^^^^^^^^

.. autoclass:: pylibecpint.TwoIndex
   :members:
   :undoc-members:

ShellPairData
^^^^^^^^^^^^^

.. autoclass:: pylibecpint.ShellPairData
   :members:
   :undoc-members:
