"""Tests for the Python convenience layer."""

import os

import numpy as np

from pylibecpint import (
    ECP,
    ECPBasis,
    ECPIntegrator,
    GaussianECP,
    GaussianShell,
    TwoIndex,
    get_data_dir,
)


def test_get_data_dir():
    d = get_data_dir()
    assert os.path.isdir(d)
    xmls = [f for f in os.listdir(d) if f.endswith(".xml")]
    assert len(xmls) > 0


def test_gaussian_shell_from_list():
    shell = GaussianShell(np.array([1.0, 2.0, 3.0]), 2)
    shell.addPrim(1.5, 0.8)
    assert shell.am() == 2
    assert shell.nprimitive() == 1
    assert shell.ncartesian() == 6
    center = shell.center()
    np.testing.assert_allclose(center, [1.0, 2.0, 3.0])


def test_ecp_construction():
    ecp = ECP(np.array([0.0, 0.0, 0.0]))
    ecp.addPrimitive(2, 0, 1.0, 0.5)
    ecp.addPrimitive(2, 1, 2.0, 0.3)
    assert ecp.getN() == 2
    assert ecp.getL() == 1
    np.testing.assert_allclose(ecp.center(), [0.0, 0.0, 0.0])


def test_gaussian_ecp():
    g = GaussianECP(2, 1, 3.5, 0.7)
    # n is stored as n-2 internally (ECP convention: n=2 means r^0)
    assert g.n == 0
    assert g.l == 1
    assert g.a == 3.5
    assert g.d == 0.7


def test_ecpbasis():
    basis = ECPBasis()
    ecp = ECP(np.array([0.0, 0.0, 0.0]))
    ecp.addPrimitive(2, 0, 1.0, 1.0)
    basis.addECP(ecp, 0)
    assert basis.getN() == 1
    assert basis.getAtom(0) == 0


def test_two_index():
    t = TwoIndex(3, 4, 0.0)
    assert t.shape == (3, 4)
    t[1, 2] = 5.0
    assert t[1, 2] == 5.0
    arr = t.to_numpy()
    assert arr.shape == (3, 4)
    assert arr[1, 2] == 5.0


def test_integrator_list_inputs():
    integrator = ECPIntegrator()
    integrator.set_gaussian_basis(
        [0.0, 0.0, 0.0],
        [1.0],
        [1.0],
        [0],
        [1],
    )
    integrator.set_ecp_basis(
        [0.0, 0.0, 0.0],
        [1.0],
        [1.0],
        [0],
        [2],
        [1],
    )
    integrator.init()
    result = integrator.compute_integrals()
    assert result.shape == (1, 1)


def test_tolerance_properties():
    integrator = ECPIntegrator()
    assert integrator.tolerance == 1e-12
    integrator.tolerance = 1e-10
    assert integrator.tolerance == 1e-10
