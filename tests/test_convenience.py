"""Tests for the Python convenience layer."""

import os

import numpy as np
import pytest

import pylibecpint
from pylibecpint import (
    ECP,
    ECPBasis,
    ECPIntegrator,
    GaussianECP,
    GaussianShell,
    TwoIndex,
    available_ecps,
    get_data_dir,
)


def test_get_data_dir():
    d = get_data_dir()
    raw = os.path.join(d, "raw")
    assert os.path.isdir(raw)
    ecps = [f for f in os.listdir(raw) if f.endswith(".ecp")]
    assert len(ecps) > 0


def test_available_ecps():
    names = available_ecps()
    assert names == sorted(names)
    for name in ["ecp10mdf", "ecp28mdf", "ecp60mdf", "lanl2dz", "def2-ecp", "crenbl", "sbkjc"]:
        assert name in names
    assert len(names) >= 91


def test_module_constants():
    assert pylibecpint.MAX_L == 8
    assert pylibecpint.__libecpint_version__.startswith("1.1")
    assert pylibecpint.__version__.startswith("1.1")


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


def test_two_index_out_of_range():
    t = TwoIndex(2, 2, 0.0)
    with pytest.raises(IndexError):
        t[2, 0]
    with pytest.raises(IndexError):
        t[0, -1] = 1.0


def test_ecp_index_out_of_range():
    ecp = ECP(np.zeros(3))
    ecp.addPrimitive(2, 0, 1.0, 1.0)
    with pytest.raises(IndexError):
        ecp.getGaussian(1)
    with pytest.raises(IndexError):
        ECPBasis().getECP(0)


def test_center_wrong_size():
    with pytest.raises(ValueError):
        GaussianShell(np.zeros(2), 0)
    with pytest.raises(ValueError):
        ECP(np.zeros(4))
