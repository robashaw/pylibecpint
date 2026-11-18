"""Angular momentum limits (port of libecpint's Api.AngularMomentumLimit)."""

import numpy as np
import pytest

from pylibecpint import ECP, MAX_L, ECPIntegral, ECPIntegrator


def _integrator(ecp_l=0):
    integrator = ECPIntegrator()
    integrator.set_gaussian_basis([0.0, 0.0, 1.0], [1.0], [1.0], [7], [1])
    integrator.set_ecp_basis([0.0, 0.0, 0.0], [2.0], [1.0], [ecp_l], [2], [1])
    return integrator


def test_max_l():
    assert MAX_L == 8


def test_basis_plus_deriv_limit():
    _integrator().init(1)  # k shell, first derivatives: L = 8
    with pytest.raises(ValueError):
        _integrator().init(2)  # second derivatives would need L = 9


def test_ecp_limit():
    with pytest.raises(ValueError):
        _integrator(ecp_l=MAX_L + 1)
    with pytest.raises(ValueError):
        ECP(np.zeros(3)).addPrimitive(2, MAX_L + 1, 1.0, 1.0)


def test_ecp_power_limit():
    with pytest.raises(ValueError):
        ECP(np.zeros(3)).addPrimitive(5, 0, 1.0, 1.0)


def test_ecpintegral_limit():
    ECPIntegral(MAX_L, MAX_L)
    with pytest.raises(ValueError):
        ECPIntegral(MAX_L, 2, deriv=1)
    with pytest.raises(ValueError):
        ECPIntegral(2, MAX_L + 1)
