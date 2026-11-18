"""Loading ECPs from the bundled MOLPRO-format library."""

import os

import numpy as np
import pytest

from pylibecpint import ECPBasis, ECPIntegrator, get_data_dir

from conftest import check_against_reference
from hbr import G_AMS, G_COEFS, G_COORDS, G_EXPS, G_LENGTHS, U_COORDS, make_integrator


def _library_integrator(names=("ecp10mdf",), charges=(35,), share_dir=None):
    integrator = ECPIntegrator()
    integrator.set_gaussian_basis(G_COORDS, G_EXPS, G_COEFS, G_AMS, G_LENGTHS)
    integrator.set_ecp_basis_from_library(U_COORDS, charges, list(names), share_dir)
    return integrator


def test_library_matches_manual_ecp():
    # The library's ECP10MDF for Br is the ECP that api_test1 enters by hand
    integrator = _library_integrator()
    integrator.init()
    integrals = integrator.compute_integrals()
    check_against_reference(integrals, "api_test1.output")
    np.testing.assert_allclose(integrals, make_integrator().compute_integrals(), rtol=0, atol=1e-14)


def test_library_explicit_share_dir():
    integrator = _library_integrator(share_dir=get_data_dir())
    integrator.init()
    assert integrator.compute_integrals().shape == (44, 44)


def test_library_unknown_name():
    with pytest.raises(RuntimeError):
        _library_integrator(names=("not-an-ecp",))


def test_library_atom_missing():
    # ECP10MDF has no ECP for hydrogen
    with pytest.raises(RuntimeError):
        _library_integrator(charges=(1,))


def test_library_names_length_mismatch():
    with pytest.raises(ValueError):
        _library_integrator(names=("ecp10mdf", "ecp10mdf"))


def test_add_ecp_from_file():
    basis = ECPBasis()
    path = os.path.join(get_data_dir(), "raw", "ecp10mdf.ecp")
    basis.addECP_from_file(35, np.zeros(3), path)
    assert basis.getN() == 1
    assert basis.getECP(0).getN() == 16


def test_add_ecp_from_missing_file(tmp_path):
    basis = ECPBasis()
    with pytest.raises(RuntimeError):
        basis.addECP_from_file(35, np.zeros(3), str(tmp_path / "missing.ecp"))
