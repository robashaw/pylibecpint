"""Port of libecpint api_test1: high-level ECPIntegrator with HBr/AVDZ-PP basis."""

from conftest import check_against_reference
from hbr import make_integrator


def test_api_integrals():
    integrator = make_integrator()
    integrals = integrator.compute_integrals()

    assert integrals.shape == (44, 44)
    check_against_reference(integrals, "api_test1.output")
