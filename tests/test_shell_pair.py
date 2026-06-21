"""Port of libecpint int_test1: low-level ECPIntegral shell-pair computation."""

import numpy as np

from pylibecpint import ECP, ECPIntegral, GaussianShell

from conftest import check_against_reference


A = np.array([0.5, 1.5, 0.9])
B = np.array([-0.5, -0.5, -0.5])


def _make_ecp():
    ecp = ECP(A)
    ecp.addPrimitive(2, 4, 1.0, 0.0)
    ecp.addPrimitive(2, 0, 2.015185434, 0.437903605)
    ecp.addPrimitive(2, 1, 1.225772391, 0.418499875)
    ecp.addPrimitive(2, 2, 3.415601371, 10.482661239)
    ecp.addPrimitive(2, 3, 3.574403408, -4.478284427, True)
    return ecp


def _make_shell_a():
    shell = GaussianShell(A, 1)
    shell.addPrim(0.15, 0.3)
    shell.addPrim(2.49, 0.9)
    shell.addPrim(10.11, 1.2)
    return shell


def _make_shell_b():
    shell = GaussianShell(B, 3)
    shell.addPrim(0.02, 1.1)
    shell.addPrim(1.41, 0.9)
    shell.addPrim(8.71, 0.4)
    return shell


def test_shell_pair_integrals():
    ecp = _make_ecp()
    shell_a = _make_shell_a()
    shell_b = _make_shell_b()
    ecpint = ECPIntegral(3, 4)

    flat_result = []

    # (A|U|A)
    result = ecpint.compute_shell_pair(ecp, shell_a, shell_a)
    flat_result.extend(result.to_numpy().ravel())

    # (A|U|B)
    result = ecpint.compute_shell_pair(ecp, shell_a, shell_b)
    flat_result.extend(result.to_numpy().ravel())

    # (B|U|B)
    result = ecpint.compute_shell_pair(ecp, shell_b, shell_b)
    flat_result.extend(result.to_numpy().ravel())

    check_against_reference(np.array(flat_result), "test1.output")
