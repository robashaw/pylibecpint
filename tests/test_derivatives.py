"""High-level derivatives: matrices, density contractions and magnetic field derivatives."""

import numpy as np
import pytest

from pylibecpint import ECPIntegrator

from hbr import (
    G_AMS,
    G_COEFS,
    G_COORDS,
    G_EXPS,
    G_LENGTHS,
    U_AMS,
    U_COEFS,
    U_COORDS,
    U_EXPS,
    U_LENGTHS,
    U_NS,
    make_integrator,
)

NATOMS = 2
NCART = 44


@pytest.fixture(scope="module")
def hess_integrator():
    return make_integrator(deriv=2)


@pytest.fixture(scope="module")
def first_derivs(hess_integrator):
    return hess_integrator.compute_first_derivs()


@pytest.fixture(scope="module")
def second_derivs(hess_integrator):
    return hess_integrator.compute_second_derivs()


@pytest.fixture(scope="module")
def density():
    rng = np.random.default_rng(42)
    p = rng.standard_normal((NCART, NCART))
    return p + p.T


def unpack_second_derivs(packed, natoms):
    """(3*natoms, 3*natoms) array of matrices from the packed order of compute_second_derivs."""
    n = 3 * natoms
    full = np.empty((n, n) + packed[0].shape)
    idx = 0
    for a in range(natoms):
        for b in range(a, natoms):
            if a == b:
                pairs = [(i, j) for i in range(3) for j in range(i, 3)]
            else:
                pairs = [(i, j) for i in range(3) for j in range(3)]
            for i, j in pairs:
                full[3 * a + i, 3 * b + j] = packed[idx]
                full[3 * b + j, 3 * a + i] = packed[idx]
                idx += 1
    assert idx == len(packed)
    return full


def test_first_derivs_shape(hess_integrator, first_derivs):
    assert hess_integrator.natoms == NATOMS
    assert len(first_derivs) == 3 * NATOMS
    assert all(d.shape == (NCART, NCART) for d in first_derivs)


def test_second_derivs_shape(second_derivs):
    n = 3 * NATOMS
    assert len(second_derivs) == n * (n + 1) // 2


def test_translational_invariance(first_derivs, second_derivs):
    # Moving every atom together leaves the integrals unchanged
    d = np.array(first_derivs)
    total = d[0::3].sum(axis=0), d[1::3].sum(axis=0), d[2::3].sum(axis=0)
    scale = np.abs(d).max()
    for t in total:
        assert np.abs(t).max() < 1e-8 * scale
    full = unpack_second_derivs(second_derivs, NATOMS)
    rows = full.reshape(NATOMS, 3, 3 * NATOMS, NCART, NCART).sum(axis=0)
    assert np.abs(rows).max() < 1e-8 * np.abs(full).max()


def test_gradient_matches_first_derivs(hess_integrator, first_derivs, density):
    gradient = hess_integrator.compute_gradient(density)
    expected = np.array([np.sum(density * d) for d in first_derivs])
    assert gradient.shape == (3 * NATOMS,)
    np.testing.assert_allclose(gradient, expected, rtol=1e-10, atol=1e-10)


def test_gradient_nonsymmetric_density(hess_integrator, first_derivs):
    rng = np.random.default_rng(7)
    p = rng.standard_normal((NCART, NCART))
    gradient = hess_integrator.compute_gradient(p)
    expected = np.array([np.sum(p * d) for d in first_derivs])
    np.testing.assert_allclose(gradient, expected, rtol=1e-10, atol=1e-10)


def test_hessian_matches_second_derivs(hess_integrator, second_derivs, density):
    hessian = hess_integrator.compute_hessian(density)
    full = unpack_second_derivs(second_derivs, NATOMS)
    expected = np.einsum("ijkl,kl->ij", full, density)
    assert hessian.shape == (3 * NATOMS, 3 * NATOMS)
    np.testing.assert_allclose(hessian, hessian.T, rtol=0, atol=1e-10)
    np.testing.assert_allclose(hessian, expected, rtol=1e-10, atol=1e-10)


def test_repeated_derivs_do_not_accumulate():
    integrator = make_integrator(deriv=2)
    first = integrator.compute_first_derivs()
    again = integrator.compute_first_derivs()
    assert len(again) == len(first) == 3 * NATOMS
    for a, b in zip(first, again):
        np.testing.assert_array_equal(a, b)
    second = integrator.compute_second_derivs()
    again = integrator.compute_second_derivs()
    assert len(again) == len(second)
    for a, b in zip(second, again):
        np.testing.assert_array_equal(a, b)


def test_bfield_derivs():
    integrator = make_integrator(deriv=1)
    derivs = integrator.compute_bfield_derivs()
    assert len(derivs) == 3
    for d in derivs:
        assert d.shape == (NCART, NCART)
        np.testing.assert_allclose(d, -d.T, rtol=0, atol=1e-12)
    assert max(np.abs(d).max() for d in derivs) > 1e-6


def test_update_coords_matches_fresh_integrator():
    shift = np.array([0.01, -0.02, 0.03, 0.0, 0.015, -0.01])  # Br, H
    atom_of_shell = (G_COORDS.reshape(-1, 3)[:, 2] != 0).astype(int)
    new_g = (G_COORDS.reshape(-1, 3) + shift.reshape(2, 3)[atom_of_shell]).ravel()
    new_u = U_COORDS + shift[:3]

    moved = make_integrator(deriv=1)
    moved.compute_integrals()
    moved.update_gaussian_basis_coords(new_g)
    moved.update_ecp_basis_coords(new_u)

    fresh = ECPIntegrator()
    fresh.set_gaussian_basis(new_g, G_EXPS, G_COEFS, G_AMS, G_LENGTHS)
    fresh.set_ecp_basis(new_u, U_EXPS, U_COEFS, U_AMS, U_NS, U_LENGTHS)
    fresh.init(1)

    np.testing.assert_array_equal(moved.compute_integrals(), fresh.compute_integrals())
    for a, b in zip(moved.compute_first_derivs(), fresh.compute_first_derivs()):
        np.testing.assert_array_equal(a, b)


# --- errors instead of undefined behaviour ---


def test_compute_before_init():
    integrator = ECPIntegrator()
    integrator.set_gaussian_basis(G_COORDS, G_EXPS, G_COEFS, G_AMS, G_LENGTHS)
    integrator.set_ecp_basis(U_COORDS, U_EXPS, U_COEFS, U_AMS, U_NS, U_LENGTHS)
    assert integrator.natoms == 0 and integrator.deriv == 0
    with pytest.raises(RuntimeError, match="init"):
        integrator.compute_integrals()


def test_init_before_bases():
    integrator = ECPIntegrator()
    integrator.set_gaussian_basis(G_COORDS, G_EXPS, G_COEFS, G_AMS, G_LENGTHS)
    with pytest.raises(RuntimeError):
        integrator.init()


@pytest.mark.parametrize(
    "deriv, method, args",
    [
        (0, "compute_first_derivs", ()),
        (0, "compute_gradient", (np.zeros((NCART, NCART)),)),
        (0, "compute_bfield_derivs", ()),
        (1, "compute_second_derivs", ()),
        (1, "compute_hessian", (np.zeros((NCART, NCART)),)),
    ],
)
def test_derivative_order_too_low(deriv, method, args):
    integrator = make_integrator(deriv=deriv)
    with pytest.raises(RuntimeError, match="deriv"):
        getattr(integrator, method)(*args)


@pytest.mark.parametrize("shape", [(NCART,), (NCART, NCART - 1), (NCART * NCART,)])
def test_density_wrong_shape(hess_integrator, shape):
    with pytest.raises(ValueError, match="density"):
        hess_integrator.compute_gradient(np.zeros(shape))
    with pytest.raises(ValueError, match="density"):
        hess_integrator.compute_hessian(np.zeros(shape))


def test_basis_length_mismatch():
    integrator = ECPIntegrator()
    with pytest.raises(ValueError, match="exponents"):
        integrator.set_gaussian_basis(G_COORDS, G_EXPS[:-1], G_COEFS, G_AMS, G_LENGTHS)
    with pytest.raises(ValueError, match="coords"):
        integrator.set_gaussian_basis(G_COORDS[:-3], G_EXPS, G_COEFS, G_AMS, G_LENGTHS)
    with pytest.raises(ValueError, match="ns"):
        integrator.set_ecp_basis(U_COORDS, U_EXPS, U_COEFS, U_AMS, U_NS[:-1], U_LENGTHS)


def test_update_coords_length_mismatch():
    integrator = make_integrator()
    with pytest.raises(ValueError):
        integrator.update_gaussian_basis_coords(G_COORDS[:-3])
    with pytest.raises(ValueError):
        integrator.update_ecp_basis_coords(np.zeros(6))
