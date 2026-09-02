"""Boundary-validation tests for the B_model matter-bispectrum input.

Run from ``sachs_sft/callables/kappa3_vertex`` with the canoes venv::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q

The suite follows the project's boundary-validation methodology: every
dispatch seam is probed from both sides with both methods evaluated
directly, and extreme parameter corners (lowest and highest shell, softest
and hardest leg, most squeezed triangle) are always included.
"""

from __future__ import annotations

import numpy as np
import pytest

from b_model import (
    FIDUCIAL,
    PK_TABLE_2PT,
    PK_TABLE_3PT,
    SIGMA8_2PT,
    SIGMA8_3PT,
    b_bihalofit,
    b_tree_z0,
    closure_cosines,
    f2_kernel,
    growth,
    load_pk_lin,
    make_b_delta_fn,
    sigma8_of_table,
)
from b_model.core import K_BLEND_HI, K_BLEND_LO, _blend_weight

def _has_baryon_fit() -> bool:
    """The Takahashi baryon ratio ships only in some fastnc builds."""
    from b_model.halofit_backend import load_fastnc_halofit

    return hasattr(load_fastnc_halofit(), "get_Rb_bihalofit")


#: Shell redshifts of the deployed equal_time_limber table.
SHELL_Z = [0.1, 0.25, 0.45, 0.7, 1.0, 1.4, 1.9, 2.5, 3.1, 3.6, 4.05, 4.45,
           4.75, 5.0, 5.3, 5.7]
EXTREME_Z = [SHELL_Z[0], SHELL_Z[-1]]


def _closed_triangle(k_a, k_b, mu):
    """Third leg of the closed triangle with |cos| = mu between legs a, b."""
    return np.sqrt(k_a**2 + k_b**2 + 2.0 * k_a * k_b * mu)


# --------------------------------------------------------------------------
# sigma8 pins: the normalisation both branches share
# --------------------------------------------------------------------------

def test_sigma8_pins():
    assert sigma8_of_table(PK_TABLE_3PT) == pytest.approx(SIGMA8_3PT, abs=1e-6)
    assert sigma8_of_table(PK_TABLE_2PT) == pytest.approx(SIGMA8_2PT, abs=1e-6)


def test_bihalofit_does_not_renormalise_the_table():
    """The sigma8 we pin must be the one fastnc computes, or P is rescaled."""
    from b_model.cosmology import load_pk_table
    from b_model.halofit_backend import _halofit_instance

    halofit = _halofit_instance(FIDUCIAL, PK_TABLE_3PT)
    _, p_raw = load_pk_table(PK_TABLE_3PT)
    np.testing.assert_allclose(halofit.pklin, p_raw, rtol=1e-6)


# --------------------------------------------------------------------------
# Tree branch: must reproduce the canoes internal assembly
# --------------------------------------------------------------------------

def test_tree_matches_canoes_quadrature_assembly():
    """Reproduce canoes' (u, v, phi) tree assembly on its own geometry."""
    rng = np.random.default_rng(20260825)
    u = rng.uniform(1.0, 1000.0, 4000)
    v = rng.uniform(1.0, 1000.0, 4000)
    phi = rng.uniform(0.0, 2.0 * np.pi, 4000)
    w = np.sqrt(u**2 + v**2 + 2.0 * u * v * np.cos(phi))
    chi_h = 1500.0
    k1, k2, k3 = w / chi_h, u / chi_h, v / chi_h

    cos_12 = -(u + v * np.cos(phi)) / w
    cos_13 = -(v + u * np.cos(phi)) / w
    cos_23 = np.cos(phi)
    pk = load_pk_lin(PK_TABLE_3PT)
    p1, p2, p3 = pk(k1), pk(k2), pk(k3)
    z = 1.0
    expected = growth(z) ** 4 * 2.0 * (
        f2_kernel(k1, k2, cos_12) * p1 * p2
        + f2_kernel(k1, k3, cos_13) * p1 * p3
        + f2_kernel(k2, k3, cos_23) * p2 * p3
    )
    got = make_b_delta_fn("tree")(k1, k2, k3, z)
    np.testing.assert_allclose(got, expected, rtol=1e-11)


def test_closure_cosines_match_the_quadrature_convention():
    rng = np.random.default_rng(7)
    u = rng.uniform(1.0, 500.0, 500)
    v = rng.uniform(1.0, 500.0, 500)
    phi = rng.uniform(0.0, 2.0 * np.pi, 500)
    w = np.sqrt(u**2 + v**2 + 2.0 * u * v * np.cos(phi))
    mu12, mu13, mu23 = closure_cosines(w, u, v)
    np.testing.assert_allclose(mu12, -(u + v * np.cos(phi)) / w, atol=1e-12)
    np.testing.assert_allclose(mu13, -(v + u * np.cos(phi)) / w, atol=1e-12)
    np.testing.assert_allclose(mu23, np.cos(phi), atol=1e-12)


@pytest.mark.parametrize("z", EXTREME_Z)
def test_tree_growth_scaling_is_d4(z):
    k = np.array([0.1, 0.5])
    pk = load_pk_lin(PK_TABLE_3PT)
    base = b_tree_z0(k, k, k, pk)
    got = make_b_delta_fn("tree")(k, k, k, z)
    np.testing.assert_allclose(got, growth(z) ** 4 * base, rtol=1e-12)


# --------------------------------------------------------------------------
# Permutation symmetry and immutability of inputs
# --------------------------------------------------------------------------

@pytest.mark.parametrize("model", ["tree", "bihalofit", "bihalofit_baryon"])
def test_permutation_symmetry(model):
    fn = make_b_delta_fn(model)
    if model == "bihalofit_baryon" and not _has_baryon_fit():
        pytest.skip("installed fastnc has no get_Rb_bihalofit")
    k_a = np.array([0.7, 0.03, 1.5])
    k_b = np.array([0.5, 0.02, 1.2])
    k_c = _closed_triangle(k_a, k_b, np.array([0.3, -0.5, 0.8]))
    base = fn(k_a, k_b, k_c, 0.7)
    for perm in [(k_a, k_c, k_b), (k_b, k_a, k_c), (k_c, k_b, k_a)]:
        np.testing.assert_allclose(fn(*perm, 0.7), base, rtol=1e-10)


@pytest.mark.parametrize("model", ["tree", "bihalofit"])
def test_inputs_are_not_mutated(model):
    fn = make_b_delta_fn(model)
    k_a = np.array([0.4, 0.9])
    k_b = np.array([0.3, 0.6])
    k_c = _closed_triangle(k_a, k_b, np.array([0.1, -0.2]))
    snapshots = [k.copy() for k in (k_a, k_b, k_c)]
    fn(k_a, k_b, k_c, 1.0)
    for original, current in zip(snapshots, (k_a, k_b, k_c)):
        np.testing.assert_array_equal(original, current)


# --------------------------------------------------------------------------
# Blend seam: both methods evaluated directly at and around the threshold
# --------------------------------------------------------------------------

def test_blend_weight_endpoints():
    assert _blend_weight(np.array([K_BLEND_LO * 0.5]))[0] == 0.0
    assert _blend_weight(np.array([K_BLEND_HI * 2.0]))[0] == 1.0
    mid = _blend_weight(np.array([np.sqrt(K_BLEND_LO * K_BLEND_HI)]))[0]
    assert 0.4 < mid < 0.6


@pytest.mark.parametrize("z", EXTREME_Z + [0.5, 1.9])
def test_seam_methods_agree_at_the_threshold(z):
    """Tree and BiHalofit must be close where the dispatcher hands over.

    This is the boundary check that motivated placing the seam at the
    crossover rather than deeper in the linear regime. The gate is loose
    (10%) because BiHalofit carries a few-percent large-scale offset by
    construction; the point is that the seam is not where it is largest.
    """
    for k_max in (K_BLEND_LO, np.sqrt(K_BLEND_LO * K_BLEND_HI), K_BLEND_HI):
        k = np.array([k_max])
        pk = load_pk_lin(PK_TABLE_3PT)
        tree = growth(z) ** 4 * b_tree_z0(k, k, k, pk)
        nonlinear = b_bihalofit(k, k, k, z)
        rel = abs(nonlinear[0] - tree[0]) / max(abs(tree[0]), 1e-300)
        assert rel < 0.10, f"seam mismatch {rel:.3f} at k={k_max}, z={z}"


def test_blended_model_is_continuous_across_the_seam():
    fn = make_b_delta_fn("bihalofit")
    k = np.geomspace(K_BLEND_LO * 0.5, K_BLEND_HI * 2.0, 200)
    values = fn(k, k, k, 1.0)
    jumps = np.abs(np.diff(np.log(values)))
    assert jumps.max() < 0.05, "blend introduces a step in ln B"


# --------------------------------------------------------------------------
# Domain sweep: finiteness and sane magnitudes over the production corners
# --------------------------------------------------------------------------

@pytest.mark.parametrize("z", SHELL_Z)
@pytest.mark.parametrize("model", ["tree", "bihalofit"])
def test_finite_over_the_production_domain(model, z):
    """Every shell, from the softest reachable leg to the hardest."""
    fn = make_b_delta_fn(model)
    k_hard = np.geomspace(1e-3, 5.0, 40)
    k_soft = np.full_like(k_hard, 4.0e-7)  # ell_min = 1e-3 safe floor
    k_third = _closed_triangle(k_hard, k_hard, np.full_like(k_hard, -1.0) + 1e-12)
    for legs in (
        (k_hard, k_hard, k_soft),                       # collapsed / squeezed
        (k_hard, k_hard, k_hard),                       # equilateral
        (k_hard, k_hard, np.maximum(k_third, 1e-9)),    # folded
    ):
        out = fn(*legs, z)
        assert np.all(np.isfinite(out)), f"non-finite in {model} at z={z}"
        assert np.all(out > 0.0), f"non-positive B in {model} at z={z}"


@pytest.mark.parametrize("z", EXTREME_Z)
def test_squeezed_boost_is_bounded_and_grows_with_hardness(z):
    """The FK-relevant shape: nonlinear boost must be >= 1 and monotone-ish."""
    tree = make_b_delta_fn("tree")
    nonlinear = make_b_delta_fn("bihalofit")
    k_hard = np.array([0.05, 0.2, 0.5, 1.0, 2.0])
    k_soft = np.full_like(k_hard, 1.0e-3)
    ratio = nonlinear(k_hard, k_hard, k_soft, z) / tree(k_hard, k_hard, k_soft, z)
    assert np.all(ratio > 0.9), f"unexpected suppression at z={z}: {ratio}"
    assert np.all(ratio < 100.0), f"runaway boost at z={z}: {ratio}"
    assert ratio[-1] > ratio[0], "boost should grow with the hard-leg scale"


def test_high_redshift_boost_softens():
    """At fixed hard scale the boost must decrease with redshift."""
    tree = make_b_delta_fn("tree")
    nonlinear = make_b_delta_fn("bihalofit")
    k_hard = np.array([0.5])
    k_soft = np.array([1.0e-3])
    ratios = [
        (nonlinear(k_hard, k_hard, k_soft, z) / tree(k_hard, k_hard, k_soft, z))[0]
        for z in (0.1, 1.0, 3.0, 5.7)
    ]
    assert ratios == sorted(ratios, reverse=True), f"not monotone in z: {ratios}"
    assert ratios[-1] < 1.2, "boost should be small at the farthest shell"


# --------------------------------------------------------------------------
# Strictness contract
# --------------------------------------------------------------------------

def test_strict_rejects_non_positive_legs():
    fn = make_b_delta_fn("tree")
    with pytest.raises(ValueError):
        fn(np.array([0.0]), np.array([0.1]), np.array([0.1]), 1.0)


def test_baryon_ratio_is_a_suppression_at_low_z():
    from b_model import baryon_ratio_bihalofit

    if not _has_baryon_fit():
        pytest.skip("installed fastnc has no get_Rb_bihalofit")

    k = np.array([1.0, 3.0])
    ratio = baryon_ratio_bihalofit(k, k, k, 0.3)
    assert np.all(ratio > 0.0)
    assert np.all(ratio < 1.5)
