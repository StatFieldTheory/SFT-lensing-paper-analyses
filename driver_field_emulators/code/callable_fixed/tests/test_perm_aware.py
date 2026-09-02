"""Correctness tests for the permutation-aware kappa3 callable.

The test that matters is not relabelling invariance: the production
callable satisfies that trivially, because it symmetrises. The test that
discriminates is whether each tensor entry equals the tabulated channel at
the geometry that entry's own field assignment requires.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perm_aware_kappa3_callable as pa  # noqa: E402

DEPLOYED = (Path(__file__).resolve().parents[4]
            / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
            / "kappa3_vertex" / "equal_time_limber"
            / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz")


def _directions(cos01, cos12, cos20):
    """Three unit vectors realising a given cosine triple, or None."""
    n0 = np.array([0.0, 0.0, 1.0])
    t1 = np.arccos(np.clip(cos01, -1, 1))
    n1 = np.array([np.sin(t1), 0.0, np.cos(t1)])
    t2 = np.arccos(np.clip(cos20, -1, 1))
    # azimuth of n2 fixed by cos12
    denom = np.sin(t1) * np.sin(t2)
    if abs(denom) < 1e-12:
        return None
    cos_phi = (cos12 - np.cos(t1) * np.cos(t2)) / denom
    if not -1.0 <= cos_phi <= 1.0:
        return None
    phi = np.arccos(cos_phi)
    n2 = np.array([np.sin(t2) * np.cos(phi), np.sin(t2) * np.sin(phi), np.cos(t2)])
    return np.stack([n0, n1, n2])


@pytest.fixture(scope="module")
def deployed_table(tmp_path_factory):
    """Point the callable at the deployed table for these tests."""
    pa.TABLE_PATH = DEPLOYED
    pa._CACHE = None
    yield np.load(DEPLOYED, allow_pickle=False)
    pa._CACHE = None


def test_permute_triple_matches_direct_cosines():
    """The index shuffle must agree with recomputing the dot products."""
    rng = np.random.default_rng(0)
    n = rng.normal(size=(3, 3))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    triple = np.array([n[0] @ n[1], n[1] @ n[2], n[2] @ n[0]])
    for perm in pa._PERMS:
        a, b, c = perm
        direct = np.array([n[a] @ n[b], n[b] @ n[c], n[c] @ n[a]])
        shuffled = pa._permute_triple(triple, perm)
        assert np.allclose(shuffled, direct, atol=1e-14), perm


def test_scalar_entry_is_the_tabulated_TTT(deployed_table):
    """K000 must be zeta_TTT at the query's own geometry."""
    d = deployed_table
    row = 400
    triple = d["cosine_triples"][row]
    n = _directions(*triple)
    if n is None:
        pytest.skip("degenerate row")
    lam = float(d["lambda_shells_Mpc"][-1])
    out = pa.coupling_fn(n, [lam] * 3)
    assert out[0, 0, 0] == pytest.approx(d["zeta_TTT"][row, -1], rel=2e-3)


def test_single_spin_leg_uses_a_permuted_geometry(deployed_table):
    """K001 and K010 must differ, and each match zeta_TTP at its own geometry.

    The production callable sets them equal. On an asymmetric triangle that
    is wrong, and this is the entry where it is most visible.
    """
    d = deployed_table
    triples = d["cosine_triples"]
    key = {tuple(np.round(t, 12)): i for i, t in enumerate(triples)}
    lam_index = -1
    lam = float(d["lambda_shells_Mpc"][lam_index])
    checked = 0
    for row, triple in enumerate(triples):
        if len(set(np.round(triple, 9))) < 3:
            continue                       # need a genuinely asymmetric row
        n = _directions(*triple)
        if n is None:
            continue
        out = pa.coupling_fn(n, [lam] * 3)
        # K001: spin on leg 2, already canonical
        want_001 = d["zeta_TTP"][row, lam_index]
        # K010: spin on leg 1, canonical ordering is (0, 2, 1)
        perm = (0, 2, 1)
        permuted = tuple(np.round(pa._permute_triple(triple, perm), 12))
        if permuted not in key:
            continue
        want_010 = d["zeta_TTP"][key[permuted], lam_index]
        assert out[0, 0, 1] == pytest.approx(want_001, rel=5e-3, abs=1e-30)
        assert out[0, 1, 0] == pytest.approx(want_010, rel=5e-3, abs=1e-30)
        checked += 1
        if checked >= 5:
            break
    assert checked >= 3, "not enough asymmetric rows exercised"


def test_symmetric_geometry_reproduces_the_production_tensor(deployed_table):
    """On an equilateral triangle the permutation cannot matter.

    Every ordering gives the same geometry, so the corrected callable must
    return exactly what the production one does. This is the check that the
    corrections are confined to asymmetric configurations.
    """
    sys.path.insert(0, str(DEPLOYED.parent))
    import equal_time_limber_kappa3_callable as prod

    prod.TABLE_PATH = DEPLOYED
    prod._CACHE = None
    cos = 0.42857142857142855                     # a mesh value, 3/7
    n = _directions(cos, cos, cos)
    assert n is not None
    lam = float(deployed_table["lambda_shells_Mpc"][-1])
    got = pa.coupling_fn(n, [lam] * 3)
    want = prod.coupling_fn(n, [lam] * 3)
    assert np.allclose(got, want, rtol=2e-3, atol=1e-30), np.abs(got - want).max()


def test_parity_zeros_are_preserved(deployed_table):
    """Entries with an odd number of Im slots must stay exactly zero."""
    rng = np.random.default_rng(3)
    n = rng.normal(size=(3, 3))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    lam = float(deployed_table["lambda_shells_Mpc"][-1])
    out = pa.coupling_fn(n, [lam] * 3)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                if sum(1 for x in (a, b, c) if x == 2) % 2 == 1:
                    assert out[a, b, c] == 0.0, (a, b, c)


def test_outside_the_shell_range_returns_zero(deployed_table):
    n = np.eye(3)
    out = pa.coupling_fn(n, [1.0] * 3)
    assert np.all(out == 0.0)


def test_the_defect_is_real_and_the_fix_removes_it(deployed_table):
    """The production callable sets K001 = K010; the table says otherwise.

    Without this test the fix could be a no-op. It pins the size of the
    disagreement on the deployed table, so a future change that quietly
    reintroduces the symmetrised fill fails here.
    """
    sys.path.insert(0, str(DEPLOYED.parent))
    import equal_time_limber_kappa3_callable as prod

    prod.TABLE_PATH = DEPLOYED
    prod._CACHE = None
    d = deployed_table
    triples = d["cosine_triples"]
    key = {tuple(np.round(t, 12)): i for i, t in enumerate(triples)}
    lam = float(d["lambda_shells_Mpc"][-1])

    spreads = []
    for row, triple in enumerate(triples):
        if len(set(np.round(triple, 9))) < 3:
            continue
        permuted = tuple(np.round(pa._permute_triple(triple, (0, 2, 1)), 12))
        if permuted not in key:
            continue
        n = _directions(*triple)
        if n is None:
            continue
        fixed = pa.coupling_fn(n, [lam] * 3)
        legacy = prod.coupling_fn(n, [lam] * 3)
        # the production tensor is symmetric by construction
        assert legacy[0, 0, 1] == pytest.approx(legacy[0, 1, 0], rel=1e-12)
        scale = max(abs(fixed[0, 0, 1]), abs(fixed[0, 1, 0]))
        if scale > 0:
            spreads.append(abs(fixed[0, 0, 1] - fixed[0, 1, 0]) / scale)
        if len(spreads) >= 40:
            break

    assert spreads, "no asymmetric rows exercised"
    median = float(np.median(spreads))
    assert median > 0.05, (
        f"the corrected callable barely distinguishes the two placements "
        f"(median spread {median:.3e}); the fix would be a no-op")
    print(f"\n  median |K001 - K010| / max on asymmetric rows: {median:.3f} "
          f"over {len(spreads)} rows (production forces this to 0)")
