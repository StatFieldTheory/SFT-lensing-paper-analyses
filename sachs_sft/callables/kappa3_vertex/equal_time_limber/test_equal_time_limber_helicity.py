"""Boundary tests for the spin-2 helicity reconstruction in the
equal_time_limber kappa3 callable.

Verifies the LOCAL-real-basis reconstruction the callable applies
(derive_kappa3_spin2_helicity.wl, 14/14 PASS):

    K011 = (zeta_TPP + Bmod)/2 ; K022 = (Bmod - zeta_TPP)/2
    K111 = (zeta_PPP + 3*Dmod)/4 ; K122 = (Dmod - zeta_PPP)/4
    odd-index-2 entries (K002,K012,K112,K222 + perms) == 0 (parity).

These are deterministic algebra tests against a synthetic 6-channel table; they
do NOT require the heavy production build.  Run with the sft-wick env:
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -m pytest \
        test_equal_time_limber_helicity.py -q
"""
from __future__ import annotations

import importlib
import json
from pathlib import Path

import numpy as np
import pytest

import equal_time_limber_kappa3_callable as cb


def _write_table(path: Path, *, with_modulus: bool) -> None:
    """A tiny synthetic table: 3 canonicalized cosine triples x 2 shells."""
    triples = np.array([[1.0, 1.0, 1.0], [0.5, 0.3, 0.2], [-0.5, 0.4, 0.1]])
    lam = np.array([1000.0, 2000.0])
    rng = np.random.default_rng(0)
    shp = (triples.shape[0], lam.size)
    data = {
        "cosine_triples": triples,
        "lambda_shells_Mpc": lam,
        "zeta_TTT": rng.standard_normal(shp),
        "zeta_TTP": rng.standard_normal(shp),
        "zeta_TPP": rng.standard_normal(shp),
        "zeta_PPP": rng.standard_normal(shp),
    }
    meta = {
        "already_R_contracted": False, "equal_time": True,
        "radial_density_measure": "lambda",
    }
    if with_modulus:
        data["zeta_Bmod"] = rng.standard_normal(shp)
        data["zeta_Dmod"] = rng.standard_normal(shp)
        data["has_modulus_channels_flag"] = np.bool_(True)
        meta["has_modulus_channels"] = True
    data["cosmo_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    np.savez(path, **data)


@pytest.fixture()
def mod_table(tmp_path):
    p = tmp_path / "synthetic_mod.npz"
    _write_table(p, with_modulus=True)
    cb.TABLE_PATH = p
    cb.coupling_fn.table_path = str(p)
    cb.reset_cache()
    yield p
    cb.reset_cache()


def _tensor_at_squeezed():
    """coupling_fn at the (1,1,1) cosine triple, shell 1000."""
    na = np.array([0.0, 0.0, 1.0])
    return cb.coupling_fn([na, na, na], [1000.0, 1000.0, 1000.0])


def test_loads_six_channels(mod_table):
    cache = cb._ensure_cache()
    assert set(("zeta_Bmod", "zeta_Dmod")).issubset(cache["per_channel"].keys())


def test_parity_zero_entries(mod_table):
    K = _tensor_at_squeezed()
    # every entry with an ODD number of index-2 legs must be exactly 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                n2 = (i == 2) + (j == 2) + (k == 2)
                if n2 % 2 == 1:
                    assert K[i, j, k] == 0.0, f"K[{i},{j},{k}] not parity-zero"


def test_reconstruction_identities(mod_table):
    """K011+K022 == Bmod ; K011-K022 == zeta_TPP ; K111,K122 vs PPP/Dmod."""
    cache = cb._ensure_cache()
    # exact-hit at the (1,1,1) triple, shell index 0 (lam=1000)
    g = cb._geom(1.0, 1.0, 1.0, 1000.0)
    TPP = cb._chan("zeta_TPP", g); PPP = cb._chan("zeta_PPP", g)
    Bmod = cb._chan("zeta_Bmod", g); Dmod = cb._chan("zeta_Dmod", g)
    K = _tensor_at_squeezed()
    assert K[0, 1, 1] == pytest.approx(0.5 * (TPP + Bmod), rel=1e-12)
    assert K[0, 2, 2] == pytest.approx(0.5 * (Bmod - TPP), rel=1e-12)
    assert K[1, 1, 1] == pytest.approx(0.25 * (PPP + 3.0 * Dmod), rel=1e-12)
    assert K[1, 2, 2] == pytest.approx(0.25 * (Dmod - PPP), rel=1e-12)
    # F-vertex SUM identity (the Sachs sigma-sigmabar vertex needs K011+K022):
    assert K[0, 1, 1] + K[0, 2, 2] == pytest.approx(Bmod, rel=1e-12)
    # helicity DIFFERENCE the old code used:
    assert K[0, 1, 1] - K[0, 2, 2] == pytest.approx(TPP, rel=1e-12)


def test_permutation_symmetry(mod_table):
    K = _tensor_at_squeezed()
    assert np.allclose([K[0, 1, 1], K[1, 0, 1], K[1, 1, 0]], K[0, 1, 1])
    assert np.allclose([K[0, 2, 2], K[2, 0, 2], K[2, 2, 0]], K[0, 2, 2])
    assert np.allclose([K[1, 2, 2], K[2, 1, 2], K[2, 2, 1]], K[1, 2, 2])


def test_scalar_equals_batch(mod_table):
    na = np.array([0.0, 0.0, 1.0])
    nb = np.array([np.sin(0.3), 0.0, np.cos(0.3)])
    nc = np.array([np.sin(0.7), 0.0, np.cos(0.7)])
    for nl in ([na, na, nb], [na, nb, nc]):
        tl = [1500.0, 1500.0, 1500.0]
        Ks = cb.coupling_fn(nl, tl)
        n_arr = np.stack([np.asarray(x)[None, :] for x in nl], axis=0)
        t_arr = np.array([[tl[0]], [tl[1]], [tl[2]]])
        Kb = cb.coupling_fn_batch(n_arr, t_arr)[0]
        assert np.array_equal(Ks, Kb), "scalar and batch reconstruction differ"


def test_legacy_table_raises(tmp_path):
    """A 4-channel (pre-fix) table must FAIL LOUD, not silently corrupt."""
    p = tmp_path / "legacy.npz"
    _write_table(p, with_modulus=False)
    cb.TABLE_PATH = p
    cb.reset_cache()
    with pytest.raises(ValueError, match="modulus channels"):
        cb._ensure_cache()
    cb.reset_cache()
