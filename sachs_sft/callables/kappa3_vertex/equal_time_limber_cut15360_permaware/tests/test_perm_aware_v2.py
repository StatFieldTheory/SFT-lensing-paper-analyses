"""Tests of version 2 of the permutation-aware kappa3 callable (paper task T-001).

Two changes against version 1 are tested:

* K_111 (three Re Psi_0 legs) is reconstructed exactly from the channels for
  every leg order; version 1 is wrong there and only there, and both versions
  give the same sum over the six leg orders, which is what sft-wick folds.
* a query that is not a tabulated row raises ValueError unless
  ALLOW_INTERPOLATION is set; the production FK placements are tabulated rows.
"""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

_HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_HERE))
import perm_aware_kappa3_callable as v1  # noqa: E402

_PKG = Path(__file__).resolve().parents[5]
R1_TABLE = (_PKG / "sachs_sft/callables/kappa3_vertex/rebuild/products/pieces_nphi512"
            / "table_permclosed_np512.npz")
R1_MAIN_CONFIG = (_PKG / "sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24_corrected_main"
                  / "main/config_main.yaml")
ARCMIN = math.pi / 180.0 / 60.0
PERMS = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))


def _fresh_v2():
    """A newly imported version 2 module with an empty cache and default settings."""
    sys.modules.pop("perm_aware_kappa3_callable_v2", None)
    return importlib.import_module("perm_aware_kappa3_callable_v2")


# --- K_111: synthetic cumulant through the channel map ---------------------------

def _channels(k: np.ndarray) -> dict[str, float]:
    """Table channels of a real cumulant K_abc over (Phi00, Re Psi0, Im Psi0).

    Psi0 = Re + i Im on every spin leg; Dmod conjugates the third leg
    (canoes kappa3.py channel definitions).
    """
    return {
        "zeta_TTT": k[0, 0, 0],
        "zeta_TTP": k[0, 0, 1],
        "zeta_TPP": k[0, 1, 1] - k[0, 2, 2],
        "zeta_Bmod": k[0, 1, 1] + k[0, 2, 2],
        "zeta_PPP": k[1, 1, 1] - k[1, 2, 2] - k[2, 1, 2] - k[2, 2, 1],
        "zeta_Dmod": k[1, 1, 1] + k[1, 2, 2] + k[2, 1, 2] - k[2, 2, 1],
    }


def _random_cumulant(seed: int) -> np.ndarray:
    """Random K with the parity zeros and with K_122, K_212, K_221 all different."""
    k = np.random.default_rng(seed).normal(size=(3, 3, 3))
    for idx in np.ndindex(3, 3, 3):
        if sum(1 for x in idx if x == 2) % 2 == 1:
            k[idx] = 0.0
    assert len({k[1, 2, 2], k[2, 1, 2], k[2, 2, 1]}) == 3
    return k


def _reconstruct(module, k: np.ndarray) -> np.ndarray:
    """The callable's tensor for a query whose true cumulant is k."""
    def channel_at(channel: str, perm) -> np.ndarray:
        # new leg i = old leg perm[i]  <=>  np.transpose(k, perm)
        return np.array([_channels(np.transpose(k, perm))[channel]])
    return module._tensor_from_channels(channel_at, 1)[0]


@pytest.mark.parametrize("seed", range(5))
def test_v2_reconstructs_every_entry_for_every_leg_order(seed):
    v2 = _fresh_v2()
    k = _random_cumulant(seed)
    for sigma in PERMS:
        k_sigma = np.transpose(k, sigma)
        np.testing.assert_allclose(_reconstruct(v2, k_sigma), k_sigma, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("seed", range(5))
def test_v1_is_wrong_at_k111_only(seed):
    k = _random_cumulant(seed)
    for sigma in PERMS:
        k_sigma = np.transpose(k, sigma)
        got = _reconstruct(v1, k_sigma)
        mask = np.ones((3, 3, 3), dtype=bool)
        mask[1, 1, 1] = False
        np.testing.assert_allclose(got[mask], k_sigma[mask], rtol=1e-12, atol=1e-12)
        expected_v1 = (k_sigma[1, 1, 1] + 0.5 * (k_sigma[1, 2, 2] + k_sigma[2, 1, 2])
                       - k_sigma[2, 2, 1])
        assert got[1, 1, 1] == pytest.approx(expected_v1, rel=1e-12)
        assert abs(got[1, 1, 1] - k_sigma[1, 1, 1]) > 1e-6


@pytest.mark.parametrize("seed", range(5))
def test_six_order_sums_agree(seed):
    v2 = _fresh_v2()
    k = _random_cumulant(seed)
    s1 = sum(_reconstruct(v1, np.transpose(k, s)) for s in PERMS)
    s2 = sum(_reconstruct(v2, np.transpose(k, s)) for s in PERMS)
    np.testing.assert_allclose(s2, s1, rtol=1e-14, atol=1e-14 * np.max(np.abs(s1)))


# --- strict lookup on the production R1 table ------------------------------------

@pytest.fixture()
def strict_v2():
    if not R1_TABLE.exists():
        pytest.skip(f"R1 table not found at {R1_TABLE}")
    v2 = _fresh_v2()
    v2.TABLE_PATH = R1_TABLE
    return v2


def _call(module, n_list, lam: float) -> np.ndarray:
    n_arr = np.asarray(n_list, dtype=float).reshape(3, 1, 3)
    return module.coupling_fn_batch(n_arr, np.full((3, 1), lam))


def _shell(module) -> float:
    cache = module._ensure_cache()
    return float(cache["lambda"][cache["lambda"].size // 2])


def _equilateral(side_arcmin: float) -> list[np.ndarray]:
    """Three unit vectors with all pairwise separations equal to side_arcmin."""
    s = side_arcmin * ARCMIN
    rho = s / math.sqrt(3.0)  # small-angle circumradius
    return [np.array([rho * math.cos(a), rho * math.sin(a), 1.0])
            for a in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)]


@pytest.mark.parametrize("side", (2.0, 5.0))
def test_strict_lookup_rejects_an_equilateral_triangle(strict_v2, side):
    lam = _shell(strict_v2)
    with pytest.raises(ValueError, match="not tabulated rows"):
        _call(strict_v2, _equilateral(side), lam)


def test_strict_lookup_rejects_an_off_grid_collapsed_triple(strict_v2):
    lam = _shell(strict_v2)
    g = 1.0 * ARCMIN
    x, y = np.array([0.0, 0.0, 1.0]), np.array([math.sin(g), 0.0, math.cos(g)])
    with pytest.raises(ValueError, match="not tabulated rows"):
        _call(strict_v2, [x, x, y], lam)


def test_opt_in_blends_without_raising(strict_v2):
    strict_v2.ALLOW_INTERPOLATION = True
    out = _call(strict_v2, _equilateral(2.0), _shell(strict_v2))
    assert np.all(np.isfinite(out)) and np.any(out != 0.0)


def test_production_fk_placements_are_tabulated_rows(strict_v2):
    """Every (x, x, y) and (x, y, y) placement of the R1 FK sweep is an exact row."""
    cfg = yaml.safe_load(R1_MAIN_CONFIG.read_text())
    grid = cfg["sweep"]["positions_grid"]
    x = np.asarray(grid["x"][0], dtype=float)
    ys = [np.asarray(y, dtype=float) for y in grid["y"]]
    assert len(ys) == 40
    lam = _shell(strict_v2)
    for y in ys:
        for legs in ([x, x, y], [x, y, x], [y, x, x], [x, y, y], [y, x, y], [y, y, x]):
            out = _call(strict_v2, legs, lam)  # raises if any ordering is off-grid
            assert np.all(np.isfinite(out))


def test_v2_equals_v1_on_tabulated_rows_except_k111(strict_v2):
    """On exact rows only the K_111 entry can differ between the versions."""
    v1_mod = importlib.reload(v1)
    v1_mod.TABLE_PATH = R1_TABLE
    cfg = yaml.safe_load(R1_MAIN_CONFIG.read_text())
    grid = cfg["sweep"]["positions_grid"]
    x = np.asarray(grid["x"][0], dtype=float)
    lam = _shell(strict_v2)
    for y in (np.asarray(grid["y"][i], dtype=float) for i in (0, 10, 20, 30)):
        a = _call(v1_mod, [x, x, y], lam)[0]
        b = _call(strict_v2, [x, x, y], lam)[0]
        mask = np.ones((3, 3, 3), dtype=bool)
        mask[1, 1, 1] = False
        np.testing.assert_array_equal(a[mask], b[mask])
