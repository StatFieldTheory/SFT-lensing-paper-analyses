"""kappa3_vertex / equal_time_limber — sft-wick non-local κ³ vertex callable.

The SECOND κ³ implementation: NON-R-contracted equal-time Limber κ³ (LOW exact
Wigner-3j ℓ≤60 + HIGH flat-sky Born-Limber ℓ>60), the source of the validated
FK κκ ``+3.086e-5`` baseline.

Exposes ``coupling_fn(n_list, t_list) -> (3, 3, 3)`` (scalar) and
``coupling_fn_batch(n_arr, t_arr) -> (n_samples, 3, 3, 3)`` (``.vectorized=True``).

Contrast with R_contracted:
  * ``already_R_contracted = False`` — sft-wick applies the response R itself.
  * indexed by COSINE-triples (cKDTree, cos γ_ij = n_i·n_j) + λ_shells, NOT
    (L,L,ℓ₃)+λ.  Each channel's ζ is interpolated in (cos-triple, λ_mean).
  * lambda-density (not chi-density).

Component order (Φ00, Re Ψ0, Im Ψ0).  The (3,3,3) tensor is filled with the
LOCAL real-basis cumulant K_{abc} reconstructed from the all-+2 helicity
channels (zeta_TTT/TTP/TPP/PPP) AND the conjugate-helicity "modulus" channels
(zeta_Bmod=⟨Φ Ψ0 Ψ0*⟩ spin (0,2,-2); zeta_Dmod=⟨Ψ0 Ψ0 Ψ0*⟩ spin (2,2,-2)):
    K000=TTT; K001(+perm)=TTP; K011(+perm)=(TPP+Bmod)/2; K022(+perm)=(Bmod-TPP)/2;
    K111=(PPP+3·Dmod)/4; K122(+perm)=(Dmod-PPP)/4.
The index-2 (σ_×) entries appear ONLY in pairs (K022, K122); entries with an ODD
count of index-2 (K002,K012,K112,K222 + perms) stay parity-zero (np.zeros init).
This is the verified spin-2 fix: see scripts/mathematica/derive_kappa3_spin2_helicity.wl
(14/14 PASS).  The modulus channels carry the d^L_{2,2}~O(1) piece that dominates
the γ^4-suppressed all-+2 channels at small angle.

Design (sachs_sft traceable callable): explicit ``TABLE_PATH`` (ONE combined
LOW+HIGH NPZ; no ``SFT_WICK_*`` env, no runtime LOW/HIGH summation), conventions
read from the NPZ ``cosmo_meta`` (never hand-typed), fail-loud guards.

Provenance — table built by ``build_equal_time_limber_table.py`` (reproducible).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Final

import numpy as np
from numpy.typing import NDArray

def _resolve_canoes_root() -> Path:
    """Locate the canoes checkout, tolerating repository moves.

    The repository was renamed (projects/canoes -> projects/AngStats ->
    projects/angular_statistics/canoes), so a single hardcoded path rots.
    Order: CANOES_ROOT env var, an importable canoes, then known locations.
    """
    import os

    env = os.environ.get("CANOES_ROOT")
    if env:
        return Path(env)
    try:
        import canoes  # already importable: use its own location
        return Path(canoes.__file__).resolve().parents[2]
    except Exception:
        pass
    for candidate in (
        Path("/Users/zzhang/projects/angular_statistics/canoes"),
        Path("/Users/zzhang/projects/canoes"),
    ):
        if (candidate / "src" / "canoes").is_dir() or (candidate / "canoes").is_dir():
            return candidate
    raise RuntimeError(
        "canoes checkout not found; set CANOES_ROOT or pip install -e it"
    )


_CANOES_ROOT = _resolve_canoes_root()
for _p in (_CANOES_ROOT / "src", _CANOES_ROOT):
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_HERE = Path(__file__).resolve().parent
TABLE_PATH = Path(r"/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/products/table_bihalofit_cut7680_r4.npz")

# --- sft-wick contract --------------------------------------------------------
QUERY_COORDINATE = "lambda_project_mpc"
QUERY_UNITS = "physical Mpc"
COMPONENT_ORDER = ("Phi_00", "Re Psi_0", "Im Psi_0")
ALREADY_R_CONTRACTED = False   # sft-wick applies the response R itself
EQUAL_TIME = True
RADIAL_MEASURE = "lambda"      # the +3.08e-5 baseline is lambda-density

_CHANNELS: Final = ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP",
                    "zeta_Bmod", "zeta_Dmod")
_CACHE: dict | None = None


def _read_meta(d) -> dict:
    if "cosmo_meta" not in d.files:
        return {}
    try:
        return json.loads(bytes(d["cosmo_meta"]).decode())
    except Exception:
        return {}


def _ensure_cache() -> dict:
    """Load the single combined LOW+HIGH ζ table + build the cKDTree (once)."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    from scipy.spatial import cKDTree

    if not TABLE_PATH.exists():
        raise FileNotFoundError(
            f"equal_time_limber κ³ table not found: {TABLE_PATH}\n"
            f"Build it with: python build_equal_time_limber_table.py."
        )
    d = np.load(TABLE_PATH, allow_pickle=False)
    meta = _read_meta(d)

    # Fail-loud guards (read from metadata; never hand-typed).
    if meta.get("already_R_contracted") is not False:
        raise ValueError(
            "equal_time_limber table must have already_R_contracted=False "
            "(sft-wick applies the response R itself). "
            f"Got {meta.get('already_R_contracted')!r}."
        )
    if meta.get("equal_time") is not True:
        raise ValueError(
            f"equal_time_limber table must have equal_time=True; got "
            f"{meta.get('equal_time')!r}."
        )
    dens = meta.get("radial_density_measure") or meta.get("radial_measure")
    if dens != RADIAL_MEASURE:
        raise ValueError(
            f"equal_time_limber table density convention is {dens!r}, expected "
            f"{RADIAL_MEASURE!r} (the +3.08e-5 baseline is lambda-density). "
            f"A chi-density table would bias the K-vertex by (dχ/dλ)^3 per shell."
        )
    # Fail-loud: the spin-2 reconstruction needs the conjugate-helicity modulus
    # channels (zeta_Bmod/zeta_Dmod). A legacy 4-channel table would make
    # Bmod=Dmod=0 and inject a spurious K022=-zeta_TPP/2 into the σ_× slots, so
    # we refuse it rather than silently corrupt the B-mode sector.
    has_mod = bool(meta.get("has_modulus_channels", False))
    if "has_modulus_channels_flag" in d.files:
        has_mod = has_mod or bool(d["has_modulus_channels_flag"])
    if not has_mod:
        raise ValueError(
            "equal_time_limber table lacks the spin-2 modulus channels "
            "(zeta_Bmod/zeta_Dmod; meta has_modulus_channels not set). A legacy "
            "4-channel table would corrupt the σ_× (index-2) K-vertex slots. "
            "Rebuild with build_equal_time_limber_table.py (which stamps "
            "has_modulus_channels=True and writes zeta_Bmod/zeta_Dmod)."
        )

    cos_triples = np.asarray(d["cosine_triples"], dtype=np.float64)
    lam = np.asarray(d["lambda_shells_Mpc"], dtype=np.float64)
    # Canonicalize the (descending-sorted) cosine triples, averaging duplicates.
    sorted_triples = np.sort(cos_triples, axis=1)[:, ::-1]
    uniq, inverse = np.unique(np.round(sorted_triples, 8), axis=0, return_inverse=True)
    n_uniq = int(uniq.shape[0])
    counts = np.bincount(inverse, minlength=n_uniq).astype(np.float64)
    per_channel: dict[str, NDArray[np.float64]] = {}
    for ch in _CHANNELS:
        if ch not in d.files:
            per_channel[ch] = np.zeros((n_uniq, lam.size))
            continue
        raw = np.nan_to_num(np.asarray(d[ch], dtype=np.float64), nan=0.0)
        avg = np.zeros((n_uniq, lam.size), dtype=np.float64)
        np.add.at(avg, inverse, raw)
        avg /= counts[:, None]
        per_channel[ch] = avg

    nonzero = [ch for ch in _CHANNELS if np.any(per_channel[ch] != 0)]
    _CACHE = {
        "uniq": uniq, "per_channel": per_channel, "tree": cKDTree(uniq),
        "lambda": lam, "lam_min": float(lam[0]), "lam_max": float(lam[-1]),
        "n_uniq": n_uniq, "meta": meta,
    }
    print(f"[etl-callable] loaded {len(nonzero)}/{len(_CHANNELS)} nonzero channels "
          f"{tuple(nonzero)} (U={n_uniq}, N_lam={lam.size}) from {TABLE_PATH.name}")
    return _CACHE


def _geom(cos12: float, cos23: float, cos31: float, lam: float):
    cache = _ensure_cache()
    if lam < cache["lam_min"] or lam > cache["lam_max"]:
        return None
    a, b, c = sorted((cos12, cos23, cos31), reverse=True)
    k = min(4, cache["n_uniq"])
    dist, idx = cache["tree"].query(np.array((a, b, c)), k=k)
    dist = np.atleast_1d(dist); idx = np.atleast_1d(idx)
    lam_grid = cache["lambda"]
    i_hi = min(int(np.searchsorted(lam_grid, lam, side="right")), lam_grid.size - 1)
    i_lo = max(i_hi - 1, 0)
    span = lam_grid[i_hi] - lam_grid[i_lo]
    w_lam = 0.0 if span <= 0 else float((lam - lam_grid[i_lo]) / span)
    return idx, 1.0 / np.maximum(dist, 1e-12), bool(dist[0] < 1e-10), i_lo, i_hi, w_lam


def _chan(channel: str, geom) -> float:
    if geom is None:
        return 0.0
    idx, w, exact, i_lo, i_hi, w_lam = geom
    z = _ensure_cache()["per_channel"][channel]
    z_at = (1.0 - w_lam) * z[idx, i_lo] + w_lam * z[idx, i_hi]
    return float(z_at[0]) if exact else float(np.sum(w * z_at) / np.sum(w))


def coupling_fn(n_list, t_list) -> NDArray[np.float64]:
    """Equal-time Limber κ³, symmetric (3, 3, 3) tensor. t_list = λ_project [Mpc]."""
    n = np.asarray(n_list, dtype=float)
    n = n / np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-300)
    c12 = float(np.clip(np.dot(n[0], n[1]), -1.0, 1.0))
    c23 = float(np.clip(np.dot(n[1], n[2]), -1.0, 1.0))
    c31 = float(np.clip(np.dot(n[2], n[0]), -1.0, 1.0))
    g = _geom(c12, c23, c31, float(np.mean(np.asarray(t_list, dtype=float))))
    v = {ch: _chan(ch, g) for ch in _CHANNELS}
    # Local real-basis reconstruction (derive_kappa3_spin2_helicity.wl, PASS).
    K011 = 0.5 * (v["zeta_TPP"] + v["zeta_Bmod"])
    K022 = 0.5 * (v["zeta_Bmod"] - v["zeta_TPP"])
    K111 = 0.25 * (v["zeta_PPP"] + 3.0 * v["zeta_Dmod"])
    K122 = 0.25 * (v["zeta_Dmod"] - v["zeta_PPP"])
    out = np.zeros((3, 3, 3), dtype=float)
    out[0, 0, 0] = v["zeta_TTT"]
    for a, b, c in ((0, 0, 1), (0, 1, 0), (1, 0, 0)):
        out[a, b, c] = v["zeta_TTP"]
    for a, b, c in ((0, 1, 1), (1, 0, 1), (1, 1, 0)):
        out[a, b, c] = K011
    for a, b, c in ((0, 2, 2), (2, 0, 2), (2, 2, 0)):
        out[a, b, c] = K022
    out[1, 1, 1] = K111
    for a, b, c in ((1, 2, 2), (2, 1, 2), (2, 2, 1)):
        out[a, b, c] = K122
    # Odd-index-2 entries (K002,K012,K112,K222 + perms) stay 0 by parity.
    return out


def coupling_fn_batch(n_arr, t_arr) -> NDArray[np.float64]:
    """Vectorized equal-time Limber κ³ (sft-wick ``vectorized=True`` fast path).

    n_arr: (3, n_samples, 3); t_arr: (3, n_samples). Returns (n_samples,3,3,3).
    """
    n_arr = np.asarray(n_arr, dtype=float)
    t_arr = np.asarray(t_arr, dtype=float)
    if n_arr.ndim != 3 or n_arr.shape[0] != 3 or n_arr.shape[2] != 3:
        raise ValueError(f"n_arr must be (3, n_samples, 3); got {n_arr.shape}")
    if t_arr.ndim != 2 or t_arr.shape[0] != 3 or t_arr.shape[1] != n_arr.shape[1]:
        raise ValueError(f"t_arr must be (3, n_samples); got {t_arr.shape}")
    ns = n_arr.shape[1]
    nn = n_arr / np.maximum(np.linalg.norm(n_arr, axis=-1, keepdims=True), 1e-300)
    n0, n1, n2 = nn[0], nn[1], nn[2]
    cos12 = np.clip(np.einsum("ij,ij->i", n0, n1), -1.0, 1.0)
    cos23 = np.clip(np.einsum("ij,ij->i", n1, n2), -1.0, 1.0)
    cos31 = np.clip(np.einsum("ij,ij->i", n2, n0), -1.0, 1.0)
    lam_mean = np.mean(t_arr, axis=0)
    queries = np.sort(np.stack([cos12, cos23, cos31], axis=-1), axis=1)[:, ::-1]

    cache = _ensure_cache()
    lam_grid = cache["lambda"]
    in_range = (lam_mean >= cache["lam_min"]) & (lam_mean <= cache["lam_max"])
    out = np.zeros((ns, 3, 3, 3), dtype=float)
    valid = np.flatnonzero(in_range)
    if valid.size == 0:
        return out
    k = min(4, cache["n_uniq"])
    dist, nn_idx = cache["tree"].query(queries[valid], k=k)
    if k == 1:
        dist = dist[:, None]; nn_idx = nn_idx[:, None]
    lam_v = lam_mean[valid]
    i_hi = np.minimum(np.searchsorted(lam_grid, lam_v, side="right"), lam_grid.size - 1)
    i_lo = np.maximum(i_hi - 1, 0)
    span = lam_grid[i_hi] - lam_grid[i_lo]
    w_lam = np.where(span > 0, (lam_v - lam_grid[i_lo]) / span, 0.0)
    w = 1.0 / np.maximum(dist, 1e-12)
    w_sum = np.sum(w, axis=1)
    exact = dist[:, 0] < 1e-10

    def _lookup(ch: str) -> NDArray[np.float64]:
        z = cache["per_channel"][ch]
        z_at = (1.0 - w_lam[:, None]) * z[nn_idx, i_lo[:, None]] + \
               w_lam[:, None] * z[nn_idx, i_hi[:, None]]
        return np.where(exact, z_at[:, 0], np.sum(w * z_at, axis=1) / w_sum)

    v = {ch: np.zeros(ns) for ch in _CHANNELS}
    for ch in _CHANNELS:
        v[ch][valid] = _lookup(ch)
    # Local real-basis reconstruction (must match coupling_fn exactly).
    K011 = 0.5 * (v["zeta_TPP"] + v["zeta_Bmod"])
    K022 = 0.5 * (v["zeta_Bmod"] - v["zeta_TPP"])
    K111 = 0.25 * (v["zeta_PPP"] + 3.0 * v["zeta_Dmod"])
    K122 = 0.25 * (v["zeta_Dmod"] - v["zeta_PPP"])
    out[:, 0, 0, 0] = v["zeta_TTT"]
    for a, b, c in ((0, 0, 1), (0, 1, 0), (1, 0, 0)):
        out[:, a, b, c] = v["zeta_TTP"]
    for a, b, c in ((0, 1, 1), (1, 0, 1), (1, 1, 0)):
        out[:, a, b, c] = K011
    for a, b, c in ((0, 2, 2), (2, 0, 2), (2, 2, 0)):
        out[:, a, b, c] = K022
    out[:, 1, 1, 1] = K111
    for a, b, c in ((1, 2, 2), (2, 1, 2), (2, 2, 1)):
        out[:, a, b, c] = K122
    # Odd-index-2 entries stay 0 by parity (np.zeros init).
    return out


coupling_fn_batch.vectorized = True  # type: ignore[attr-defined]
coupling_fn.query_coordinate = QUERY_COORDINATE  # type: ignore[attr-defined]
coupling_fn.component_order = COMPONENT_ORDER    # type: ignore[attr-defined]
coupling_fn.already_R_contracted = ALREADY_R_CONTRACTED  # type: ignore[attr-defined]
coupling_fn.equal_time = EQUAL_TIME              # type: ignore[attr-defined]
coupling_fn.radial_measure = RADIAL_MEASURE      # type: ignore[attr-defined]
coupling_fn.table_path = str(TABLE_PATH)         # type: ignore[attr-defined]


def reset_cache() -> None:
    global _CACHE
    _CACHE = None


__all__ = ["coupling_fn", "coupling_fn_batch", "reset_cache",
           "TABLE_PATH", "QUERY_COORDINATE", "COMPONENT_ORDER"]


if __name__ == "__main__":  # self-check
    c = _ensure_cache()
    print(f"[equal_time_limber coupling_fn] table: {TABLE_PATH.name}")
    print(f"  already_R_contracted={ALREADY_R_CONTRACTED} equal_time={EQUAL_TIME} "
          f"measure={RADIAL_MEASURE}")
    print(f"  cosmology={c['meta'].get('cosmology')} "
          f"λ-grid n={c['lambda'].size} [{c['lam_min']:.1f}..{c['lam_max']:.1f}] Mpc")
    na = np.array([0.0, 0.0, 1.0]); th = 0.3
    nb = np.array([np.sin(th), 0.0, np.cos(th)])
    t = float(c["lambda"][len(c["lambda"]) // 2])
    K = coupling_fn([na, na, nb], [t, t, t])
    print(f"  coupling_fn(squeezed, t={t:.1f}) shape={K.shape}, "
          f"TTT={K[0,0,0]:.3e}, max|K|={np.max(np.abs(K)):.3e}")
