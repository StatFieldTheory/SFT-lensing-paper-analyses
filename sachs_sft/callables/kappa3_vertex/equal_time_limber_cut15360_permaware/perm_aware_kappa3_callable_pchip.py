"""Candidate signed-log PCHIP radial interpolation for native kappa3 channels.

Pending matched actual-FK/native-channel reference and fold review; no adoption
claim. PCHIP interpolates log(abs(channel)) in physical lambda on contiguous
same-sign nonzero runs. Zero and sign-crossing intervals retain the linear rule.
Angular blending follows radial evaluation, as in the loglam candidate.

**Inherited version 2 changes (2026-10, paper task T-001).** Version 2 differs from
``perm_aware_kappa3_callable.py`` (version 1, whose sha256 is pinned by
``r1_sft061/permaware_ladder/guard.py`` and recorded in the R1 run records)
in two ways, retained by this candidate:

* the entry with three ``Re Psi_0`` legs, K_111, averages the conjugated leg
  of ``zeta_Dmod`` over the three legs:
  ``K_111 = (zeta_PPP + D[0] + D[1] + D[2]) / 4`` with ``D[k]`` the Dmod
  channel conjugated on leg k. Version 1 used ``(zeta_PPP + 3 D[2]) / 4``,
  which equals ``K_111 + (K_122 + K_212)/2 - K_221`` and is K_111 only when
  the three legs are equivalent. A fold that sums over the six leg
  assignments (sft-wick >= 0.5) gets the same result from both versions; a
  consumer that reads one leg order does not.
* a query that is not a tabulated row raises ``ValueError`` (strict lookup).
  Version 1 answered it silently with an inverse-distance blend of the four
  nearest rows. Set ``ALLOW_INTERPOLATION = True`` on the module to blend on
  purpose.

Version 1 corrected two defects in the production callable
(`sachs_sft/callables/kappa3_vertex/equal_time_limber/
equal_time_limber_kappa3_callable.py`), both measured in
`driver_field_emulators/note`.

**Defect 1, leg permutations.** The vertex is a function of three
(field, direction) pairs, and its tabulated channels fix which legs carry
the spin-2 field: `zeta_TTP` has spin weights (0, 0, 2), so the spin leg is
the third; `zeta_TPP` has (0, 2, 2), so the scalar leg is the first. The
production callable sorts the query triple and then writes one value into
every index permutation of the coupling tensor. That is exact for the fully
symmetric `zeta_TTT`, and wrong for the rest: asking for the spin field on
a different leg is a question about a PERMUTED geometry, not about the same
geometry with a relabelled channel. Measured cost at the collapsed
configurations: a median 23% on `zeta_TPP` and 70% on `zeta_Dmod`.

Here every tensor entry is evaluated at the geometry its own field
assignment requires. Concretely, with legs written 0,1,2:

* no spin leg: `zeta_TTT`, any ordering.
* one spin leg at p: permute p into slot 2 and read `zeta_TTP`.
* two spin legs, scalar at q: permute q into slot 0 and read the
  `(zeta_TPP, zeta_Bmod)` pair.
* three spin legs, all `Re`: `zeta_PPP` at the given ordering and
  `zeta_Dmod` at the three cyclic orderings, which put the conjugated leg on
  each leg in turn (version 2; see above).
* three spin legs, one `Re` and two `Im`: this is the case the table cannot
  answer with a single lookup. `(zeta_Dmod - zeta_PPP)/4` is the AVERAGE of
  two of the three placements of the `Re` leg, so three cyclically permuted
  evaluations are combined to isolate one placement. With a symmetric
  geometry the three are equal and the combination collapses to the
  production value, which is the consistency check in
  `tests/test_perm_aware.py`.

**Defect 2, the mesh key.** The production callable deduplicates triples
after rounding to 8 decimals, so two collapsed configurations closer than
`1 - cos g = 5e-9` become one node however finely the table samples. That is
a hard floor at `g = 0.344` arcmin, and the draft plots down to 0.5. Here
the key is 12 decimals, a floor at 0.0034 arcmin. The production mesh
deduplicates to the same 351 rows at either precision, so nothing else
changes.

The table must be closed under the leg permutations that get queried. The
production mesh already is, being a symmetric meshgrid; the densified
collapsed family is made so by
`code/rebuild/make_perm_closed_triples.py`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Final

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from scipy.interpolate import PchipInterpolator

_HERE = Path(__file__).resolve().parent
TABLE_PATH = _HERE / "table_permclosed.npz"

_CHANNELS: Final = ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP",
                    "zeta_Bmod", "zeta_Dmod")
#: Decimals used to canonicalise the mesh key. See "Defect 2" above.
_KEY_DECIMALS: Final = 12
#: All six orderings of the three legs.
_PERMS: Final = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
#: The three cyclic orderings; slot 2 of each holds a different original leg.
_CYCLIC: Final = ((0, 1, 2), (1, 2, 0), (2, 0, 1))
#: A query counts as a tabulated row when its nearest row is closer than this
#: in cosine space. The production FK placements sit at 9.7e-15 to 3.4e-11
#: (the margin comes from the 10 decimals of the config direction vectors).
_EXACT_TOL: Final = 1e-10
#: Opt-in for deliberate interpolation between rows; off means a query that is
#: not a tabulated row raises ValueError.
ALLOW_INTERPOLATION = False

_CACHE: dict | None = None


def _ensure_cache() -> dict:
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    d = np.load(TABLE_PATH, allow_pickle=False)
    meta = json.loads(np.asarray(d["cosmo_meta"], dtype=np.uint8).tobytes().decode())
    if not bool(d.get("has_modulus_channels_flag", False)):
        raise RuntimeError(f"{TABLE_PATH.name} lacks the modulus channels")
    if str(meta.get("radial_measure")) != "lambda":
        raise RuntimeError(f"unexpected radial measure {meta.get('radial_measure')!r}")

    triples = np.asarray(d["cosine_triples"], dtype=np.float64)
    lam = np.asarray(d["lambda_shells_Mpc"], dtype=np.float64)

    # Deduplicate WITHOUT sorting: the ordering of a row is physical here,
    # because it says which leg carries which field.
    key = np.round(triples, _KEY_DECIMALS)
    uniq, inverse = np.unique(key, axis=0, return_inverse=True)
    counts = np.bincount(inverse, minlength=uniq.shape[0]).astype(np.float64)
    per_channel = {}
    for ch in _CHANNELS:
        raw = np.nan_to_num(np.asarray(d[ch], dtype=np.float64), nan=0.0)
        acc = np.zeros((uniq.shape[0], lam.size), dtype=np.float64)
        np.add.at(acc, inverse, raw)
        per_channel[ch] = acc / counts[:, None]
        # Avoid summation changing the sign bit of zero on unique native rows.
        single = counts[inverse] == 1
        per_channel[ch][inverse[single]] = raw[single]

    _CACHE = {"uniq": uniq, "tree": cKDTree(uniq), "per_channel": per_channel,
              "lambda": lam, "lam_min": float(lam[0]), "lam_max": float(lam[-1]),
              "n_uniq": int(uniq.shape[0]), "meta": meta, "pchippers": {}}
    closure = permutation_closure(uniq)
    print(f"[perm-aware] {uniq.shape[0]} unordered rows, {lam.size} shells, "
          f"key {_KEY_DECIMALS} decimals, permutation closure "
          f"{100.0 * closure:.1f}%, from {TABLE_PATH.name}")
    # A handful of rows sit exactly on the triangle-inequality boundary,
    # where one angle is 180 degrees and the three sum to 360, and the
    # build's filter evaluates asymmetrically there in floating point. Those
    # are fully degenerate configurations that no fold queries, so the
    # threshold allows for them and only a real gap trips the warning.
    if closure < 0.99:
        print("[perm-aware] WARNING: the table is not closed under leg "
              "permutation. Strict mode rejects missing orderings; explicit "
              "ALLOW_INTERPOLATION permits a nearest-row blend. Rebuild with "
              "make_perm_closed_triples.py.")
    return _CACHE


def permutation_closure(uniq: NDArray[np.float64]) -> float:
    """Fraction of rows whose every leg reordering is also tabulated.

    A permutation-aware lookup asks for the vertex at reordered geometries.
    If those rows are absent the nearest-neighbour blend silently answers
    with a different configuration, so this is checked once at load rather
    than discovered in a result.
    """
    present = {tuple(row) for row in np.round(uniq, _KEY_DECIMALS)}
    if not present:
        return 1.0
    complete = 0
    for row in uniq:
        if all(tuple(np.round(_permute_triple(row, perm), _KEY_DECIMALS)) in present
               for perm in _PERMS):
            complete += 1
    return complete / uniq.shape[0]


def _permute_triple(triple: NDArray[np.float64], perm) -> NDArray[np.float64]:
    """Cosine triple of the same three directions in a new leg order.

    The triple is (cos01, cos12, cos20) for legs 0, 1, 2, so the entry for a
    pair is looked up by which two legs it joins, not by position.
    """
    a, b, c = perm
    pair = {(0, 1): 0, (1, 0): 0, (1, 2): 1, (2, 1): 1, (2, 0): 2, (0, 2): 2}
    idx = [pair[(a, b)], pair[(b, c)], pair[(c, a)]]
    return triple[..., idx]


def _locate(queries: NDArray[np.float64], lam: NDArray[np.float64]) -> dict:
    """Neighbour indices and weights for one leg ordering.

    Separated from the channel read because the tree query dominates the
    cost and depends only on the ordering, while up to six channels are
    read at each ordering. Doing it the other way round, which is the
    obvious way to write it, costs about four times as much.
    """
    cache = _ensure_cache()
    grid = cache["lambda"]
    k = min(4, cache["n_uniq"])
    dist, idx = cache["tree"].query(queries, k=k, workers=-1)
    if k == 1:
        dist, idx = dist[:, None], idx[:, None]
    i_hi = np.minimum(np.searchsorted(grid, lam, side="right"), grid.size - 1)
    i_lo = np.maximum(i_hi - 1, 0)
    span = grid[i_hi] - grid[i_lo]
    w_lam = np.where(span > 0, (lam - grid[i_lo]) / span, 0.0)
    exact = dist[:, 0] < _EXACT_TOL
    if not ALLOW_INTERPOLATION and not np.all(exact):
        worst = int(np.argmax(dist[:, 0]))
        raise ValueError(
            f"{int(np.sum(~exact))} of {exact.size} queries are not tabulated rows of "
            f"{TABLE_PATH.name} (nearest-row distance up to {dist[worst, 0]:.3e}, "
            f"tolerance {_EXACT_TOL:g}; first offending cosine triple "
            f"{np.round(queries[worst], 12).tolist()}). Set ALLOW_INTERPOLATION = True "
            "to blend the four nearest rows on purpose.")
    w = 1.0 / np.maximum(dist, 1e-12)
    return {"idx": idx, "i_lo": i_lo[:, None], "i_hi": i_hi[:, None],
            "w_lam": w_lam[:, None], "w": w, "w_sum": np.sum(w, axis=1),
            "exact": exact}


def _row_pchippers(cache: dict, channel: str, row: int) -> dict:
    """Lazy interval-to-spline map owned by this loaded table cache.

    A run never crosses a zero or a sign change. Two-point runs are valid
    PCHIP inputs and reduce to the signed-log linear interpolation rule.
    """
    key = (channel, int(row))
    splines = cache["pchippers"]
    if key not in splines:
        grid = cache["lambda"]
        values = cache["per_channel"][channel][row]
        intervals = {}
        start = 0
        while start < values.size:
            if values[start] == 0.0:
                start += 1
                continue
            stop = start + 1
            while (stop < values.size and values[stop] != 0.0
                   and np.signbit(values[stop]) == np.signbit(values[start])):
                stop += 1
            if stop - start >= 2:
                spline = PchipInterpolator(
                    grid[start:stop], np.log(np.abs(values[start:stop])),
                    extrapolate=False)
                for interval in range(start, stop - 1):
                    intervals[interval] = spline
            start = stop
        splines[key] = intervals
    return splines[key]


def _read(channel: str, at: dict) -> NDArray[np.float64]:
    """Evaluate each native row radially before an explicit angular blend."""
    cache = _ensure_cache()
    z = cache["per_channel"][channel]
    idx = at["idx"]
    lo_idx = np.broadcast_to(at["i_lo"], idx.shape)
    hi_idx = np.broadcast_to(at["i_hi"], idx.shape)
    lo = z[idx, lo_idx]
    hi = z[idx, hi_idx]
    weight = np.broadcast_to(at["w_lam"], idx.shape)
    z_at = (1.0 - weight) * lo + weight * hi
    interior = (weight != 0.0) & (weight != 1.0)
    same_sign = (lo != 0.0) & (hi != 0.0) & (np.signbit(lo) == np.signbit(hi))
    selected = interior & same_sign
    # Group by native row and interval; each spline uses the entire sign run.
    for row in np.unique(idx[selected]):
        runs = _row_pchippers(cache, channel, int(row))
        row_mask = selected & (idx == row)
        for interval in np.unique(lo_idx[row_mask]):
            mask = row_mask & (lo_idx == interval)
            grid = cache["lambda"]
            query = grid[lo_idx[mask]] + weight[mask] * (
                grid[hi_idx[mask]] - grid[lo_idx[mask]])
            z_at[mask] = np.sign(lo[mask]) * np.exp(runs[int(interval)](query))
    # Bypass all interpolation arithmetic at shells, including signed zeros.
    z_at = np.where(weight == 0.0, lo, np.where(weight == 1.0, hi, z_at))
    blended = np.sum(at["w"] * z_at, axis=1) / at["w_sum"]
    return np.where(at["exact"], z_at[:, 0], blended)


def _tensor_from_channels(channel_at, ns: int) -> NDArray[np.float64]:
    """Fill the (3, 3, 3) coupling tensor, each entry at its own geometry.

    ``channel_at(channel, perm)`` returns the channel evaluated on the legs
    reordered by ``perm``, as an array over samples.

    Component order is (Phi_00, Re Psi_0, Im Psi_0). Entries with an odd
    number of Im slots vanish by parity: under a reflection about the
    separation axis Re Psi_0 is even and Im Psi_0 is odd, so such a
    correlator equals minus itself. They are left at zero.
    """
    out = np.zeros((ns, 3, 3, 3), dtype=np.float64)
    identity = (0, 1, 2)

    for a in range(3):
        for b in range(3):
            for c in range(3):
                abc = (a, b, c)
                if sum(1 for x in abc if x == 2) % 2 == 1:
                    continue                       # parity zero
                spin = [i for i, x in enumerate(abc) if x != 0]

                if not spin:
                    out[:, a, b, c] = channel_at("zeta_TTT", identity)
                    continue

                if len(spin) == 1:
                    p = spin[0]
                    perm = tuple([i for i in range(3) if i != p] + [p])
                    out[:, a, b, c] = channel_at("zeta_TTP", perm)
                    continue

                if len(spin) == 2:
                    q = next(i for i in range(3) if i not in spin)
                    perm = (q, spin[0], spin[1])
                    tpp = channel_at("zeta_TPP", perm)
                    bmod = channel_at("zeta_Bmod", perm)
                    # Both spin legs carry the same component in any entry
                    # that survives parity, so their order does not matter.
                    out[:, a, b, c] = (0.5 * (tpp + bmod) if abc[spin[0]] == 1
                                       else 0.5 * (bmod - tpp))
                    continue

                # three spin legs, all Re: zeta_Dmod conjugates its slot-2 leg,
                # so the three cyclic orderings conjugate each leg once and
                # (PPP + D[0] + D[1] + D[2]) / 4 = Re Re Re (sympy check:
                # T-001 evidence checks_V6-refute-math/r7_small_symbolic)
                if all(x == 1 for x in abc):
                    out[:, a, b, c] = 0.25 * (channel_at("zeta_PPP", identity)
                                              + sum(channel_at("zeta_Dmod", p)
                                                    for p in _CYCLIC))
                    continue

                # one Re and two Im. (Dmod - PPP)/4 on an ordering (o0,o1,o2)
                # is the mean of the Re-on-o0 and Re-on-o1 placements, so the
                # three cyclic orderings give three pairwise means and the
                # wanted placement is mean(A) - mean(B) + mean(C).
                plus = next(i for i, x in enumerate(abc) if x == 1)
                others = [i for i in range(3) if i != plus]
                cycles = [(plus, others[0], others[1]),
                          (others[0], others[1], plus),
                          (others[1], plus, others[0])]
                means = [0.25 * (channel_at("zeta_Dmod", p)
                                 - channel_at("zeta_PPP", p)) for p in cycles]
                out[:, a, b, c] = means[0] - means[1] + means[2]
    return out


def coupling_fn_batch(n_arr, t_arr) -> NDArray[np.float64]:
    """Vectorised vertex. n_arr: (3, n, 3); t_arr: (3, n). Returns (n,3,3,3)."""
    n_arr = np.asarray(n_arr, dtype=float)
    t_arr = np.asarray(t_arr, dtype=float)
    if n_arr.ndim != 3 or n_arr.shape[0] != 3 or n_arr.shape[2] != 3:
        raise ValueError(f"n_arr must be (3, n_samples, 3); got {n_arr.shape}")
    if t_arr.shape != (3, n_arr.shape[1]):
        raise ValueError(f"t_arr must be (3, n_samples); got {t_arr.shape}")

    ns = n_arr.shape[1]
    nn = n_arr / np.maximum(np.linalg.norm(n_arr, axis=-1, keepdims=True), 1e-300)
    triple = np.stack([
        np.clip(np.einsum("ij,ij->i", nn[0], nn[1]), -1.0, 1.0),
        np.clip(np.einsum("ij,ij->i", nn[1], nn[2]), -1.0, 1.0),
        np.clip(np.einsum("ij,ij->i", nn[2], nn[0]), -1.0, 1.0),
    ], axis=-1)
    lam = np.mean(t_arr, axis=0)

    cache = _ensure_cache()
    inside = np.flatnonzero((lam >= cache["lam_min"]) & (lam <= cache["lam_max"]))
    out = np.zeros((ns, 3, 3, 3), dtype=np.float64)
    if inside.size == 0:
        return out

    # Evaluate every (channel, ordering) pair once and reuse it across the
    # tensor entries that need it; there are at most six orderings.
    located: dict[tuple[int, int, int], dict] = {}
    memo: dict[tuple[str, tuple[int, int, int]], NDArray[np.float64]] = {}

    def channel_at(channel: str, perm) -> NDArray[np.float64]:
        perm = tuple(perm)
        if perm not in located:
            located[perm] = _locate(_permute_triple(triple[inside], perm),
                                    lam[inside])
        key = (channel, perm)
        if key not in memo:
            memo[key] = _read(channel, located[perm])
        return memo[key]

    out[inside] = _tensor_from_channels(channel_at, inside.size)
    return out


def coupling_fn(n_list, t_list) -> NDArray[np.float64]:
    """Scalar entry point, expressed through the batch one."""
    n_arr = np.asarray(n_list, dtype=float).reshape(3, 1, 3)
    t_arr = np.asarray(t_list, dtype=float).reshape(3, 1)
    return coupling_fn_batch(n_arr, t_arr)[0]


coupling_fn_batch.vectorized = True  # type: ignore[attr-defined]
