"""Is the plotted FK curve the interpolation weight?

This probe needs only the DEPLOYED table and the DEPLOYED folded result.
It builds no second table, so it cannot be answered with "your rebuild
differs from theirs". It reconstructs the callable's own lookup, reads off
the weight the four-neighbour inverse-distance rule puts on the (1, 1, 1)
node at each plotted separation, and compares that weight with the shape of
the published FK curve.

If the two agree, the plotted angular dependence is the interpolator's
weight function and carries no information about the vertex.

Usage::

    python probe_weight_identity.py --out ../../products/weight_identity.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from observables import load_observables

_REPO = Path(__file__).resolve().parents[3]
DEPLOYED_TABLE = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
                  / "kappa3_vertex" / "equal_time_limber"
                  / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz")
DEPLOYED_FK = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "sftwick_outputs"
               / "2PCF" / "C_corr_op_K_limber_FK"
               / "xi_C_corr_op_K_limber_FK.npz")
ARCMIN = np.pi / (180.0 * 60.0)


def build_lookup(table: Path):
    """The callable's canonicalisation and tree, reconstructed verbatim.

    Mirrors `equal_time_limber_kappa3_callable.py:153-204`: sort each triple
    descending, round to 8 decimals, deduplicate, and build a cKDTree over
    the unique rows.
    """
    d = np.load(table, allow_pickle=False)
    triples = np.asarray(d["cosine_triples"], float)
    sorted_desc = np.sort(triples, axis=1)[:, ::-1]
    uniq = np.unique(np.round(sorted_desc, 8), axis=0)
    return uniq, cKDTree(uniq), np.asarray(d["zeta_TTT"], float)


def corner_weight(uniq, tree, gamma_arcmin: np.ndarray):
    """Normalised inverse-distance weight on the (1, 1, 1) node."""
    corner = int(np.argmin(np.linalg.norm(uniq - 1.0, axis=1)))
    if np.linalg.norm(uniq[corner] - 1.0) > 1e-12:
        raise SystemExit("the deployed table has no (1, 1, 1) row")
    out, neighbours = [], []
    for g in gamma_arcmin:
        c = np.cos(g * ARCMIN)
        query = np.array([1.0, c, c])
        query = np.sort(query)[::-1]
        distance, index = tree.query(query, k=min(4, uniq.shape[0]))
        distance = np.atleast_1d(distance)
        index = np.atleast_1d(index)
        weight = 1.0 / np.maximum(distance, 1e-12)
        share = weight / weight.sum()
        hit = np.flatnonzero(index == corner)
        out.append(float(share[hit[0]]) if hit.size else 0.0)
        neighbours.append([uniq[i].tolist() for i in index])
    return np.asarray(out), neighbours


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, default=DEPLOYED_TABLE)
    ap.add_argument("--fk", type=Path, default=DEPLOYED_FK)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    gamma, obs = load_observables(args.fk, order=2)
    curve = obs["xi_kappa"]
    shape = curve / curve[0]

    uniq, tree, _ = build_lookup(args.table)
    weight, neighbours = corner_weight(uniq, tree, gamma)

    print(f"deployed table: {uniq.shape[0]} unique sorted triples")
    header = ("gamma[arcmin]".rjust(13) + "FK/FK(0.5')".rjust(14)
              + "w(1,1,1)".rjust(12) + "difference".rjust(13) + "rel".rjust(11))
    print(header)
    for g, s, w in zip(gamma, shape, weight):
        difference = s - w
        rel = difference / w if w > 0 else np.nan
        print(f"{g:13.2f} {s:14.6f} {w:12.6f} {difference:+13.2e} {rel:+11.1e}")

    inside = weight > 1e-6
    residual = np.abs(shape[inside] - weight[inside])
    print(f"\nmax |FK shape - corner weight| where the corner is used: "
          f"{residual.max():.2e}")
    print(f"max relative           : {np.max(residual / weight[inside]):.2e}")
    zero_at = gamma[~inside]
    if zero_at.size:
        print(f"the corner leaves the four-neighbour set at gamma = "
              f"{zero_at.min():.1f} arcmin = {zero_at.min() / 60:.1f} degrees")
        print(f"  neighbours there: {neighbours[int(np.argmax(~inside))]}")

    if args.out:
        np.savez(args.out, gamma_arcmin=gamma, fk_shape=shape,
                 corner_weight=weight, fk_kappa=curve)
        print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
