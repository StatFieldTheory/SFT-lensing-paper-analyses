"""Triple set closed under the leg permutations the fold actually queries.

The vertex is a function of three (field, direction) pairs. Its channels are
labelled by which legs carry the spin-2 field: zeta_TTP has spin weights
(0, 0, 2), so the spin leg is the THIRD one, and zeta_TPP has (0, 2, 2), so
the scalar leg is the FIRST. A lookup that needs the spin leg in a
different slot must therefore evaluate the vertex at a PERMUTED geometry,
not at the same geometry with a relabelled channel.

The production mesh is already closed under permutation, being a full
meshgrid with a symmetric triangle filter. The densified collapsed family
is not: it contains only (1, c, c). This script adds the other two
placements of the offset leg, (c, 1, c) and (c, c, 1), so a
permutation-aware callable can find every row it needs.

Usage::

    python make_perm_closed_triples.py --refine 4 \
        --out ../../products/triples_permclosed_r4.npz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from make_dense_triples import (GAMMA_ARCMIN, collapsed_cosines,  # noqa: E402
                                production_triples)


def collapsed_orderings(cos_values: np.ndarray) -> np.ndarray:
    """(1,c,c), (c,1,c) and (c,c,1) for every collapsed cosine.

    In the callable's convention the triple is (cos12, cos23, cos31), so
    cos23 = 1 means legs 2 and 3 coincide and leg 1 carries the offset.
    The three rows below are the three choices of which leg is offset.
    """
    one = np.ones_like(cos_values)
    return np.vstack([
        np.column_stack([one, cos_values, cos_values]),   # leg 3 offset
        np.column_stack([cos_values, one, cos_values]),   # leg 1 offset
        np.column_stack([cos_values, cos_values, one]),   # leg 2 offset
    ])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refine", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    prod = production_triples()
    cos_values = collapsed_cosines(args.refine)
    collapsed = collapsed_orderings(cos_values)
    stacked = np.vstack([prod, collapsed])
    _, keep = np.unique(np.round(stacked, 12), axis=0, return_index=True)
    triples = stacked[np.sort(keep)]

    np.savez(args.out, cosine_triples=triples,
             n_production=np.int64(prod.shape[0]),
             refine=np.int64(args.refine),
             gamma_arcmin_production=GAMMA_ARCMIN,
             note=np.asarray(
                 "production mesh plus the collapsed family in all three "
                 "placements of the offset leg, so a permutation-aware "
                 "lookup finds every row it needs", dtype=object))
    print(f"{triples.shape[0]} triples -> {args.out}")
    print(f"  production mesh          : {prod.shape[0]}")
    print(f"  collapsed, three orderings: {collapsed.shape[0]}")
    print(f"  added by the union       : {triples.shape[0] - prod.shape[0]}")


if __name__ == "__main__":
    main()
