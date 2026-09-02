"""Build a kappa3 cosine-triple set: production mesh plus a densified
collapsed family.

The deployed vertex table is indexed by triangle-filtered cosine triples
drawn from a uniform 15-point mesh in cos, spacing 1/7. The FK two-point
fold queries only the collapsed family (1, cos g, cos g), and the plotted
separations 0.5' to 5000' map to 1 - cos g between 1e-8 and 8.8e-1, so
almost the whole plotted range lies inside the single mesh cell touching
the (1, 1, 1) corner.

`make_collapsed_triples.py` writes the collapsed family alone, which is
enough to fold the FK two-point function and nothing else. This script
writes the UNION of

  * the production triangle-filtered mesh (n_cosine = 15 gives 1671 rows),
    so the same table still serves the order-1 three-point vertex, the
    open-triangle checks and the multi-redshift figure, and
  * the collapsed family refined in log(1 - cos g),

so that one table serves every consumer and the only change relative to
production is added resolution.

Usage::

    python make_dense_triples.py --refine 2 --out products/triples_full_r2.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

#: The 40 production separations, arcmin, log-spaced over [0.5', 5000'].
GAMMA_ARCMIN = np.geomspace(0.5, 5000.0, 40)
ARCMIN = np.pi / (180.0 * 60.0)

#: Cosine mesh size of the deployed table (covgrid16 -> 15 uniform values).
N_COSINE_PRODUCTION = 15


def production_triples(n_cosine: int = N_COSINE_PRODUCTION) -> np.ndarray:
    """Triangle-filtered (cos g12, cos g23, cos g31) mesh.

    Verbatim reimplementation of `_triangle_filtered_triples` in
    `build_equal_time_limber_table.py`, so the production rows in the union
    are bit-identical to the ones the deployed table was built on.
    """
    c = np.linspace(-1.0, 1.0, n_cosine)
    c12, c23, c31 = np.meshgrid(c, c, c, indexing="ij")
    tr = np.stack([c12.ravel(), c23.ravel(), c31.ravel()], axis=1)
    g = np.arccos(np.clip(tr, -1, 1))
    g12, g23, g31 = g[:, 0], g[:, 1], g[:, 2]
    valid = ((g12 + g23 >= g31) & (g23 + g31 >= g12) & (g31 + g12 >= g23)
             & (g12 + g23 + g31 <= 2 * np.pi))
    return tr[valid]


def collapsed_cosines(refine: int) -> np.ndarray:
    """Cosines of the collapsed family: production points plus refinement.

    Refinement is uniform in log(1 - cos g) because the production
    separations are log-spaced in g, which makes them very nearly
    log-spaced in 1 - cos g as well. `refine = 1` reproduces the 40
    production points and adds nothing.
    """
    cos_production = np.cos(GAMMA_ARCMIN * ARCMIN)
    log_grid = np.log(1.0 - cos_production)
    if refine > 1:
        fine = np.concatenate(
            [np.linspace(log_grid[i], log_grid[i + 1], refine + 1)[:-1]
             for i in range(log_grid.size - 1)] + [log_grid[-1:]]
        )
    else:
        fine = log_grid
    cos_values = 1.0 - np.exp(fine)
    cos_values = np.append(cos_values, 1.0)          # the exact corner
    return np.sort(np.unique(np.clip(cos_values, -1.0, 1.0)))[::-1]


def union_triples(refine: int, n_cosine: int = N_COSINE_PRODUCTION):
    """Production mesh plus collapsed family, deduplicated."""
    prod = production_triples(n_cosine)
    cos_values = collapsed_cosines(refine)
    collapsed = np.column_stack(
        [np.ones_like(cos_values), cos_values, cos_values])
    stacked = np.vstack([prod, collapsed])
    # Deduplicate on a tolerance far below the refinement spacing: rows are
    # compared at 1e-12 in cos, while the finest refinement step near the
    # corner is about 1e-9, so no distinct refined row is merged away.
    _, keep = np.unique(np.round(stacked, 12), axis=0, return_index=True)
    return stacked[np.sort(keep)], prod, collapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--refine", type=int, default=2,
                        help="Sub-intervals per production gap (1 = none).")
    parser.add_argument("--n-cosine", type=int, default=N_COSINE_PRODUCTION)
    args = parser.parse_args()

    triples, prod, collapsed = union_triples(args.refine, args.n_cosine)
    one_minus = 1.0 - collapsed[:, 1]
    positive = one_minus[one_minus > 0]
    np.savez(
        args.out,
        cosine_triples=triples,
        n_production=np.int64(prod.shape[0]),
        n_collapsed=np.int64(collapsed.shape[0]),
        refine=np.int64(args.refine),
        gamma_arcmin_production=GAMMA_ARCMIN,
        note=np.asarray(
            "union of the production triangle-filtered mesh and the "
            "collapsed (1, c, c) family refined in log(1 - cos g); serves "
            "both the FK two-point fold and the open-triangle consumers",
            dtype=object),
    )
    print(f"refine={args.refine}: {triples.shape[0]} triples -> {args.out}")
    print(f"  production mesh      : {prod.shape[0]}")
    print(f"  collapsed family     : {collapsed.shape[0]}")
    print(f"  added by the union   : {triples.shape[0] - prod.shape[0]}")
    print(f"  finest 1-cos step    : {np.diff(np.sort(positive)).min():.3e}")
    print(f"  1-cos range (collapsed): [{positive.min():.3e}, "
          f"{one_minus.max():.3e}]")


if __name__ == "__main__":
    main()
