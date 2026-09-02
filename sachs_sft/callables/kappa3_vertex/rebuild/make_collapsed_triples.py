"""Generate the densified collapsed-family cosine grid for the FK vertex.

The production FK two-point fold queries the kappa3 table ONLY on the
collapsed family (1, cos gamma, cos gamma), at the 40 gamma values of the
production grid. On the deployed covgrid16 mesh (15 uniform cosines,
spacing 1/7) the callable's inverse-distance lookup puts at least 96% of
its weight on the single (1, 1, 1) corner for 27 of those 40 points, so the
table cannot resolve the vertex's gamma dependence below about 232 arcmin.

This script writes a cosine grid that contains the production family
directly, plus a log-refined set in (1 - cos gamma) so that interpolation
between queried points is along the family rather than across unrelated
open triangles.

Usage::

    python make_collapsed_triples.py --out collapsed_triples_dense.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

#: The 40 production separations, arcmin, log-spaced over [0.5', 5000'].
GAMMA_ARCMIN = np.geomspace(0.5, 5000.0, 40)
ARCMIN = np.pi / (180.0 * 60.0)


def collapsed_grid(refine: int = 2) -> np.ndarray:
    """Cosine values of the collapsed family, production points plus refinement.

    Refinement is uniform in log(1 - cos gamma), the variable in which the
    production points are (nearly) uniform, so no interval is left coarse.
    """
    cos_production = np.cos(GAMMA_ARCMIN * ARCMIN)
    one_minus = 1.0 - cos_production
    log_grid = np.log(one_minus)
    if refine > 1:
        fine = np.concatenate(
            [
                np.linspace(log_grid[i], log_grid[i + 1], refine + 1)[:-1]
                for i in range(log_grid.size - 1)
            ]
            + [log_grid[-1:]]
        )
    else:
        fine = log_grid
    cos_values = 1.0 - np.exp(fine)
    # The exact corner is what the smallest separations approach; keep it.
    cos_values = np.unique(np.clip(np.append(cos_values, 1.0), -1.0, 1.0))
    return np.sort(cos_values)[::-1]


def collapsed_triples(cos_values: np.ndarray) -> np.ndarray:
    """Rows (1, c, c), already in the descending-sorted canonical order."""
    return np.column_stack(
        [np.ones_like(cos_values), cos_values, cos_values]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--refine", type=int, default=2,
                        help="Sub-intervals per production gap (1 = none).")
    args = parser.parse_args()

    cos_values = collapsed_grid(args.refine)
    triples = collapsed_triples(cos_values)
    np.savez(
        args.out,
        cosine_triples=triples,
        cos_values=cos_values,
        gamma_arcmin_production=GAMMA_ARCMIN,
        note=np.asarray(
            "collapsed family (1, c, c) only: valid for the FK two-point fold, "
            "NOT for the order-1 three-point vertex, which needs open triangles",
            dtype=object,
        ),
    )
    print(f"{triples.shape[0]} collapsed triples -> {args.out}")
    print(f"  cos range [{cos_values.min():.6f}, {cos_values.max():.6f}]")
    print(f"  1 - cos range [{1 - cos_values.max():.3e}, {1 - cos_values.min():.3e}]")


if __name__ == "__main__":
    main()
