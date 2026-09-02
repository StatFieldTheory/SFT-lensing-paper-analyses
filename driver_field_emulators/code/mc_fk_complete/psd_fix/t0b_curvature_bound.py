"""Second-order bound on the mid-cell diagonal interpolation error, table-only.

For a SMOOTH f(x, y) = C(lam1, lam2) the bilinear interpolant on the cell
[a,b]^2 evaluated at the cell-diagonal midpoint m = (a+b)/2 satisfies

    B(m, m) - f(m, m) = (h^2 / 8) (f_xx + f_yy) = (h^2 / 4) f_xx        (1)

(the mixed term cancels exactly on the diagonal), and the two table-visible
combinations are

    F_ii - 2 F_ij + F_jj = h^2 f_xy    (the interpolant's diagonal curvature
                                        is 2 f_xy, NOT the true 2 f_xx + 2 f_xy)

so the diagonal interpolant misses exactly the f_xx part of the curvature.
f_xx is measurable from the production table by differencing along ROW lam2 =
const over three consecutive nodes (unequal spacing -> divided differences).

Caveat printed with the result: C(lam1, lam2) may be non-smooth ACROSS the
diagonal (the theta(lam - chi) window cuts at min(lam1, lam2)), in which case
(1) underestimates the error.  That is why a dense rebuild is the real test;
this is the cheap bound.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    load_table_2d_from_npz, make_C_fn_lambda_project,
)

PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
            "sachs_sft/callables/C_propagator/corr_op/"
            "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz")
N1 = np.array([0.0, 0.0, 1.0])


def n2_of_cos(c: float) -> np.ndarray:
    c = float(np.clip(c, -1.0, 1.0))
    return np.array([np.sqrt(max(0.0, 1.0 - c * c)), 0.0, c])


def second_div_diff(x, y):
    """f'' from three unequally spaced points (Newton divided difference * 2)."""
    x0, x1, x2 = x
    y0, y1, y2 = y
    return 2.0 * (((y2 - y1) / (x2 - x1)) - ((y1 - y0) / (x1 - x0))) / (x2 - x0)


def main() -> None:
    tab = load_table_2d_from_npz(str(PROD))
    lam = np.asarray(tab.cl_table.lambda_grid, float)
    C = make_C_fn_lambda_project(tab)
    n2 = n2_of_cos(1.0)
    comp = (0, 0)

    print("cell-midpoint bilinear error estimate  eps = (h^2/4) f_xx , "
          "f_xx from row differencing at fixed lam2 (one node BELOW the diagonal)")
    print(f"{'cell':>16} {'h':>6} {'f(m,m)~':>11} {'f_xx':>11} "
          f"{'eps=(h^2/4)fxx':>15} {'eps/f  [%]':>11}")
    for i in range(1, lam.size - 2):
        a, b = float(lam[i]), float(lam[i + 1])
        h = b - a
        m = 0.5 * (a + b)
        # f_xx along lam1 at fixed lam2 = a (strictly inside the lam1 >= lam2
        # half-plane for the two upper nodes; the row is smooth there).
        xs = [float(lam[i]), float(lam[i + 1]), float(lam[i + 2])]
        ys = [C(N1, x, n2, a)[comp] for x in xs]
        fxx = second_div_diff(xs, ys)
        fmm = 0.5 * (C(N1, a, n2, a)[comp] + C(N1, b, n2, b)[comp])
        eps = 0.25 * h * h * fxx
        print(f"{a:7.1f}-{b:7.1f} {h:6.1f} {fmm:11.4e} {fxx:11.4e} "
              f"{eps:15.4e} {100*eps/fmm:11.3f}")


if __name__ == "__main__":
    main()
