"""Table-only proof that the production diagonal interpolant is unphysical.

On the cell [a,b] the bilinear interpolant restricted to the diagonal is the
quadratic Bezier with control points (F_ii, F_ij, F_jj), so its slope at the
cell start is  2 (F_ij - F_ii) / h.  Because C(lam1,lam2) has a RIDGE on the
diagonal (F_ij < F_ii even though F_jj >> F_ii), that slope is NEGATIVE: the
interpolated C(lam,lam) decreases just past every table node, then over-steepens
to catch up.  The cell-mean slope is exactly (F_jj - F_ii)/h -- correct -- which
is why an integral of dC/dlam over a whole cell telescopes and looks healthy
while the pointwise density is a sawtooth.

No rebuild is needed for this; it is arithmetic on the deployed table.
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


def main() -> None:
    tab = load_table_2d_from_npz(str(PROD))
    lam = np.asarray(tab.cl_table.lambda_grid, float)
    C = make_C_fn_lambda_project(tab)
    n2 = N1.copy()
    print(f"{'cell':>16} {'h':>6} {'slope(u=0)':>12} {'slope(u=1)':>12} "
          f"{'cell mean':>12} {'u* (slope=0)':>13} {'B(0.5)/F_ii':>12}")
    for i in range(lam.size - 1):
        a, b = float(lam[i]), float(lam[i + 1])
        h = b - a
        Fii = C(N1, a, n2, a)[0, 0]
        Fij = C(N1, a, n2, b)[0, 0]
        Fjj = C(N1, b, n2, b)[0, 0]
        s0 = 2.0 * (Fij - Fii) / h
        s1 = 2.0 * (Fjj - Fij) / h
        sm = (Fjj - Fii) / h
        ustar = s0 / (s0 - s1) if s1 != s0 else np.nan
        B = 0.25 * Fii + 0.5 * Fij + 0.25 * Fjj
        print(f"{a:7.1f}-{b:7.1f} {h:6.1f} {s0:12.4e} {s1:12.4e} {sm:12.4e} "
              f"{ustar:13.3f} {B/Fii:12.4f}")
    print("\nslope(u=0) < 0 in every cell -> the interpolated C(lam,lam) DECREASES "
          "just past each table node")
    print("u* is the fraction of the cell over which it decreases before turning")


if __name__ == "__main__":
    main()
