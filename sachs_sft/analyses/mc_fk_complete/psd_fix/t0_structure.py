"""What the production table's lambda interpolation actually is, from the data.

RegularGridInterpolator((lam, lam), M, method='linear') is BILINEAR.  Restricted
to the diagonal t1 = t2 = lam_i + u*h it is the quadratic Bezier

    B(u) = (1-u)^2 F_ii + 2 u(1-u) F_{i,i+1} + u^2 F_{i+1,i+1},

so the mid-cell diagonal value is 0.25 F_ii + 0.50 F_{i,i+1} + 0.25 F_{i+1,i+1}
-- it leans on the OFF-DIAGONAL table entry.  This script prints the three
control values per cell, the implied second derivative, and the size of the
"Bezier pull" 2*(F_{i,i+1} - (F_ii+F_{i+1,i+1})/2), which is exactly what a
midpoint truth value has to be compared against.
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


def main() -> None:
    tab = load_table_2d_from_npz(str(PROD))
    lam = np.asarray(tab.cl_table.lambda_grid, float)
    C = make_C_fn_lambda_project(tab)
    print(f"table {PROD.name}\n  lambda nodes ({lam.size}): {lam.tolist()}")
    print(f"  interp: RegularGridInterpolator((lam,lam), method='linear') "
          f"-> BILINEAR in (lam1,lam2), LINEAR in lambda (not log)")
    for cosg, tag in ((1.0, "cos=1"),):
        n2 = n2_of_cos(cosg)
        print(f"\n== component (0,0) = <Phi00 Phi00>, {tag} ==")
        print(f"{'cell':>16} {'h':>7} {'F_ii':>12} {'F_i,i+1':>12} {'F_jj':>12} "
              f"{'B(0.5)':>12} {'pull/B':>9} {'d2 (Bezier)':>12}")
        for i in range(lam.size - 1):
            a, b = float(lam[i]), float(lam[i + 1])
            h = b - a
            Fii = C(N1, a, n2, a)[0, 0]
            Fij = C(N1, a, n2, b)[0, 0]
            Fjj = C(N1, b, n2, b)[0, 0]
            B = 0.25 * Fii + 0.5 * Fij + 0.25 * Fjj
            pull = 2.0 * (Fij - 0.5 * (Fii + Fjj))
            d2 = 2.0 * (Fii - 2.0 * Fij + Fjj) / h ** 2
            print(f"{a:7.1f}-{b:7.1f} {h:7.1f} {Fii:12.5e} {Fij:12.5e} "
                  f"{Fjj:12.5e} {B:12.5e} {0.25*pull/B:9.2e} {d2:12.5e}")


if __name__ == "__main__":
    main()
