"""Validate the candidate (lam1, lam2) readings against the dense rebuilds.

The four `psd_fix/dense/corr_op_dense_cell?5.npz` tables are 5-node rebuilds of
ONE production cell each (same build arguments, same chi_min = 50), so their
own nodes are exact C values inside a cell where the production table has to
interpolate.  For each scheme we report

  * the pointwise error at the three interior dense nodes, on the equal-time
    diagonal (where the ridge lives) and at the off-diagonal interior points;
  * the error of the 2-D cell integral (tensor trapezoid on the dense 5x5
    grid, which is exact for any of the piecewise-planar schemes).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ridge_C import RidgeC  # noqa: E402
sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    load_table_2d_from_npz, make_C_fn_lambda_project)

HERE = Path(__file__).resolve().parent
DENSE = HERE.parent / "dense"
PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
            "sachs_sft/callables/C_propagator/corr_op/"
            "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz")
N1 = np.array([0.0, 0.0, 1.0])
SCHEMES = ("bilinear", "cubic", "tri", "tri_log", "ridge")
# a small-gamma and a large-gamma direction (the angular channels differ)
DIRS = {
    "gamma=0.5'":   np.array([0.0001454441, 0.0, 0.9999999894]),
    "gamma=17.3'":  np.array([0.0050252287, 0.0, 0.9999873735]),
    "gamma=293.9'": np.array([0.0853881735, 0.0, 0.9963477605]),
}


def main() -> None:
    servers = {s: RidgeC(str(PROD), s) for s in SCHEMES}
    for tag in ("cellA5", "cellB5", "cellC5", "cellD5"):
        p = DENSE / f"corr_op_dense_{tag}.npz"
        if not p.exists():
            continue
        dt = load_table_2d_from_npz(str(p))
        lam = np.asarray(dt.cl_table.lambda_grid, float)
        Cd = make_C_fn_lambda_project(dt)
        for dtag, n2 in DIRS.items():
            truth = np.array([[Cd(N1, float(a), n2, float(b))[0, 0] for b in lam]
                              for a in lam])
            Itr = np.trapezoid(np.trapezoid(truth, lam, axis=1), lam)
            interior = [1, 2, 3]
            print(f"\n=== {tag}  lam in [{lam[0]:.1f}, {lam[-1]:.1f}]  {dtag} ===")
            print(f"{'scheme':>10} {'diag mid':>10} {'diag q1':>9} {'diag q3':>9} "
                  f"{'off-diag worst':>15} {'cell integral':>14}")
            for s, srv in servers.items():
                M = np.array([[srv.C_fn(N1, float(a), n2, float(b))[0, 0]
                               for b in lam] for a in lam])
                rel = M / truth - 1.0
                dg = np.diag(rel)
                off = rel.copy()
                np.fill_diagonal(off, 0.0)
                offw = off[np.unravel_index(np.argmax(np.abs(off)), off.shape)]
                Iap = np.trapezoid(np.trapezoid(M, lam, axis=1), lam)
                print(f"{s:>10} {dg[2]*100:9.2f}% {dg[1]*100:8.2f}% {dg[3]*100:8.2f}% "
                      f"{offw*100:14.2f}% {(Iap/Itr-1)*100:13.2f}%")


if __name__ == "__main__":
    main()
