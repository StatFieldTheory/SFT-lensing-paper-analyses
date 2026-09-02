"""Does the ridge also bias a TWO-DIMENSIONAL integral over the cell?

The Sigma2 consumer queries the diagonal, where the bilinear error is largest
and one-signed (-10 to -24%).  The sft-wick folds instead integrate C over
(lam1, lam2), where the +4% far corners partly offset the -16% ridge.  This
measures the residual by trapezoid-integrating the production and dense values
over the SAME cell on the dense 5x5 grid (a 4x4-panel trapezoid, so the ridge is
resolved at half the production spacing).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/zzhang/projects/angular_statistics/canoes/src")
from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    load_table_2d_from_npz, make_C_fn_lambda_project)

HERE = Path(__file__).resolve().parent
PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
            "sachs_sft/callables/C_propagator/corr_op/"
            "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz")
N1 = np.array([0.0, 0.0, 1.0])
Cp = make_C_fn_lambda_project(load_table_2d_from_npz(str(PROD)))

print(f"{'cell':>16} {'diag error at mid':>18} {'2-D cell-integral error':>24}")
for tag in ("cellA5", "cellB5", "cellC5", "cellD5"):
    p = HERE / "dense" / f"corr_op_dense_{tag}.npz"
    if not p.exists():
        continue
    tab = load_table_2d_from_npz(str(p))
    lam = np.asarray(tab.cl_table.lambda_grid, float)
    Cd = make_C_fn_lambda_project(tab)
    Md = np.array([[Cd(N1, float(a), N1, float(b))[0, 0] for b in lam] for a in lam])
    Mp = np.array([[Cp(N1, float(a), N1, float(b))[0, 0] for b in lam] for a in lam])
    Id = np.trapezoid(np.trapezoid(Md, lam, axis=1), lam)
    Ip = np.trapezoid(np.trapezoid(Mp, lam, axis=1), lam)
    m = len(lam) // 2
    print(f"{lam[0]:7.1f}-{lam[-1]:7.1f} {100*(Mp[m,m]/Md[m,m]-1):17.2f}% "
          f"{100*(Ip/Id-1):23.2f}%")
