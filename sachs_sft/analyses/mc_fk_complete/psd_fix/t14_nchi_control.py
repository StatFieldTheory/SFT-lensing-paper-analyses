"""Is the dense 'truth' itself converged in the inner chi quadrature?

Same 5-node lambda grid, n_chi_per_pair 256 vs 512 (the builder snaps to
2^j+1, and the snapped grids NEST, so the difference is a genuine Cauchy step).
If the two agree far better than the 10-16% interpolation error being reported,
the dense rebuild is a valid yardstick.
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

t256 = load_table_2d_from_npz(str(HERE / "dense" / "corr_op_dense_cellA5.npz"))
t512 = load_table_2d_from_npz(str(HERE / "dense" / "corr_op_dense_cellA5_nchi512.npz"))
C256 = make_C_fn_lambda_project(t256)
C512 = make_C_fn_lambda_project(t512)
Cp = make_C_fn_lambda_project(load_table_2d_from_npz(str(PROD)))
lam = np.asarray(t256.cl_table.lambda_grid, float)
print(f"{'lambda':>9} {'n_chi=256':>14} {'n_chi=512':>14} {'512/256-1 [%]':>14} "
      f"{'prod/512-1 [%]':>15}")
for t in lam:
    a = C256(N1, float(t), N1, float(t))[0, 0]
    b = C512(N1, float(t), N1, float(t))[0, 0]
    p = Cp(N1, float(t), N1, float(t))[0, 0]
    print(f"{t:9.2f} {a:14.6e} {b:14.6e} {100*(b/a-1):14.4f} {100*(p/b-1):15.3f}")
