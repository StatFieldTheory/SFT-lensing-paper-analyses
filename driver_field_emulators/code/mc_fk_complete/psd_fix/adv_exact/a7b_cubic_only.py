import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
from canoes.sachs.sft_input.corr_op.table import make_C_fn_lambda_project
cub = make_C_fn_lambda_project(co._TABLE, interp_method="cubic", output_h_power=co.OUTPUT_H_POWER)
n2 = ds._n2_of_cos(float(np.cos(np.deg2rad(1.0/60.0))))
M = np.asarray(cub(ds._N1, 1373.5, n2, 1373.5), float)
print("cubic-only, fresh process, full (3,3) at t=1373.5:")
print(M)
print("all zero:", bool(np.all(M == 0.0)))
# and at a table node
print("at node 1300:", np.asarray(cub(ds._N1, 1300.0, n2, 1300.0), float)[0,0])
