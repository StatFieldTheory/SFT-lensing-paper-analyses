import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
from canoes.sachs.sft_input.corr_op.table import make_C_fn_lambda_project
cubic = make_C_fn_lambda_project(co._TABLE, interp_method="cubic", output_h_power=co.OUTPUT_H_POWER)
n2 = ds._n2_of_cos(float(np.cos(np.deg2rad(1.0/60.0))))
for t in (1373.5, 1300.0, 2160.5, 800.0):
    a = np.asarray(co.C_fn(ds._N1, t, n2, t), float)
    b = np.asarray(cubic(ds._N1, t, n2, t), float)
    print("t=%8.1f  linear %.10e  cubic %.10e  reldiff %.3e"
          % (t, a[0,0], b[0,0], abs(b[0,0]/a[0,0]-1)))
