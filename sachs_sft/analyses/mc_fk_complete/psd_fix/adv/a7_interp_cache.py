"""Two things:
 (1) is corr_op's interp_method silently ignored once a (cos,psi1,psi2) triple
     has been touched?  _get_or_build_channels keys its cache on that triple
     ONLY -- interp_method is not part of the key.
 (2) does R5Exact's structural guard fire on a deliberately non-bilinear C?
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
import candidates as K
from canoes.sachs.sft_input.corr_op.table import make_C_fn_lambda_project

ORDER = sys.argv[1] if len(sys.argv) > 1 else "cubic_first"
n2 = ds._n2_of_cos(float(np.cos(np.deg2rad(1.0/60.0))))
TS = (800.0, 1373.5, 2160.5)

cub = make_C_fn_lambda_project(co._TABLE, interp_method="cubic",
                               output_h_power=co.OUTPUT_H_POWER)
if ORDER == "cubic_first":
    c = [np.asarray(cub(ds._N1, t, n2, t), float)[0,0] for t in TS]
    l = [np.asarray(co.C_fn(ds._N1, t, n2, t), float)[0,0] for t in TS]
else:
    l = [np.asarray(co.C_fn(ds._N1, t, n2, t), float)[0,0] for t in TS]
    c = [np.asarray(cub(ds._N1, t, n2, t), float)[0,0] for t in TS]
print(f"=== order = {ORDER} ===")
for t, a, b in zip(TS, l, c):
    print("  t=%8.1f  linear %.10e  cubic %.10e  reldiff %.3e"
          % (t, a, b, abs(b/a - 1)))
