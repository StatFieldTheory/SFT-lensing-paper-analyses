"""V4: verify the table-resolution claim from the dense rebuilds already on disk.
prod bilinear C00(t,t) vs the dense-rebuild truth at mid-cell lambdas."""
import sys
from pathlib import Path
import numpy as np
_CO = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/C_propagator/corr_op")
sys.path.insert(0, str(_CO))
import corr_op_C_callable as prod
from canoes.sachs.sft_input.corr_op.table import load_table_2d_from_npz, make_C_fn_lambda_project

N1 = np.array([0.,0.,1.])
def n2_of_cos(c):
    c = float(np.clip(c,-1,1)); s = float(np.sqrt(max(0.,1-c*c)))
    return np.array([s,0.,c])

DENSE = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix/dense")
for name in ("cellC5","cellA5","cellB5","cellD5"):
    t = load_table_2d_from_npz(str(DENSE/f"corr_op_dense_{name}.npz"))
    fn = make_C_fn_lambda_project(t, output_h_power=0.0)
    lam = np.load(DENSE/f"corr_op_dense_{name}.npz", allow_pickle=True)["lambda_grid"].astype(float)
    print(f"\n=== {name}  nodes {lam}")
    for c, tag in ((1.0,"cos=1"),):
        n2 = n2_of_cos(c)
        for x in lam:
            d = float(fn(N1, float(x), n2, float(x))[0,0])
            p = float(prod.C_fn(N1, float(x), n2, float(x))[0,0])
            mark = "NODE(shared)" if (x==lam[0] or x==lam[-1]) else "mid-cell"
            print(f"  lam={x:8.2f} {mark:>12}  dense {d:.6e}  prod {p:.6e}  prod/dense-1 {p/d-1:+.4%}")
