"""V7: does halving the table cell halve the mid-cell error (O(h)) or quarter it (O(h^2))?
Uses the 3-node (h=73.5) and 5-node (h=36.75, truth) dense rebuilds already on disk."""
import sys
from pathlib import Path
import numpy as np
_CO = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/C_propagator/corr_op")
sys.path.insert(0, str(_CO))
import corr_op_C_callable as prod
from canoes.sachs.sft_input.corr_op.table import load_table_2d_from_npz, make_C_fn_lambda_project
N1 = np.array([0.,0.,1.]); n2 = np.array([0.,0.,1.])
D = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix/dense")
fn3 = make_C_fn_lambda_project(load_table_2d_from_npz(str(D/"corr_op_dense_calib3.npz")), output_h_power=0.0)
fn5 = make_C_fn_lambda_project(load_table_2d_from_npz(str(D/"corr_op_dense_cellA5.npz")), output_h_power=0.0)
fn5b = make_C_fn_lambda_project(load_table_2d_from_npz(str(D/"corr_op_dense_cellA5_nchi512.npz")), output_h_power=0.0)
def c(fn,x): return float(fn(N1,float(x),n2,float(x))[0,0])
print("independent-build agreement at shared node 1373.5:")
print(f"  3-node {c(fn3,1373.5):.6e}   5-node {c(fn5,1373.5):.6e}   rel {c(fn3,1373.5)/c(fn5,1373.5)-1:+.4%}")
print(f"  n_chi 256 vs 512 at 1336.75: {c(fn5,1336.75):.6e} vs {c(fn5b,1336.75):.6e} rel {c(fn5b,1336.75)/c(fn5,1336.75)-1:+.4%}")
print("\nmid-cell error vs the 5-node truth, at lam=1336.75 (the midpoint of [1300,1373.5]):")
truth = c(fn5,1336.75)
e_prod = c(prod.C_fn.__self__ if False else prod, 1336.75) if False else None
p = float(prod.C_fn(N1,1336.75,n2,1336.75)[0,0])
t3 = c(fn3,1336.75)
print(f"  h=147.0 (production 21pt)  C={p:.6e}  err {p/truth-1:+.4%}")
print(f"  h= 73.5 (3-node rebuild)   C={t3:.6e}  err {t3/truth-1:+.4%}")
print(f"  error ratio on halving h: {abs(p/truth-1)/abs(t3/truth-1):.2f}x   (O(h)=2.0, O(h^2)=4.0)")
