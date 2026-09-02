import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import corr_op_C_callable as co
n1 = ds._N1; n2 = ds._n2_of_cos(np.cos(np.deg2rad(1.0/60.0)))
t0=time.time()
for t in np.linspace(500,2300,50): co.C_fn(n1,float(t),n2,float(t))
print("scalar C_fn per call (s):", (time.time()-t0)/50)
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE+5.0)
t0=time.time(); b = ds.Sigma2Builder(background=bg, apply_c0=False); b.matrix(1.0, 1000.0)
print("stock builder cos=1 build (s):", time.time()-t0)
import sigma2_repaired as R
t0=time.time(); r = R.R1Knots(background=bg, apply_c0=False); r.matrix(1.0, 1000.0)
print("R1 builder cos=1 build (s):", time.time()-t0)
print("TABLE_LAM:", R.TABLE_LAM)
print("breaks:", R._cells(406.0, 2328.0))
