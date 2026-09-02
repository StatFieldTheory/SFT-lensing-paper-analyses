"""Does a MULTI-POINT batch call perturb the scalar C_fn at the same cos?"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
N1 = np.array([0.0, 0.0, 1.0])
def n2c(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1-c*c)))
    return np.array([s, 0.0, c])
TS = np.linspace(420.0, 2320.0, 180)
for gam in (0.5, 1.0, 5.0, 17.3, 2.0, 10.0):
    cg = float(np.cos(np.deg2rad(gam/60.0))); n2 = n2c(cg)
    before = np.array([float(corr_op.C_fn(N1, float(t), n2, float(t))[0, 0])
                       for t in (700.5, 1234.5, 2100.7)])
    B = corr_op.C_fn_batch(N1, TS, n2, TS)                 # the 180-point batch
    after = np.array([float(corr_op.C_fn(N1, float(t), n2, float(t))[0, 0])
                      for t in (700.5, 1234.5, 2100.7)])
    bat = np.array([float(corr_op.C_fn_batch(N1, np.array([t]), n2, np.array([t]))[0, 0, 0])
                    for t in (700.5, 1234.5, 2100.7)])
    print(f"gamma={gam:6.2f}'  max|after/before-1| = {np.abs(after/before-1).max():.3e}   "
          f"max|batch/scalar-1| = {np.abs(bat/before-1).max():.3e}")

print("\nsame, but comparing the 180-point batch itself to the scalar loop:")
for gam in (0.5, 1.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam/60.0))); n2 = n2c(cg)
    B = corr_op.C_fn_batch(N1, TS, n2, TS)[:, 0, 0]
    S = np.array([float(corr_op.C_fn(N1, float(t), n2, float(t))[0, 0]) for t in TS])
    print(f"  gamma={gam:6.2f}'  max rel |batch-scalar| = {np.abs(B/S-1).max():.3e}")
