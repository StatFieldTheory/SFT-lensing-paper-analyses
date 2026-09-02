"""Does calling corr_op.C_fn_batch change what the SCALAR corr_op.C_fn returns
afterwards at the same cos?  (Shared angular-channel cache with rounding.)"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

N1 = np.array([0.0, 0.0, 1.0])
def n2c(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1-c*c)))
    return np.array([s, 0.0, c])

print(f"{'gamma':>8} {'scalar BEFORE':>18} {'scalar AFTER batch':>20} {'rel diff':>11} "
      f"{'batch':>18} {'batch/scalar-1':>15}")
for gam in (0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    t = 1234.5
    a = float(corr_op.C_fn(N1, t, n2c(cg), t)[0, 0])          # scalar first
    b = float(corr_op.C_fn_batch(N1, np.array([t]), n2c(cg), np.array([t]))[0, 0, 0])
    c = float(corr_op.C_fn(N1, t, n2c(cg), t)[0, 0])          # scalar again
    print(f"{gam:>8.2f} {a:>18.10e} {c:>20.10e} {c/a-1:>11.2e} {b:>18.10e} "
          f"{b/a-1:>15.2e}")

print("\nsame test but BATCH FIRST on a fresh cos (never queried scalar before):")
for gam in (3.7, 7.9, 23.1):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    t = 1234.5
    b = float(corr_op.C_fn_batch(N1, np.array([t]), n2c(cg), np.array([t]))[0, 0, 0])
    a = float(corr_op.C_fn(N1, t, n2c(cg), t)[0, 0])
    print(f"  gamma={gam:6.2f}'  batch {b:.10e}  scalar-after {a:.10e}  rel {a/b-1:+.2e}")

print("\nOrder-0 with a builder whose C comes from scalar vs from batch:")
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
for gam in (0.5, 1.0, 5.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    b1 = ds.Sigma2Builder(background=bg, apply_c0=False)
    v1 = float(ds.order0_mc(cg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b1, n_gauss=256))
    _ = corr_op.C_fn_batch(N1, np.linspace(500, 2300, 40), n2c(cg),
                           np.linspace(500, 2300, 40))     # warm the batch path
    b2 = ds.Sigma2Builder(background=bg, apply_c0=False)    # FRESH builder
    v2 = float(ds.order0_mc(cg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b2, n_gauss=256))
    print(f"  gamma={gam:6.2f}'  O0 before batch {v1:.8e}   after batch {v2:.8e}   "
          f"rel {v2/v1-1:+.3e}")
