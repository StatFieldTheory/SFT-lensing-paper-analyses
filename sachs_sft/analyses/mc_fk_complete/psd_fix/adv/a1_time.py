import sys, time, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as CO
N1 = ds._N1; n2 = ds._n2_of_cos(float(np.cos(np.deg2rad(1.0/60))))
CO.C_fn(N1, 1000.0, n2, 1000.0)
t0=time.time()
for t in np.linspace(500,2300,200): CO.C_fn(N1, float(t), n2, float(t))
print(f"C_fn scalar: {(time.time()-t0)/200*1e3:.3f} ms/call (same cos)")
n2b = ds._n2_of_cos(float(np.cos(np.deg2rad(5.0/60))))
t0=time.time(); CO.C_fn(N1, 1000.0, n2b, 1000.0); print(f"first call new cos: {(time.time()-t0)*1e3:.1f} ms")
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE+5.0)
t0=time.time()
for _ in range(20): ds.order0_window(bg, np.linspace(410,2310,256), bgm.LAM_SOURCE_BASELINE, n_gauss=256)
print(f"order0_window(256 pts, n_gauss=256): {(time.time()-t0)/20:.3f} s")
