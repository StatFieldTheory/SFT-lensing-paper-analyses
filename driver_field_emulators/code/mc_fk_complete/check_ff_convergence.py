"""Is the FF Monte-Carlo converged in the lambda step?

The FF vertex is a DRIFT term in the Euler recursion, so its discretisation
error is O(dlam), unlike the noise term.  The published markers used
n_lambda = 4000; this asks whether that is converged.
"""
import time
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from pathlib import Path

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3 = np.asarray(d["g"], float); ff_a3 = np.asarray(d["ff_mom"], float)
gm = np.asarray(d["g_mc"], float); ffc = np.asarray(d["ffc_mc"], float)

GAM = 1.0
ana = float(np.interp(GAM, g_a3, ff_a3))
pub = float(ffc[np.argmin(np.abs(gm - GAM))])
print(f"gamma = {GAM}'   analytic FF = {ana:.4e}   published marker = {pub:.4e} "
      f"(ratio {pub/ana:.4f})")
print(f"\n{'n_lambda':>9} {'dlam[Mpc]':>10} {'FF (anchor OFF)':>17} {'+/-':>10} "
      f"{'/analytic':>10} {'sec':>6}")
prev = None
for nl in (500, 1000, 2000, 4000, 8000):
    t0 = time.time(); acc = []
    for s in range(4):
        cfg = mc.MCConfig(n_real=12000, batch_size=3000, n_lambda=nl,
                          use_f_vertex=True, apply_anchor=False, seed=11 + 7 * s)
        acc.append(mc.simulate_ff_crn(cfg, GAM).ff_moment[0, 0])
    v = np.array(acc); m = v.mean(); se = v.std(ddof=1) / np.sqrt(v.size)
    dl = (bgm.LAM_SOURCE_BASELINE - 406.0) / (nl - 1)
    print(f"{nl:>9} {dl:>10.3f} {m:>17.4e} {se:>10.1e} {m/ana:>10.4f} "
          f"{time.time()-t0:>6.0f}", flush=True)
    prev = m
