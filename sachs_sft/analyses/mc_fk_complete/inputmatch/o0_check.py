"""Order-0: the MC's driving variance against the workflow's, at the FF gammas.

Order-0 is LINEAR in Sigma2; FF is QUADRATIC.  This fixes the linear ratio c so
the FF prediction c^2 can be tested.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g = np.asarray(d["g"], float); o0 = np.asarray(d["o0"], float)
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
print(f"ANCHOR_C0 = {ds.ANCHOR_C0:.6f}   ANCHOR_C0^2 = {ds.ANCHOR_C0**2:.6f}")
print(f"{'gamma':>8} {'O0_mc(OFF)':>13} {'O0_workflow':>13} {'c = ratio':>10} "
      f"{'c^2':>8} {'c*anchor':>9} {'(c*anc)^2':>10}")
for gam in (1.0, 2.6, 6.7, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    a3 = float(np.interp(gam, g, o0))
    b = ds.Sigma2Builder(background=bg, apply_c0=False)
    v = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b, n_gauss=256))
    c = v / a3
    print(f"{gam:>8.2f} {v:>13.6e} {a3:>13.6e} {c:>10.5f} {c*c:>8.5f} "
          f"{c*ds.ANCHOR_C0:>9.5f} {(c*ds.ANCHOR_C0)**2:>10.5f}")
