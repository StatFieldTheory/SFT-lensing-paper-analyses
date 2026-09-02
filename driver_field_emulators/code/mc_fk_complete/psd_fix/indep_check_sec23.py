"""Bisect: which call in section 2/3 changes stock's Order-0?"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, _bg, _ds
import sigma2_repaired as R
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)
r0 = R.R0Reference(background=bgd, apply_c0=False)
GAM = [0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3]
def vec():
    return np.array([float(_ds.order0_mc(float(np.cos(np.deg2rad(g/60.0))),
        lam_f=LAMF, builder=stock, n_gauss=256)) for g in GAM])
v0 = vec(); print("baseline:", " ".join(f"{x:.6e}" for x in v0))
LAM = np.linspace(410.0, 2310.0, 600)
steps = [("ev.matrix", lambda cc: [ev.matrix(cc, float(l)) for l in LAM]),
         ("stock.matrix", lambda cc: [stock.matrix(cc, float(l)) for l in LAM]),
         ("r1.matrix", lambda cc: [r1.matrix(cc, float(l)) for l in LAM]),
         ("r0.matrix", lambda cc: [r0.matrix(cc, float(l)) for l in LAM])]
for gam in (1.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    for cc, blab in ((1.0, "within"), (cg, "cross")):
        for nm, fn in steps:
            fn(cc)
            v = vec(); d = np.abs(v/v0 - 1)
            if d.max() > 1e-12:
                print(f"  CHANGED after {nm}(cos={cc:.9f}, gamma={gam}, {blab}): "
                      f"max rel {d.max():.3e} at gamma={GAM[int(np.argmax(d))]}")
                print("     ", " ".join(f"{x:+.2e}" for x in (v/v0-1)))
                v0 = v
print("done")
