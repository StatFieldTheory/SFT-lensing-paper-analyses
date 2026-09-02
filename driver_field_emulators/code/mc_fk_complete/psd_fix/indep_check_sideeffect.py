import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, _bg, _ds
import sigma2_repaired as R
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"), allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5,1.0,2.0,5.0,10.0,17.3,44.4,114.3,293.9,596.9])

def med(b):
    return float(np.median([float(_ds.order0_mc(float(np.cos(np.deg2rad(g/60.0))),
        lam_f=LAMF, builder=b, n_gauss=256))/float(np.interp(g,g_a3,o0_a3)) for g in GAM]))

s1 = _ds.Sigma2Builder(background=bgd, apply_c0=False)
print("fresh stock builder                     :", f"{med(s1):.5f}")
r0 = R.R0Reference(background=bgd, apply_c0=False)
_ = [r0.matrix(1.0, float(l)) for l in np.linspace(410, 2310, 50)]
print("stock builder AFTER using R0Reference   :", f"{med(s1):.5f}")
s2 = _ds.Sigma2Builder(background=bgd, apply_c0=False)
print("NEW stock builder after R0Reference     :", f"{med(s2):.5f}")
print("R0Reference itself                      :", f"{med(r0):.5f}")
print("bg identity check: r0.bg is bgd ->", r0.bg is bgd,
      " ; D(1300) =", float(bgd.D(1300.0)))
