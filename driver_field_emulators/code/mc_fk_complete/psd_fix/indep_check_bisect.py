import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, assemble6, _bg, _ds
import sigma2_repaired as R
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
ev = IndepDensity(background=bgd, apply_c0=False)
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
r1 = R.R1Knots(background=bgd, apply_c0=False)
r0 = R.R0Reference(background=bgd, apply_c0=False)
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"), allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5,1.0,2.0,5.0,10.0,17.3,44.4,114.3,293.9,596.9])
nod, wt = leggauss(256)
la = 0.5*(LAMF-406.0)*nod + 0.5*(LAMF+406.0); w = 0.5*(LAMF-406.0)*wt
G2 = _ds.order0_window(bgd, la, LAMF, n_gauss=64)**2
def med(b):
    r = []
    for g in GAM:
        cg = float(np.cos(np.deg2rad(float(g)/60.0)))
        s = np.array([b.matrix(cg, float(x))[0,0] for x in la])
        r.append(float(np.sum(w*s*G2))/float(np.interp(g, g_a3, o0_a3)))
    return float(np.median(r))
print("A) before anything else            :", f"{med(stock):.5f}")
_ = [assemble6(ev, float(np.cos(np.deg2rad(0.5/60))), float(l)) for l in np.linspace(410,2310,50)]
print("B) after using ev/assemble6        :", f"{med(stock):.5f}")
LAM = np.linspace(410.0, 2310.0, 600)
_ = [stock.matrix(1.0, float(l)) for l in LAM]
print("C) after stock.matrix(1.0, 600 lam):", f"{med(stock):.5f}")
cg = float(np.cos(np.deg2rad(17.3/60.0)))
_ = [stock.matrix(cg, float(l)) for l in LAM]
print("D) after stock.matrix(cg17.3, ...) :", f"{med(stock):.5f}")
_ = [r1.matrix(1.0, float(l)) for l in LAM]; _ = [r0.matrix(1.0, float(l)) for l in LAM[:50]]
print("E) after r1/r0 use                 :", f"{med(stock):.5f}")
print("   stock cache keys:", sorted(stock._cache.keys()))
