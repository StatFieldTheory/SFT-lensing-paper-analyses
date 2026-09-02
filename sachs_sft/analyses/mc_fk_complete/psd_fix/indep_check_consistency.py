import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
from indep_density import IndepDensity, NODES, _bg, _ds
bgd = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
LAMF = _bg.LAM_SOURCE_BASELINE
stock = _ds.Sigma2Builder(background=bgd, apply_c0=False)
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"), allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5,1.0,2.0,5.0,10.0,17.3,44.4,114.3,293.9,596.9])
nod, wt = leggauss(256)
la = 0.5*(LAMF-406.0)*nod + 0.5*(LAMF+406.0); w = 0.5*(LAMF-406.0)*wt
for nwin in (48, 64, 128, 256, 512):
    G2 = _ds.order0_window(bgd, la, LAMF, n_gauss=nwin)**2
    rr = []
    for g in GAM:
        cg = float(np.cos(np.deg2rad(g/60.0)))
        s = np.array([stock.matrix(cg, float(l))[0,0] for l in la])
        rr.append(float(np.sum(w*s*G2))/float(np.interp(g, g_a3, o0_a3)))
    print(f"  window n_gauss={nwin:4d}: median ratio = {np.median(rr):.5f}")
print("  ds.order0_mc(n_gauss=256) median:",
      f"{np.median([float(_ds.order0_mc(float(np.cos(np.deg2rad(g/60.0))), lam_f=LAMF, builder=stock, n_gauss=256))/float(np.interp(g,g_a3,o0_a3)) for g in GAM]):.5f}")
