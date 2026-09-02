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
la = 0.5*(LAMF-406)*nod+0.5*(LAMF+406); w = 0.5*(LAMF-406)*wt
G2 = _ds.order0_window(bgd, la, LAMF, n_gauss=64)**2
b = stock
A = [float(np.sum(w*np.array([b.matrix(float(np.cos(np.deg2rad(g/60))), float(l))[0,0]
     for l in la])*G2))/float(np.interp(g, g_a3, o0_a3)) for g in GAM]      # summary form
B = [float(_ds.order0_mc(float(np.cos(np.deg2rad(g/60.0))), lam_f=LAMF,
     builder=b, n_gauss=256))/float(np.interp(g, g_a3, o0_a3)) for g in GAM]  # canonical
print(f"{'gamma':>8} {'summary form':>14} {'canonical':>14} {'ratio':>10}")
for g, a, c in zip(GAM, A, B):
    print(f"{g:>8.2f} {a:>14.5f} {c:>14.5f} {a/c:>10.6f}")
print("medians:", f"{np.median(A):.5f}", f"{np.median(B):.5f}")
