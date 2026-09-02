"""V6: does R4's diagonal-aware C degrade the gamma-SHAPE agreement with
analysis-3 O0 (which is folded from the SAME bilinear table)?"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/psd_fix")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import sigma2_repaired as R
from sigma2_repaired import TABLE_LAM
from t12_diag_aware_builder import R4DiagAware
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE+5.0)
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"), allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"],float), np.asarray(d["o0"],float)
GAM = g_a3[(g_a3 >= 0.5) & (g_a3 <= 200.0)]
print("gamma grid (analysis-3 native):", np.round(GAM,3))
LF, LO = bgm.LAM_SOURCE_BASELINE, 406.0
br = np.concatenate([[LO], TABLE_LAM[(TABLE_LAM>LO+1e-9)&(TABLE_LAM<LF-1e-9)], [LF]])
x,w = leggauss(48); xs,ws = [],[]
for lo,hi in zip(br[:-1],br[1:]):
    xs.append(0.5*(hi-lo)*x+0.5*(hi+lo)); ws.append(0.5*(hi-lo)*w)
la, wt = np.concatenate(xs), np.concatenate(ws)
G2 = ds.order0_window(bg, la, LF, n_gauss=256)**2
BLD = {"stock": ds.Sigma2Builder(background=bg,apply_c0=False),
       "R1": R.R1Knots(background=bg,apply_c0=False),
       "R4": R4DiagAware(background=bg,apply_c0=False)}
rat = {}
for k,b in BLD.items():
    r=[]
    for gam in GAM:
        cg=float(np.cos(np.deg2rad(gam/60.0)))
        s=np.array([b.matrix(cg,float(l))[0,0] for l in la])
        r.append(float(np.sum(wt*s*G2))/float(np.interp(gam,g_a3,o0_a3)))
    rat[k]=np.array(r); print(f"  {k} done",flush=True)
print(f"\n{'builder':>8} {'median':>9} {'std(r/med)':>11} {'p2p(r/med)':>11}")
for k,r in rat.items():
    n=r/np.median(r)
    print(f"{k:>8} {np.median(r):>9.5f} {np.std(n):>11.4f} {n.max()-n.min():>11.4f}")
print("\nper-gamma normalised ratio r/median:")
for i,g in enumerate(GAM):
    print(f"  {g:>8.3f}' " + " ".join(f"{k}={rat[k][i]/np.median(rat[k]):.4f}" for k in rat))
