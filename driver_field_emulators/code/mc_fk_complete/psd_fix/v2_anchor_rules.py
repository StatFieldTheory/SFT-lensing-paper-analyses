"""V3: reproduce the d6 anchor, then re-integrate the SAME densities with
jump-aware and uniform rules.  Question: is the +1.61% physics or quadrature?"""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R
from sigma2_repaired import TABLE_LAM

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),
            allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])
LF = bgm.LAM_SOURCE_BASELINE
LO = 406.0

def nodes_gl(n):
    x, w = leggauss(int(n))
    return 0.5*(LF-LO)*x + 0.5*(LF+LO), 0.5*(LF-LO)*w

def nodes_panel(npp):
    br = np.concatenate([[LO], TABLE_LAM[(TABLE_LAM > LO+1e-9) & (TABLE_LAM < LF-1e-9)], [LF]])
    x, w = leggauss(int(npp)); xs, ws = [], []
    for lo, hi in zip(br[:-1], br[1:]):
        xs.append(0.5*(hi-lo)*x + 0.5*(hi+lo)); ws.append(0.5*(hi-lo)*w)
    return np.concatenate(xs), np.concatenate(ws)

def nodes_trap(N):
    la = np.linspace(LO, LF, int(N)); h = la[1]-la[0]
    w = np.full(N, h); w[0] *= 0.5; w[-1] *= 0.5
    return la, w

RULES = {"GL48": nodes_gl(48), "GL256": nodes_gl(256), "GL512": nodes_gl(512),
         "GL2048": nodes_gl(2048), "PANEL24": nodes_panel(24),
         "PANEL48": nodes_panel(48), "PANEL96": nodes_panel(96),
         "TRAP1000": nodes_trap(1000), "TRAP4000": nodes_trap(4000)}

WIN = {k: ds.order0_window(bg, la, LF, n_gauss=256)**2 for k, (la, w) in RULES.items()}
BLD = {"stock": ds.Sigma2Builder(background=bg, apply_c0=False),
       "R1":    R.R1Knots(background=bg, apply_c0=False)}

res = {}
for bk, b in BLD.items():
    for rk, (la, w) in RULES.items():
        meds = []
        for gam in GAM:
            cosg = float(np.cos(np.deg2rad(gam/60.0)))
            s = np.array([b.matrix(cosg, float(l))[0, 0] for l in la])
            I = float(np.sum(w * s * WIN[rk]))
            meds.append(I / float(np.interp(gam, g_a3, o0_a3)))
        res[(bk, rk)] = (float(np.median(meds)), np.array(meds))
        print(f"  {bk:>6} {rk:>9}  median {res[(bk,rk)][0]:.6f}", flush=True)

print(f"\n{'rule':>10} {'stock':>10} {'R1':>10} {'R1/stock-1':>12}")
for rk in RULES:
    a, bb = res[("stock", rk)][0], res[("R1", rk)][0]
    print(f"{rk:>10} {a:>10.6f} {bb:>10.6f} {bb/a-1:>+11.4%}")

print("\nper-gamma R1/stock under GL256 vs PANEL48:")
for i, gam in enumerate(GAM):
    r256 = res[("R1","GL256")][1][i]/res[("stock","GL256")][1][i]
    rp   = res[("R1","PANEL48")][1][i]/res[("stock","PANEL48")][1][i]
    print(f"  {gam:>7.2f}'  GL256 {r256:+.5%}   PANEL48 {rp:+.5%}"
          .replace("+1","+1").replace("GL256 1","GL256 "))
