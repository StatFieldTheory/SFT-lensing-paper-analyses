"""V5: does the diagonal-aware builder move the Order-0 anchor by ~+6%?
Same harness as v2, adding R4DiagAware.  Plus the R0Reference endpoint bug."""
import sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/psd_fix")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import sigma2_repaired as R
from sigma2_repaired import TABLE_LAM
from t12_diag_aware_builder import R4DiagAware

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"), allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)
GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3])
LF, LO = bgm.LAM_SOURCE_BASELINE, 406.0

def gl(n):
    x, w = leggauss(int(n)); return 0.5*(LF-LO)*x+0.5*(LF+LO), 0.5*(LF-LO)*w
def panel(npp):
    br = np.concatenate([[LO], TABLE_LAM[(TABLE_LAM>LO+1e-9)&(TABLE_LAM<LF-1e-9)], [LF]])
    x, w = leggauss(int(npp)); xs, ws = [], []
    for lo, hi in zip(br[:-1], br[1:]):
        xs.append(0.5*(hi-lo)*x+0.5*(hi+lo)); ws.append(0.5*(hi-lo)*w)
    return np.concatenate(xs), np.concatenate(ws)
def trap(N):
    la = np.linspace(LO, LF, int(N)); h = la[1]-la[0]
    w = np.full(N, h); w[0]*=0.5; w[-1]*=0.5; return la, w

RULES = {"GL256": gl(256), "PANEL48": panel(48), "TRAP4000": trap(4000)}
WIN = {k: ds.order0_window(bg, la, LF, n_gauss=256)**2 for k,(la,w) in RULES.items()}
BLD = {"stock": ds.Sigma2Builder(background=bg, apply_c0=False),
       "R1": R.R1Knots(background=bg, apply_c0=False),
       "R4": R4DiagAware(background=bg, apply_c0=False)}
out = {}
for bk,b in BLD.items():
    for rk,(la,w) in RULES.items():
        m=[]
        for gam in GAM:
            cg=float(np.cos(np.deg2rad(gam/60.0)))
            s=np.array([b.matrix(cg,float(l))[0,0] for l in la])
            m.append(float(np.sum(w*s*WIN[rk]))/float(np.interp(gam,g_a3,o0_a3)))
        out[(bk,rk)]=(float(np.median(m)),np.array(m))
        print(f"  {bk:>5} {rk:>9} median {out[(bk,rk)][0]:.6f}",flush=True)
print(f"\n{'rule':>9} {'stock':>9} {'R1':>9} {'R4':>9} | {'R1/stock':>10} {'R4/stock':>10}")
for rk in RULES:
    a,b1,b4 = out[("stock",rk)][0], out[("R1",rk)][0], out[("R4",rk)][0]
    print(f"{rk:>9} {a:>9.5f} {b1:>9.5f} {b4:>9.5f} | {b1/a-1:>+9.4%} {b4/a-1:>+9.4%}")
print("\nR4/stock per gamma (PANEL48):")
for i,g in enumerate(GAM):
    print(f"  {g:>7.2f}'  {out[('R4','PANEL48')][1][i]/out[('stock','PANEL48')][1][i]-1:+.4%}")

# --- R0Reference endpoint bug ---
r0 = R.R0Reference(background=bg, apply_c0=False)
print("\n[R0 endpoint bug] Sigma00 at lam = lam_lo and just inside:")
for x in (406.0, 406.0+1e-6, 407.0, 410.0):
    print(f"  lam={x:<12.6f} R0 {r0.matrix(1.0,x)[0,0]:+.6e}   R1 {BLD['R1'].matrix(1.0,x)[0,0]:+.6e}")
