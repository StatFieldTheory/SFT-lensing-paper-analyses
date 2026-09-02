"""ATTACK (a), sharpened: the paper's published external anchor is
"the Order-0 comparison against analysis-3 over a 290x gamma lever arm, which
pins the response exponent to p = 2 with a SPREAD OF 0.029, against 1.29 and
1.82 for p = 1 and 3"  (PAPER_PROPOSAL s.0 / review/ADJUDICATION.md B).

That statement is a claim about the gamma-SHAPE of order0_mc / O0_analysis3.
Here it is recomputed on the analysis-3 gamma grid with the stock and repaired
builders, under the rule the code actually uses (GL n_gauss) and under a
converged rule (Gauss-Legendre panelised at the corr_op table nodes).
The p = 1 / 3 counterfactuals are rebuilt from the same Sigma2 by changing the
response exponent in the LOS window G(la) = int (D_la/D_t)^p dt.
"""
import sys, math
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
anc = np.load("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
              "/sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz",
              allow_pickle=True)
m = (anc["a"]==0)&(anc["b"]==0)&(anc["order"]==0)
gam, cos = [], []
for xi, yi in zip(anc["x"][m], anc["y"][m]):
    xi=np.asarray(xi,float); yi=np.asarray(yi,float)
    c=float(np.clip(np.dot(xi/np.linalg.norm(xi), yi/np.linalg.norm(yi)),-1,1))
    gam.append(math.degrees(math.acos(c))*60.0); cos.append(c)
o=np.argsort(gam)
gam=np.asarray(gam)[o]; cos=np.asarray(cos)[o]
vs=np.asarray(anc["value"][m],float)[o]
LAM_F = float(anc["t_final"][0]); LAM_LO = 406.0
sel = (vs > 0) & (gam < 200.0)
print(f"analysis-3 grid: {sel.sum()} usable gammas, "
      f"{gam[sel].min():.3f}' -> {gam[sel].max():.3f}'  "
      f"(lever arm {gam[sel].max()/gam[sel].min():.0f}x), lam_f={LAM_F:.3f}")

BREAKS = R._cells(LAM_LO, LAM_F)
def Gp(la, p):
    nd, wt = leggauss(96); Dla = np.asarray(bg.D(la), float)
    out = np.empty_like(np.asarray(la, float))
    for i, l in enumerate(np.atleast_1d(la)):
        t = 0.5*(LAM_F-l)*nd + 0.5*(LAM_F+l); jw = 0.5*(LAM_F-l)*wt
        out[i] = float(np.sum(jw*(Dla[i]/np.asarray(bg.D(t)))**p))
    return out

def o0_gl(b, c, p, n):
    nd, wt = leggauss(n)
    la = 0.5*(LAM_F-LAM_LO)*nd + 0.5*(LAM_F+LAM_LO); jw = 0.5*(LAM_F-LAM_LO)*wt
    s = np.array([b.matrix(c, float(l))[0,0] for l in la])
    return float(np.sum(jw*s*Gp(la,p)**2))

def o0_panel(b, c, p, npt=64):
    nd, wt = leggauss(npt); tot = 0.0
    for lo, hi in zip(BREAKS[:-1], BREAKS[1:]):
        la = 0.5*(hi-lo)*nd + 0.5*(hi+lo); jw = 0.5*(hi-lo)*wt
        s = np.array([b.matrix(c, float(l))[0,0] for l in la])
        tot += float(np.sum(jw*s*Gp(la,p)**2))
    return tot

def stats(r):
    r = r/np.median(r)
    return r.max()-r.min(), float(np.std(r))

print(f"\n{'builder':>7} {'rule':>8} {'p':>3} {'median ratio':>13} "
      f"{'p-t-p spread':>13} {'std':>8}")
for bname, cls in (("stock", ds.Sigma2Builder), ("R1", R.R1Knots)):
    b = cls(background=bg, apply_c0=False)
    for rule, fn in (("GL48", lambda bb,c,p: o0_gl(bb,c,p,48)),
                     ("GL256", lambda bb,c,p: o0_gl(bb,c,p,256)),
                     ("panel", o0_panel)):
        for p in (2, 1, 3):
            o0 = np.array([fn(b, float(c), p) for c in cos[sel]])
            ratio = o0/vs[sel]
            ptp, sd = stats(ratio)
            print(f"{bname:>7} {rule:>8} {p:>3} {np.median(ratio):>13.5f} "
                  f"{ptp:>13.4f} {sd:>8.4f}")
