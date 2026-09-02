"""A8 -- is the (D_la/D_t)^2 response exponent testable against something OUTSIDE
this folder, or is it an assumption shared by everything?

The equal-time collapse O0(gamma) = int Sigma2_00(cos;la) G_p(la)^2 dla with
G_p(la) = int_la^lam_f (D_la/D_t)^p dt is compared to the sft-wick analysis-3
C_corr_op_O0 sweep.  Sigma2_00(cos; la) has a genuinely gamma-DEPENDENT profile in
la, so a wrong exponent cannot be absorbed into one constant: it must show up as a
gamma-dependent ratio.  The published anchor is the claim that the ratio is a single
constant (0.881) over 0.5'-150'.  Repeat that with p = 1, 2, 3, 4 and see which p
gives a gamma-independent ratio."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from numpy.polynomial.legendre import leggauss
import fk_expect_exact as ex
DS, BG = ex._DS, ex._BG

ANC = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
       "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")
ad = np.load(ANC, allow_pickle=True)
m = (ad["a"] == 0) & (ad["b"] == 0) & (ad["order"] == 0)
gam, cos = [], []
for x, y in zip(ad["x"][m], ad["y"][m]):
    x = np.asarray(x, float); y = np.asarray(y, float)
    c = float(np.clip(np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1))
    gam.append(math.degrees(math.acos(c))*60.0); cos.append(c)
o = np.argsort(gam)
gam = np.asarray(gam)[o]; cos = np.asarray(cos)[o]
vs = np.asarray(ad["value"][m], float)[o]
tfin = float(ad["t_final"][0])
grid = ex.build_grid(n_lambda=200, apply_anchor=False)   # raw Sigma2, no anchor
b = grid.builder
bg = b.bg
nod, wt = leggauss(64)
edges = np.linspace(406.0, tfin, 25)


def window(la, p, ng=192):
    n2, w2 = leggauss(ng)
    D_la = np.asarray(bg.D(la), float)
    out = np.empty_like(np.asarray(la, float))
    for i, l in enumerate(np.atleast_1d(la)):
        t = 0.5*(tfin-l)*n2 + 0.5*(tfin+l)
        jw = 0.5*(tfin-l)*w2
        out[i] = float(np.sum(jw * (D_la[i]/np.asarray(bg.D(t), float))**p))
    return out


def o0_p(cg, p):
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        la = 0.5*(hi-lo)*nod + 0.5*(hi+lo)
        jw = 0.5*(hi-lo)*wt
        sig = np.array([b.matrix(cg, float(l))[0, 0] for l in la])
        tot += float(np.sum(jw*sig*window(la, p)**2))
    return tot


sel = [i for i in range(len(gam)) if gam[i] < 200.0 and vs[i] > 0]
sel = sel[::max(1, len(sel)//10)]
print(f"{'gamma':>10} {'analysis3':>13}" + "".join(f"{'p='+str(p):>12}" for p in (1,2,3,4)))
ratios = {p: [] for p in (1, 2, 3, 4)}
for i in sel:
    row = []
    for p in (1, 2, 3, 4):
        r = o0_p(float(cos[i]), p) / vs[i]
        ratios[p].append(r); row.append(r)
    print(f"{gam[i]:10.3f} {vs[i]:13.5e}" + "".join(f"{r:12.4f}" for r in row))
print()
print(f"{'exponent p':>12} {'median ratio':>14} {'rel. spread (max/min - 1)':>28}")
for p in (1, 2, 3, 4):
    r = np.array(ratios[p])
    print(f"{p:>12} {np.median(r):14.4f} {r.max()/r.min()-1:28.4f}")
print("\nThe published anchor is p = 2 with a single constant 0.881; a wrong p must")
print("show as a gamma-dependent ratio, i.e. a large spread.")
