"""Is the ~1% residual of A2.4 in the DISCRETE kernel or in driver_stats.order0_mc's
48-node Gauss-Legendre quadrature?  Refine the quadrature and see which one moves."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
DS, BG, CORE = ex._DS, ex._BG, ex._CORE

ANC = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
       "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")
ad = np.load(ANC, allow_pickle=True)
m = (ad["a"] == 0) & (ad["b"] == 0) & (ad["order"] == 0)
gam_a, cos_a = [], []
for xi, yi in zip(ad["x"][m], ad["y"][m]):
    xi = np.asarray(xi, float); yi = np.asarray(yi, float)
    c = float(np.clip(np.dot(xi/np.linalg.norm(xi), yi/np.linalg.norm(yi)), -1, 1))
    gam_a.append(math.degrees(math.acos(c))*60.0); cos_a.append(c)
o = np.argsort(gam_a); gam_a = np.asarray(gam_a)[o]; cos_a = np.asarray(cos_a)[o]
vs = np.asarray(ad["value"][m], float)[o]
tfin = float(ad["t_final"][0])

grid4 = ex.build_grid(n_lambda=4000, apply_anchor=True)
b = grid4.builder
print(f"{'gamma':>8} {'ng=48':>13} {'ng=96':>13} {'ng=192':>13} {'ng=384':>13} {'ng=768':>13}")
res = {}
for gt in (1.0, 5.0, 17.3):
    i = int(np.argmin(np.abs(gam_a-gt))); cg = float(cos_a[i])
    vals = [DS.order0_mc(cg, lam_f=tfin, builder=b, n_gauss=ng) for ng in (48,96,192,384,768)]
    res[gt] = (cg, vals, i)
    print(f"{gam_a[i]:8.3f}" + "".join(f"{v:13.6e}" for v in vals))

print()
print("ratio to the ng=768 value:")
for gt, (cg, vals, i) in res.items():
    print(f"  gamma={gt:6.2f}': " + "  ".join(f"{v/vals[-1]:.5f}" for v in vals))

print()
print("Composite-panel Gauss-Legendre reference (32 panels x 64 nodes) of the SAME formula,")
print("and the discrete white-noise Order-0 at N=4000 / Richardson, against it:")
from numpy.polynomial.legendre import leggauss
nod, wt = leggauss(64)
edges = np.linspace(406.0, tfin, 33)


def o0_panel(cg):
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        la = 0.5*(hi-lo)*nod + 0.5*(hi+lo)
        jw = 0.5*(hi-lo)*wt
        sig = np.array([b.matrix(cg, float(l))[0, 0] for l in la])
        G = DS.order0_window(b.bg, la, tfin, n_gauss=192)
        tot += float(np.sum(jw*sig*G**2))
    return tot


def o0_white(grid, cg):
    S2 = np.array([CORE._assemble_6x6(grid.builder, cg, float(l)) for l in grid.lam])
    return np.einsum("j,jAB->AB", grid.dlam*grid.Wd**2, S2, optimize=True)[0, 3]


grid2 = ex.build_grid(n_lambda=2000, apply_anchor=True)
for gt, (cg, vals, i) in res.items():
    ref = o0_panel(cg)
    w4 = o0_white(grid4, cg); w2 = o0_white(grid2, cg)
    rich = 2*w4 - w2
    print(f"  gamma={gam_a[i]:7.3f}'  panel_ref={ref:.6e}  ng48/ref={vals[0]/ref:.5f}"
          f"  ng768/ref={vals[-1]/ref:.5f}"
          f"  disc(N=4000)/ref={w4/ref:.5f}  Richardson/ref={rich/ref:.5f}"
          f"  analysis3/ref={vs[i]/ref:.5f}")
