"""Boundary validation: extreme gamma, extreme lambda, the narrowest cells, and
both one-sided limits AT every table node.  Plus the mechanism of R1's residual."""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
import sigma2_repaired as R
import candidates as K

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
L = K.TABLE_LAM
r5 = K.R5Exact(background=bg, apply_c0=False)

# ---- A. extreme gamma x extreme lambda ------------------------------------
print("=== A. max |scheme/R5 - 1| over lambda, per gamma (incl. extremes) ===")
GAMS = [0.5, 1.0, 5.0, 17.3, 114.3, 596.9, 1800.0]
xs = np.sort(np.concatenate([
    np.linspace(406.05, 2327.95, 500),
    np.array([406.0+d for d in (0.01,0.1,1,3)]),
    np.array([2328.0-d for d in (0.01,0.1,1,3)]),
    # inside the four narrowest cells (2300-2309-2314-2318-2328)
    np.linspace(2300.1, 2317.9, 60),
]))
xs = xs[np.min(np.abs(xs[:,None]-L[None,:]), axis=1) > 1e-7]
SCH = [("stock", ds.Sigma2Builder), ("R1 knots", R.R1Knots),
       ("R2 FIXED", K.R2Fixed), ("R7 akima", K.R7Akima), ("R8 hermite", K.R8HermiteNodes)]
print(f"{'gamma':>8} " + " ".join(f"{n:>22}" for n,_ in SCH))
print(f"{'':>8} " + " ".join(f"{'max %':>10}{'median %':>12}" for _ in SCH))
for gam in GAMS:
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    ref = np.array([r5.matrix(cg, float(x)) for x in xs])
    sc = np.maximum(np.abs(ref).max(axis=(1,2)), 1e-300)
    cells = []
    for n, cls in SCH:
        b = cls(background=bg, apply_c0=False)
        got = np.array([b.matrix(cg, float(x)) for x in xs])
        e = np.abs(got-ref).max(axis=(1,2))/sc
        cells.append(f"{100*e.max():>10.3f}{100*np.median(e):>12.2e}")
    print(f"{gam:>8.1f} " + " ".join(cells))

# ---- B. AT the nodes, both sides ------------------------------------------
print("\n=== B. exactly at interior table nodes, gamma=1', Sigma00 ===")
print("  the true density is TWO-VALUED there; a scheme can at best pick a side")
cg = float(np.cos(np.deg2rad(1.0/60.0)))
b1 = R.R1Knots(background=bg, apply_c0=False)
bs = ds.Sigma2Builder(background=bg, apply_c0=False)
print(f"{'node':>7} {'R5 left':>13} {'R5 right':>13} {'R1 knots':>13} {'stock':>13}"
      f" {'R1/left':>9} {'R1/right':>9} {'stock/left':>11}")
for nd in L[1:-1]:
    l_ = r5.matrix_side(cg, float(nd), "left")[0,0]
    r_ = r5.matrix_side(cg, float(nd), "right")[0,0]
    v1 = b1.matrix(cg, float(nd))[0,0]
    vs = bs.matrix(cg, float(nd))[0,0]
    print(f"{nd:>7.0f} {l_:>13.5e} {r_:>13.5e} {v1:>13.5e} {vs:>13.5e}"
          f" {v1/l_:>9.4f} {v1/r_:>9.4f} {vs/l_:>11.4f}")

# ---- C. mechanism of R1's residual: it is the D^4 factor, not C -----------
print("\n=== C. why R1 is not exact: it splines D^4 C, and D^4 is not polynomial ===")
n2 = ds._n2_of_cos(cg)
for lo, hi in [(1300.0, 1447.0), (2314.0, 2318.0), (406.0, 550.0)]:
    n = max(6, int(np.ceil((hi-lo)/12.0))+1)
    g = np.linspace(lo, hi, n)
    Cg = np.array([co.C_fn(ds._N1, float(t), n2, float(t))[0,0] for t in g])
    D4 = np.asarray(bg.D(g), float)**4
    sC, sF = CubicSpline(g, Cg), CubicSpline(g, D4*Cg)
    q = np.linspace(lo+1e-6, hi-1e-6, 41)
    trueC = np.array([co.C_fn(ds._N1, float(t), n2, float(t))[0,0] for t in q])
    eC = np.max(np.abs(sC(q)/trueC-1))
    d4q = np.asarray(bg.D(q), float)**4
    eF = np.max(np.abs(sF(q)/(d4q*trueC)-1))
    print(f"  cell [{lo:6.0f},{hi:6.0f}]  n={n:2d}:  spline of C  max rel err {eC:.3e}"
          f"   |  spline of D^4 C  max rel err {eF:.3e}")

# ---- D. cost in the real MC pattern ---------------------------------------
print("\n=== D. cost of one production precompute (N=1000 nodes, 2 cos) ===")
grid = np.linspace(406.0, bgm.LAM_SOURCE_BASELINE, 1000)
for n, cls in [("stock", ds.Sigma2Builder), ("R1 knots", R.R1Knots),
               ("R2 FIXED", K.R2Fixed), ("R5 EXACT", K.R5Exact)]:
    b = cls(background=bg, apply_c0=False)
    t0 = time.time()
    for x in grid:
        core._assemble_6x6(b, cg, float(x))
    print(f"  {n:<10} {time.time()-t0:>7.3f} s")
