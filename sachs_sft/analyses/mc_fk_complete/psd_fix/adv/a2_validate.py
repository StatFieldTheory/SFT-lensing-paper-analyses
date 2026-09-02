"""Validate the closed-form density R5Exact against an INDEPENDENT high-order
finite difference confined strictly inside each cell, plus the two assumptions
it rests on (transpose symmetry; the background theta = D'/D)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
import candidates as K

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
L = K.TABLE_LAM
N1 = ds._N1

# ---- 1. transpose symmetry C_ab(t1,t2) == C_ba(t2,t1) ---------------------
print("=== 1. transpose symmetry of the off-diagonal table entries ===")
worst = 0.0
for gam in (1.0, 17.3):
    n2 = ds._n2_of_cos(np.cos(np.deg2rad(gam/60.0)))
    for i in range(L.size-1):
        up = np.asarray(co.C_fn(N1, float(L[i]), n2, float(L[i+1])), float)
        dn = np.asarray(co.C_fn(N1, float(L[i+1]), n2, float(L[i])), float)
        sc = max(np.max(np.abs(up)), np.max(np.abs(dn)), 1e-300)
        worst = max(worst, float(np.max(np.abs(up.T - dn))/sc))
print("  max rel |C(t1,t2)^T - C(t2,t1)| over 20 cells x 2 gammas: %.3e" % worst)

# ---- 2. theta_sa vs a direct FD of the D spline ---------------------------
print("\n=== 2. theta_sa == dlnD/dlam ? ===")
xs = np.array([420., 700., 1300., 1900., 2300., 2325.])
for h in (1.0, 0.25):
    fd = (np.log(bg.D(xs+h)) - np.log(bg.D(xs-h)))/(2*h)
    th = np.asarray(bg.theta_sa(xs), float)
    sp = bg._D_spline; an = sp(xs,1)/sp(xs)
    print("  h=%4.2f  |FD/theta_sa-1| max %.3e   |FD/(D'/D)-1| max %.3e"
          % (h, np.max(np.abs(fd/th-1)), np.max(np.abs(fd/an-1))))
print("  per-x  theta_sa vs D'/D :", " ".join("%.1e"%v for v in np.abs(np.asarray(bg.theta_sa(xs))/(bg._D_spline(xs,1)/bg._D_spline(xs))-1)))

# ---- 3. R5Exact vs independent Richardson FD, strictly inside cells --------
print("\n=== 3. R5Exact vs Richardson FD of D^4 C (independent) ===")
r5 = K.R5Exact(background=bg, apply_c0=False)
rng = np.random.default_rng(7)
for gam in (1.0, 17.3):
    cg = float(np.cos(np.deg2rad(gam/60.0)))
    n2 = ds._n2_of_cos(cg)
    def f(t):
        return (float(bg.D(t))**4) * np.asarray(co.C_fn(N1, float(t), n2, float(t)), float)
    errs = []
    for i in range(L.size-1):
        lo, hi = L[i], L[i+1]; w = hi-lo
        for u in (0.2, 0.5, 0.8):
            x = lo + u*w
            h = min(0.4*min(u, 1-u)*w, 2.0)   # x +- 2h stays strictly inside
            # 5-point stencil, all nodes inside the cell
            d = (f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h)
            ref = d/float(bg.D(x))**4
            got = r5.matrix(cg, x)
            den = np.maximum(np.abs(0.5*(ref+ref.T)), 1e-300)
            errs.append(np.max(np.abs(got-0.5*(ref+ref.T))/den))
    errs = np.array(errs)
    print("  gamma=%5.2f'  n=%d  max %.3e  median %.3e  p90 %.3e"
          % (gam, errs.size, errs.max(), np.median(errs), np.percentile(errs,90)))

# ---- 4. is the density continuous at the nodes? ---------------------------
print("\n=== 4. one-sided density limits at interior table nodes (gamma=1', Sigma00) ===")
cg = float(np.cos(np.deg2rad(1.0/60.0)))
print("   node        left            right        (right/left - 1)")
for nd in L[1:-1]:
    lft = r5.matrix_side(cg, float(nd), "left")[0,0]
    rgt = r5.matrix_side(cg, float(nd), "right")[0,0]
    print("  %7.1f  %14.6e  %14.6e   %+9.2f%%" % (nd, lft, rgt, (rgt/lft-1)*100))
