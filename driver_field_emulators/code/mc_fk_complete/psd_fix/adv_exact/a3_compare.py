"""Head-to-head: stock, R0, R1, R2, R2Fixed, R7Akima, R8Hermite, R5Exact.

Reference = R5Exact (closed form; validated against an independent Richardson
FD in a2 and by the h-convergence block below).
Metrics: max/median relative error on the density, PSD violations on the 6x6,
corr_op evaluations at construction, seconds per matrix() call.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
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
N1 = ds._N1

# ---- 0. FD -> R5 convergence (proves R5 is the exact one, not the FD) ------
print("=== 0. 5-point FD of D^4 C converges to R5Exact as h^4 ===")
r5 = K.R5Exact(background=bg, apply_c0=False)
cg1 = float(np.cos(np.deg2rad(1.0/60.0))); n2 = ds._n2_of_cos(cg1)
def f(t):
    return (float(bg.D(t))**4)*np.asarray(co.C_fn(N1, float(t), n2, float(t)), float)
x = 1373.5   # mid-cell of 1300-1447
ex = r5.matrix(cg1, x)[0,0]
print("   R5Exact Sigma00(1373.5) = %.12e" % ex)
prev=None
for h in (8.0, 4.0, 2.0, 1.0, 0.5, 0.25):
    d = (f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h)/float(bg.D(x))**4
    e = abs(d[0,0]/ex - 1)
    print("   h=%5.2f  FD=%.12e  rel err %.3e%s" % (h, d[0,0], e, "" if prev is None else "   ratio %.2f"%(prev/e)))
    prev=e

# ---- 1. build every candidate, counting corr_op calls ---------------------
CAND = [
    ("stock",      lambda: ds.Sigma2Builder(background=bg, apply_c0=False)),
    ("R0 ref",     lambda: R.R0Reference(background=bg, apply_c0=False)),
    ("R1 knots",   lambda: R.R1Knots(background=bg, apply_c0=False)),
    ("R2 density", lambda: R.R2Density(background=bg, apply_c0=False)),
    ("R2 FIXED",   lambda: K.R2Fixed(background=bg, apply_c0=False)),
    ("R7 akima",   lambda: K.R7Akima(background=bg, apply_c0=False)),
    ("R8 hermite", lambda: K.R8HermiteNodes(background=bg, apply_c0=False)),
    ("R5 EXACT",   lambda: K.R5Exact(background=bg, apply_c0=False)),
]

GAMS = [1.0, 17.3]
COSL = [1.0] + [float(np.cos(np.deg2rad(g/60.0))) for g in GAMS]

# 600 evaluation points, avoiding exact table nodes
rng = np.random.default_rng(11)
LO, HI = 406.0, 2328.0
xs = np.sort(rng.uniform(LO+0.05, HI-0.05, 600))
xs = xs[np.min(np.abs(xs[:,None]-L[None,:]), axis=1) > 1e-6]
print("\n   %d evaluation lambda in [%.0f, %.0f]" % (xs.size, LO, HI))

# reference density on every cos
ref = {}
for c in COSL:
    ref[c] = np.array([r5.matrix(c, float(x)) for x in xs])

print("\n=== 1. construction cost + per-call cost + accuracy vs R5Exact ===")
print(f"{'scheme':<11} {'corrop/cos':>10} {'build s':>8} {'us/call':>8} "
      f"{'max err %':>10} {'p99 %':>8} {'median %':>9}   {'argmax lambda':>13}")
results = {}
for name, mk in CAND:
    K.COUNT.n = 0
    b = mk()
    # count corr_op calls for ONE cos_gamma by monkey-patching the stock path too
    orig = co.C_fn
    cnt = {"n": 0}
    def counted(n1, t1, nn2, t2, _o=orig, _c=cnt):
        _c["n"] += 1
        return _o(n1, t1, nn2, t2)
    co.C_fn = counted
    ds._corr_op.C_fn = counted
    K._corr_op.C_fn = counted
    R._corr_op.C_fn = counted
    t0 = time.time()
    _ = b.matrix(1.0, 1000.0)          # triggers construction for cos=1
    build_s = time.time()-t0
    n_build = cnt["n"]
    co.C_fn = orig; ds._corr_op.C_fn = orig; K._corr_op.C_fn = orig; R._corr_op.C_fn = orig
    # per-call cost (post-construction)
    t0 = time.time()
    for x in xs[:200]: b.matrix(1.0, float(x))
    us = (time.time()-t0)/200*1e6
    # accuracy
    errs = []
    argmax_lam = float("nan"); mx = -1.0
    for c in COSL:
        got = np.array([b.matrix(c, float(x)) for x in xs])
        sc = np.maximum(np.abs(ref[c]).max(axis=(1,2)), 1e-300)
        e = np.abs(got-ref[c]).max(axis=(1,2))/sc
        errs.append(e)
        if e.max() > mx: mx = e.max(); argmax_lam = xs[int(np.argmax(e))]
    e = np.concatenate(errs)
    results[name] = (b, e)
    print(f"{name:<11} {n_build:>10d} {build_s:>8.2f} {us:>8.1f} "
          f"{100*e.max():>10.3f} {100*np.percentile(e,99):>8.3f} "
          f"{100*np.median(e):>9.2e}   {argmax_lam:>13.2f}")

# ---- 2. PSD on the 6x6 ----------------------------------------------------
print("\n=== 2. PSD violations of the 6x6 (min eig < 0), 600 lambda x 2 gammas ===")
print(f"{'scheme':<11} " + " ".join(f"{'g=%.1f nviol'%g:>14}" for g in GAMS)
      + f" {'worst min/max eig':>20}")
for name,_ in CAND:
    b = results[name][0]
    row = []; worst = np.inf
    for g in GAMS:
        c = float(np.cos(np.deg2rad(g/60.0)))
        nv = 0
        for x in xs:
            M = core._assemble_6x6(b, c, float(x))
            w = np.linalg.eigvalsh(0.5*(M+M.T))
            r = w[0]/max(abs(w[-1]), 1e-300)
            if w[0] < 0: nv += 1
            worst = min(worst, r)
        row.append(nv)
    print(f"{name:<11} " + " ".join(f"{v:>14d}" for v in row) + f" {worst:>20.3e}")
