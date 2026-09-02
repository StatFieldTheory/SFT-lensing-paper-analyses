"""A7 -- close the loop: my own Monte-Carlo at an EXACT fold node, against the fold."""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc
import _bootstrap

d = np.load(_bootstrap.FOLD, allow_pickle=True)
m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
gs, vsf = [], np.asarray(d["value"], float)[m]
for x, y in zip(d["x"][m], d["y"][m]):
    x = np.asarray(x, float); y = np.asarray(y, float)
    gs.append(math.degrees(math.acos(float(np.clip(
        np.dot(x/np.linalg.norm(x), y/np.linalg.norm(y)), -1, 1)))) * 60)
gs = np.asarray(gs); o = np.argsort(gs); gs = gs[o]; vsf = vsf[o]
GAM, FOLD = float(gs[3]), float(vsf[3])
N, NREAL, NSEED = 1000, 24_000, 16
cosg = float(np.cos(np.deg2rad(GAM/60.0)))
grid = ex.build_grid(n_lambda=N)
print(f"gamma = {GAM:.5f}' (exact fold node)   fold = {FOLD:.6e}   N={N}  "
      f"n_real={NREAL} x {NSEED} seeds")
res = {}
for sig in (8.0, 4.0):
    nodes = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    V, A, Q, rho = nodes
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    exact = float(((T1+T2)+(T1+T2).T)[0, 3])
    vr = float((T1+T1.T)[0, 3])
    t0 = time.time(); vals = []
    for s in range(NSEED):
        r = fc.simulate_fk_complete(GAM, sig, n_lambda=N, n_real=NREAL,
                                    batch_size=2000, seed=771000+s,
                                    grid=grid, nodes=nodes)
        vals.append(r.kk)
    v = np.array(vals); mu = v.mean(); se = v.std(ddof=1)/math.sqrt(NSEED)
    res[sig] = (mu, se, exact)
    print(f"  sigma={sig:4.1f}  exact={exact:.6e}  MC={mu:.6e} +- {se:.2e}"
          f"   MC/exact={mu/exact:.4f}+-{se/exact:.4f}"
          f"   MC/fold={mu/FOLD:.4f}+-{se/FOLD:.4f}"
          f"   exact/fold={exact/FOLD:.4f}"
          f"   [T1-only/fold = {vr/FOLD:.4f}]   ({time.time()-t0:.0f}s)")
m8, s8, e8 = res[8.0]; m4, s4, e4 = res[4.0]
ex0 = 2*m4 - m8; se0 = math.sqrt((2*s4)**2 + s8**2)
print(f"\n  linear sigma->0 extrapolation (independent seeds, so errors add in quadrature):")
print(f"    MC    -> {ex0:.6e} +- {se0:.2e}   MC/fold    = {ex0/FOLD:.4f} +- {se0/FOLD:.4f}")
print(f"    exact -> {2*e4-e8:.6e}            exact/fold = {(2*e4-e8)/FOLD:.4f}")
