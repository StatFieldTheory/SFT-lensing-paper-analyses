"""Does the Monte-Carlo converge to the exact expectation of the same estimator?

Runs both on a small grid where the closed-form Wick contraction is cheap.  This
separates implementation error (MC vs exact, same estimator) from physics error
(exact vs analytic FK, different limits).
"""
import sys
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc

N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
sig = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
gam = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
nreal = int(sys.argv[4]) if len(sys.argv) > 4 else 40_000
nseed = int(sys.argv[5]) if len(sys.argv) > 5 else 4

cosg = float(np.cos(np.deg2rad(gam / 60.0)))
grid = ex.build_grid(n_lambda=N, lam_min=406.0)
nodes = ex.node_stats(grid, cosg, sig)
V, A, Q, rho = nodes
T1x, T2x = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
exact = ((T1x + T2x) + (T1x + T2x).T)[0, 3]
exact_vr = (T1x + T1x.T)[0, 3]
print(f"N={N} dlam={grid.dlam:.2f} sigma={sig} sig/dlam={sig/grid.dlam:.2f} "
      f"gamma={gam}' n_real={nreal} seeds={nseed}")
print(f"exact  T1+T2 (kk) = {exact:.6e}   T1 only = {exact_vr:.6e} "
      f"(share {exact_vr/exact:.4f})")

vals, vals_vr, vals_full, ses = [], [], [], []
for s in range(nseed):
    r = fc.simulate_fk_complete(gam, sigma_lambda=sig, n_lambda=N, n_real=nreal,
                                batch_size=min(4000, nreal), seed=1000 + s,
                                grid=grid, nodes=nodes)
    vals.append(r.kk); vals_vr.append(float(r.fk_vr_like[0, 0]))
    vals_full.append(float(r.fk_full[0, 0])); ses.append(float(r.fk_batch_se[0, 0]))
    print(f"  seed {s}: MC={r.kk:.6e} (batch SE {ses[-1]:.2e}) "
          f"ratio={r.kk/exact:.4f} | vr-like={vals_vr[-1]:.4e} "
          f"ratio_vr={vals_vr[-1]/exact_vr:.4f} | full={vals_full[-1]:.4e} "
          f"n_good={r.n_good}")
v = np.array(vals)
sem = v.std(ddof=1) / np.sqrt(len(v))
print(f"\nMC mean = {v.mean():.6e} +/- {sem:.2e} (seed SEM)   "
      f"ratio to exact = {v.mean()/exact:.4f} +/- {sem/abs(exact):.4f}")
print(f"seed scatter (std/mean) = {v.std(ddof=1)/abs(v.mean()):.3%}")
vv = np.array(vals_vr)
print(f"T1-only  ratio to exact T1-only = {vv.mean()/exact_vr:.4f}")
vf = np.array(vals_full)
print(f"nonlinear-arm flavour / perturbative = {vf.mean()/v.mean():.4f}")
