"""Where does the MC sit relative to its own exact expectation, and on what?

Same discretisation on both sides, so any gap is a sampling/estimator effect.
Scans n_lambda (grid) and batch_size (the in-sample centering's M).
"""
import sys, time
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc

gam = float(sys.argv[1]); sig = float(sys.argv[2])
seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 16
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
print(f"gamma={gam}' sigma={sig} seeds={seeds}")
print(f"{'N':>6} {'batch':>7} {'n_real':>8} {'exact':>12} {'MC':>12} "
      f"{'MC/exact':>9} {'sem':>8} {'T1 mc/ex':>9} {'T2 mc/ex':>9}")
for N, batch, nreal in [(250, 2000, 24000), (500, 2000, 24000),
                        (1000, 2000, 24000), (2000, 2000, 24000),
                        (1000, 500, 24000), (1000, 8000, 24000),
                        (1000, 2000, 96000)]:
    grid = ex.build_grid(n_lambda=N)
    nodes = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    V, A, Q, rho = nodes
    T1x, T2x = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    exact = ((T1x + T2x) + (T1x + T2x).T)[0, 3]
    vals, t1s, t2s = [], [], []
    for s in range(seeds):
        r = fc.simulate_fk_complete(gam, sigma_lambda=sig, n_lambda=N,
                                    n_real=nreal, batch_size=batch,
                                    seed=20260827 + 991 * s, grid=grid, nodes=nodes)
        vals.append(r.kk); t1s.append(r.t1[0, 3]); t2s.append(r.t2[0, 3])
    v = np.array(vals)
    print(f"{N:>6} {batch:>7} {nreal:>8} {exact:>12.5e} {v.mean():>12.5e} "
          f"{v.mean()/exact:>9.4f} {v.std(ddof=1)/np.sqrt(v.size)/abs(exact):>8.4f} "
          f"{np.mean(t1s)/T1x[0,3]:>9.4f} {np.mean(t2s)/T2x[0,3]:>9.4f}", flush=True)
