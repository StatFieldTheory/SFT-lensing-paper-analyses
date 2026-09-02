"""Exact (noise-free) sigma_lambda dependence of the estimator's expectation."""
import sys, time
import numpy as np
import fk_expect_exact as ex

gam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
pairs = sys.argv[2] if len(sys.argv) > 2 else "1000:8"
calib = sys.argv[3] if len(sys.argv) > 3 else "smeared"
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
print(f"gamma = {gam}'   calibration = {calib}")
print(f"{'N':>6} {'dlam':>7} {'sigma':>7} {'s/dl':>6} {'ref_inj':>12} "
      f"{'T1+T2':>12} {'ratio':>8} {'share':>7} {'vr/ref':>8} {'ref_true':>12} {'sec':>5}")
for item in pairs.split(","):
    Ns, sigs = item.split(":")
    N, sig = int(Ns), float(sigs)
    t0 = time.time()
    grid = ex.build_grid(n_lambda=N)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
    M = ex.calib_matrix(grid, V, A, rho, sig, calib)
    zi = ex.zeta_injected(M, Q)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    ref_inj = ex.fk_reference(grid, zi)[0, 3]
    ref_true = ex.fk_reference(grid, zt)[0, 3]
    blk = 128 if N <= 2000 else 64
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho, block=blk)
    tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
    vr = (T1 + T1.T)[0, 3]
    print(f"{N:>6} {grid.dlam:>7.3f} {sig:>7.1f} {sig/grid.dlam:>6.2f} "
          f"{ref_inj:>12.5e} {tot:>12.5e} {tot/ref_inj:>8.4f} "
          f"{vr/tot:>7.4f} {vr/ref_inj:>8.4f} {ref_true:>12.5e} "
          f"{time.time()-t0:>5.0f}", flush=True)
