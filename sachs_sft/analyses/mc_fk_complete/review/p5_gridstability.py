"""PART 4b -- grid stability at FIXED sigma_lambda.

A calibration that is a well-posed prescription must give a grid-independent
answer once the grid resolves the correlation length.  This separates "the
regulator is not yet removed" (a limit statement) from "the prescription is
singular" (no limit at all).
"""
import sys
import time
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                             # noqa: E402
from fold_ref import fk_table                            # noqa: E402

fold = fk_table()
fg = np.array(sorted(fold)); fv = np.array([fold[g] for g in fg])
GAM = float(sys.argv[1]) if len(sys.argv) > 1 else 17.3
SIG = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
NS = [500, 707, 1000, 1414, 2000, 2828, 4000]
f = float(np.exp(np.interp(np.log(GAM), np.log(fg), np.log(fv))))
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
print(f"gamma = {GAM}'  sigma_lambda = {SIG} FIXED;  paper fold = {f:.6e}")
print(f"{'calib':>8} {'N':>6} {'dlam':>7} {'s/dl':>6} {'exact FK':>13} "
      f"{'exact/fold':>11} {'local ref':>13} {'sec':>5}")
for calib in ("nominal", "smeared"):
    for N in NS:
        t0 = time.time()
        grid = ex.build_grid(n_lambda=N)
        V, A, Q, rho = ex.node_stats(grid, cosg, SIG, calibrate=calib)
        zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
        ref = ex.fk_reference(grid, zt)[0, 3]
        T1, T2 = ex.expectation(grid, cosg, SIG, V=V, A=A, Q=Q, rho=rho,
                                block=128 if N <= 2000 else 64)
        tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
        print(f"{calib:>8} {N:>6} {grid.dlam:>7.4f} {SIG/grid.dlam:>6.2f} "
              f"{tot:>13.5e} {tot/f:>11.4f} {ref:>13.5e} {time.time()-t0:>5.0f}",
              flush=True)
