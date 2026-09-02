"""PART 4 -- the decisive test of the "you fitted it" charge.

At a separation where the two-ray covariance is well conditioned the nominal
calibration should not be pathological, and BOTH calibrations must converge to
the SAME sigma_lambda -> 0 limit, equal to the local (white-noise) reference on
the same grid.

sigma_lambda is reduced at FIXED sigma/dlam (= 4.19), so the AR(1) chain is
always resolved by the same number of grid steps; that is the only way both
requirements for B -> Sigma2 (dlam/sigma -> 0 and sigma -> 0) are approached
together.
"""
import sys
import time
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                              # noqa: E402
from fold_ref import fk_table                             # noqa: E402

PAIRS = [(500, 16.0), (1000, 8.0), (1999, 4.0), (3997, 2.0)]
GAMMAS = [float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 \
    else [17.3, 30.0]

fold = fk_table()
fg = np.array(sorted(fold))
fv = np.array([fold[g] for g in fg])


def fold_at(gam):
    return float(np.exp(np.interp(np.log(gam), np.log(fg), np.log(fv))))


rows = []
for gam in GAMMAS:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print("=" * 108)
    print(f"gamma = {gam}'      paper fold FK(kk) = {fold_at(gam):.6e}"
          f"{'  (interpolated)' if round(gam,4) not in fold else ''}")
    print(f"{'calib':>8} {'N':>6} {'dlam':>7} {'sigma':>6} {'s/dl':>6} "
          f"{'|B-S2|/|S2| mid':>15} {'cond(S2) mid':>12} "
          f"{'exact FK':>13} {'local ref':>13} {'exact/ref':>10} "
          f"{'ref/fold':>9} {'exact/fold':>11} {'sec':>5}")
    for calib in ("nominal", "smeared"):
        for N, sig in PAIRS:
            t0 = time.time()
            grid = ex.build_grid(n_lambda=N)
            V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
            S2 = V * (2.0 * sig)
            B = ex.smeared_covariance(grid, A, rho)
            k = N // 2
            dev = float(np.abs(B[k] - S2[k]).max() / np.abs(S2[k]).max())
            cnd = float(np.linalg.cond(S2[k]))
            zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
            ref = ex.fk_reference(grid, zt)[0, 3]
            blk = 128 if N <= 2000 else 64
            T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho,
                                    block=blk)
            tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
            f = fold_at(gam)
            print(f"{calib:>8} {N:>6} {grid.dlam:>7.4f} {sig:>6.1f} "
                  f"{sig/grid.dlam:>6.2f} {dev:>15.4f} {cnd:>12.4e} "
                  f"{tot:>13.5e} {ref:>13.5e} {tot/ref:>10.4f} "
                  f"{ref/f:>9.4f} {tot/f:>11.4f} {time.time()-t0:>5.0f}",
                  flush=True)
            rows.append((gam, calib, N, sig, dev, cnd, tot, ref, f))
np.savez("_p4_decisive.npz", rows=np.array(rows, dtype=object))
