"""PART 2 -- the calibration is fixed by an observable with no F vertex in it.

<kappa1^3> is the third moment of the LINEAR observable.  Its local (white-noise)
value is sum_k dlam Wd_k^3 zeta[k].  The exact O(Q) value is computed here with
no approximation.  A calibration that is right must reproduce the local value as
the regulator is removed; the nominal one does not, the smeared one does.
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                    # noqa: E402
from kappa3_norm import (kappa1_third_moment, local_prediction,  # noqa: E402
                         smeared_B)

CHANNELS = [(0, 0, 0), (0, 0, 3), (0, 3, 3), (3, 3, 3), (0, 1, 1)]

# --- first: the O(N) recursion agrees with the production O(N^2) smearing ----
g = ex.build_grid(n_lambda=300)
cosg = float(np.cos(np.deg2rad(1.0 / 60.0)))
V, A, Q, rho = ex.node_stats(g, cosg, 8.0, calibrate="nominal")
B1 = ex.smeared_covariance(g, A, rho)
B2 = smeared_B(g, A, rho)
print(f"[check] O(N) smearing recursion vs production smeared_covariance: "
      f"rel {np.abs(B1 - B2).max() / np.abs(B1).max():.3e}")
print()

N = 1000
grid = ex.build_grid(n_lambda=N)
for gam in (1.0, 5.0, 17.3, 30.0):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    Zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    loc_true = local_prediction(grid, Zt)
    S2mid = None
    print("=" * 96)
    print(f"gamma = {gam}'    N = {N}, dlam = {grid.dlam:.4f}   "
          f"channel ratios  <k1^3>_exact / local prediction")
    print(f"{'calib':>8} {'sigma':>6} {'|B-S2|/|S2|':>12} " +
          " ".join(f"{str(c):>14}" for c in CHANNELS))
    for calib in ("nominal", "smeared"):
        for sig in (16.0, 8.0, 4.0, 2.0):
            V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
            S2 = V * (2.0 * sig)
            if S2mid is None:
                S2mid = np.linalg.cond(S2[N // 2])
            B = smeared_B(grid, A, rho)
            dev = np.abs(B - S2).max() / np.abs(S2).max()
            m3 = kappa1_third_moment(grid, A, Q, rho)
            row = " ".join(f"{m3[c] / loc_true[c]:>14.4f}" for c in CHANNELS)
            print(f"{calib:>8} {sig:>6.1f} {dev:>12.4f} {row}")
    print(f"  cond(Sigma2) at mid-ray = {S2mid:.3e}")
