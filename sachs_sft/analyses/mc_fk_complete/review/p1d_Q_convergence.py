"""PART 1d -- Q_smeared -> Q_nominal as the regulator is removed.

Along the fixed sigma/dlam sequence (both dlam/sigma -> 0 and sigma -> 0), B ->
Sigma2, so the two calibrations must produce the same Q.  Reported as a
node-wise relative distance: the median node (where Q is well defined) and the
worst node (where Sigma2 nearly vanishes and Q_nominal is singular, so no
agreement is possible or meaningful).
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                # noqa: E402
from kappa3_norm import smeared_B                           # noqa: E402

_DS = ex._DS
PAIRS = [(500, 16.0), (1000, 8.0), (1999, 4.0), (3997, 2.0)]
for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    print(f"gamma = {gam}'")
    print(f"{'N':>6} {'sigma':>6} {'|B-S2|/|S2| med':>16} "
          f"{'|Qs-Qn|/|Qn| med':>17} {'p90':>10} {'worst node':>11} {'lam':>9}")
    for N, sig in PAIRS:
        grid = ex.build_grid(n_lambda=N)
        V, A, Qn, rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
        S2 = V * (2.0 * sig)
        B = smeared_B(grid, A, rho)
        Qs = np.array([_DS.solve_Q(B[k], _DS.zeta6(cosg, float(grid.lam[k])))
                       for k in range(N)])
        num = np.abs(Qs - Qn).max(axis=(1, 2, 3))
        den = np.maximum(np.abs(Qn).max(axis=(1, 2, 3)), 1e-300)
        d = num / den
        dv = np.abs(B - S2).max(axis=(1, 2)) / np.abs(S2).max(axis=(1, 2))
        w = int(np.argmax(d))
        print(f"{N:>6} {sig:>6.1f} {np.median(dv):>16.4f} "
              f"{np.median(d):>17.4f} {np.percentile(d, 90):>10.4f} "
              f"{d[w]:>11.3e} {grid.lam[w]:>9.1f}", flush=True)
    print()
