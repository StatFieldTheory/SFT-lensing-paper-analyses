"""PART 2b -- separate the calibration error from the vertex table's leg asymmetry.

solve_Q fully symmetrises its target, so the antisymmetric part of the tabulated
zeta6 under leg exchange cannot be injected by ANY Q.  The honest test of a
calibration is against the symmetrised target; the residual against the raw
table is an injection-fidelity floor shared by both calibrations.
"""
import sys
import numpy as np

sys.path.insert(0, "..")
import fk_expect_exact as ex                                        # noqa: E402
from kappa3_norm import kappa1_third_moment, local_prediction, smeared_B  # noqa: E402

_DS = ex._DS
N, grid = 1000, None
grid = ex.build_grid(n_lambda=1000)
CH = [(0, 0, 0), (0, 0, 3), (0, 3, 3), (0, 1, 1)]


def sym3(T):
    return sum(T.transpose(0, *p) if False else T.transpose(p) for p in
               [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]) / 6.0


for gam in (1.0, 17.3, 30.0):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    Zt = np.array([_DS.zeta6(cosg, float(l)) for l in grid.lam])
    Zs = np.array([sym3(Zt[k]) for k in range(N)])
    asym = np.abs(Zt - Zs).max() / np.abs(Zt).max()
    Ltrue = local_prediction(grid, Zt)
    Lsym = local_prediction(grid, Zs)
    print(f"gamma = {gam}'   max leg-asymmetry of zeta6 = {asym*100:.2f}% ; "
          "channel ratios sym/raw = " +
          " ".join(f"{c}:{Lsym[c]/Ltrue[c]:.4f}" for c in CH))
    print(f"{'calib':>8} {'sigma':>6} " +
          " ".join(f"{str(c)+' /sym':>16}" for c in CH))
    for calib in ("nominal", "smeared"):
        for sig in (8.0, 2.0):
            V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
            m3 = kappa1_third_moment(grid, A, Q, rho)
            print(f"{calib:>8} {sig:>6.1f} " +
                  " ".join(f"{m3[c]/Lsym[c]:>16.4f}" for c in CH))
    print()
