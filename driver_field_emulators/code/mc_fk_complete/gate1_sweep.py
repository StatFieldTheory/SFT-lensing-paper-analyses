"""Gate 1 across configurations: the identity must be structural."""
import numpy as np
import fk_expect_exact as ex
from gate1_placement_share import discrete_truth

print(f"{'N':>4} {'gamma':>7} {'sigma':>7} {'sig/dlam':>9} {'vr share':>9} "
      f"{'new/truth':>10} {'max rel dev':>12}")
for N, gam, sig in [(24, 1.0, 60.0), (24, 0.5, 60.0), (32, 5.0, 40.0),
                    (16, 1.0, 120.0), (40, 17.3, 30.0), (24, 1.0, 200.0)]:
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=N, lam_min=406.0)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig)
    truth, *_ = discrete_truth(grid, A, Q, rho, ex.F6)
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    new = (T1 + T2) + (T1 + T2).T
    vr = T1 + T1.T
    dev = np.abs(new - truth).max() / np.abs(truth).max()
    print(f"{N:>4} {gam:>7.1f} {sig:>7.1f} {sig/grid.dlam:>9.2f} "
          f"{vr[0,3]/truth[0,3]:>9.4f} {new[0,3]/truth[0,3]:>10.6f} {dev:>12.2e}")
