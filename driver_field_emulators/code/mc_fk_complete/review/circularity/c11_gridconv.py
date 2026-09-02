"""Is the baseline 1.00 a discretisation coincidence?  N-refinement at fixed
sigma/dlam ratio and at fixed sigma, at the fold's own gamma node."""
import _env, numpy as np
import fk_expect_exact as ex
from _fold import fold_nodes
from pathlib import Path
PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products")
gg, vv = fold_nodes(PROD / "table_permclosed_cut15360_permfix_xi.npz")
j = int(np.argmin(np.abs(gg - 1.0)))
fold = float(vv[j]); cosg = float(np.cos(np.deg2rad(gg[j] / 60.0)))
print(f"gamma node {gg[j]:.4f}'  fold {fold:.6e}")
print(f"{'N':>7}{'dlam':>8}{'s=8':>9}{'s=4':>9}{'r(s->0)':>10}")
for N in (500, 1000, 2000, 4000):
    grid = ex.build_grid(n_lambda=N)
    r = {}
    for sig in (8.0, 4.0):
        V, A, Q, rho = ex.node_stats(grid, cosg, sig)
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        r[sig] = float(((T1 + T2) + (T1 + T2).T)[0, 3]) / fold
    print(f"{N:>7}{grid.dlam:>8.2f}{r[8.0]:>9.4f}{r[4.0]:>9.4f}"
          f"{2*r[4.0]-r[8.0]:>10.4f}", flush=True)
