"""Is the reference's gamma-growing excess over the fold the vertex table's leg
asymmetry, and does solve_Q's symmetrisation absorb it?

zeta6 assembles the (6,6,6) two-ray cumulant by evaluating the tabulated vertex
at each ORDERED ray triple.  The physical vertex is symmetric under leg
exchange, so any dependence on which slot carries n2 is a table artefact.
solve_Q symmetrises its target, so the Monte-Carlo injects the symmetrised
object; the reference must therefore be built from the same.
"""
import numpy as np, fk_expect_exact as ex

N = 1000
grid = ex.build_grid(n_lambda=N)
print(f"{'gamma':>7} {'ref single':>13} {'ref sym':>13} {'sym/single':>10} "
      f"{'ref_inj':>13} {'inj/sym':>8} {'legspread':>10}")
for gam in (0.5, 1.0, 2.0, 5.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    zs = (zt + zt.transpose(0, 1, 3, 2) + zt.transpose(0, 2, 1, 3)
          + zt.transpose(0, 2, 3, 1) + zt.transpose(0, 3, 1, 2)
          + zt.transpose(0, 3, 2, 1)) / 6.0
    r_single = ex.fk_reference(grid, zt)[0, 3]
    r_sym = ex.fk_reference(grid, zs)[0, 3]
    V, A, Q, rho = ex.node_stats(grid, cosg, 8.0, calibrate="smeared")
    M = ex.calib_matrix(grid, V, A, rho, 8.0, "smeared")
    r_inj = ex.fk_reference(grid, ex.zeta_injected(M, Q))[0, 3]
    # leg spread: how much the contraction depends on which slot carries ray 1
    mid = N // 2
    z = zt[mid]
    legs = [-sum(z[c, c, 3] for c in range(3)),
            -sum(z[c, 3, c] for c in range(3)),
            -sum(z[3, c, c] for c in range(3))]
    spread = (max(legs) - min(legs)) / abs(np.mean(legs))
    print(f"{gam:>7.2f} {r_single:>13.5e} {r_sym:>13.5e} "
          f"{r_sym/r_single:>10.5f} {r_inj:>13.5e} {r_inj/r_sym:>8.5f} "
          f"{spread:>10.4f}")
