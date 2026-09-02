"""Normalisation of the injected three-point cumulant, both calibrations."""
import numpy as np, fk_expect_exact as ex
gam, N = 1.0, 1000
cosg = float(np.cos(np.deg2rad(gam/60.0)))
grid = ex.build_grid(n_lambda=N); idx = np.arange(N)
print(f"{'calib':>9} {'sigma':>6} {'<k^3> smear':>14} {'<k^3> local':>14} {'ratio':>8}")
for calib in ("nominal", "smeared"):
    for sig in (2.0, 4.0, 8.0, 16.0, 32.0):
        V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
        M = ex.calib_matrix(grid, V, A, rho, sig, calib)
        zi = ex.zeta_injected(M, Q)
        R = (rho ** np.abs(idx[:, None]-idx[None, :]))[:, :, None, None] \
            * A[np.minimum(idx[:, None], idx[None, :])]
        wd = grid.dlam * grid.Wd
        Psi = np.einsum("l,klib->kib", wd, R, optimize=True)
        t1 = np.einsum("k,kaij,kib,kjc->abc", wd, Q, Psi, Psi, optimize=True)
        sm = t1 + t1.transpose(1,0,2) + t1.transpose(2,1,0)
        lo = np.einsum("k,kabc->abc", grid.dlam*grid.Wd**3, zi, optimize=True)
        print(f"{calib:>9} {sig:>6.1f} {sm[0,0,3]:>14.5e} {lo[0,0,3]:>14.5e} "
              f"{sm[0,0,3]/lo[0,0,3]:>8.4f}")
