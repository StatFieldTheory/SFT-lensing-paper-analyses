"""SUPERSEDED by probe_norm2.py (which sweeps both calibrations).

Kept because it is the probe that located defect 3.2.  Uses the pre-fix
node_stats signature, so run it against the "nominal" calibration only.

Normalisation test that does NOT involve the F vertex.

The third moment of the LINEAR observable, <k1_a k1_b k1_c> = sum_klm dlam^3
Wd Wd Wd Xi[k,l,m], must equal its local (white-noise) prediction
sum_k dlam Wd_k^3 zeta[k].  If the injected cumulant's smearing normalisation is
right this ratio is 1 regardless of the FK kernel.
"""
import numpy as np, fk_expect_exact as ex

gam, N = 1.0, 1000
cosg = float(np.cos(np.deg2rad(gam/60.0)))
grid = ex.build_grid(n_lambda=N)
idx = np.arange(N)
print(f"{'sigma':>7} {'<k^3> smeared':>15} {'<k^3> local':>15} {'ratio':>8} "
      f"{'Psi/WdS2':>10} {'B/S2':>8}")
for sig in (1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="nominal")
    zi = ex.zeta_injected(ex.calib_matrix(grid, V, A, rho, sig, "nominal"), Q)          # what the deformation injects, local
    R = (rho ** np.abs(idx[:, None]-idx[None, :]))[:, :, None, None] \
        * A[np.minimum(idx[:, None], idx[None, :])]
    wd = grid.dlam * grid.Wd
    Psi = np.einsum("l,klib->kib", wd, R, optimize=True)      # (N,6,6)
    B = np.einsum("l,klib->kib", np.full(N, grid.dlam), R, optimize=True)
    t1 = np.einsum("k,kaij,kib,kjc->abc", wd, Q, Psi, Psi, optimize=True)
    smeared = t1 + t1.transpose(1, 0, 2) + t1.transpose(2, 1, 0)
    local = np.einsum("k,kabc->abc", grid.dlam*grid.Wd**3, zi, optimize=True)
    a, b, c = 0, 0, 3
    S2mid = A[N//2]*(2*sig)
    print(f"{sig:>7.1f} {smeared[a,b,c]:>15.5e} {local[a,b,c]:>15.5e} "
          f"{smeared[a,b,c]/local[a,b,c]:>8.4f} "
          f"{np.abs(Psi[N//2]).max()/(grid.Wd[N//2]*np.abs(S2mid).max()):>10.4f} "
          f"{np.abs(B[N//2]).max()/np.abs(S2mid).max():>8.4f}")
