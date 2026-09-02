"""Independent re-derivation of the identity that makes the sigma->0 limit provable.

Built ONLY from SPEC.md sections 1-2 (kernels re-derived here from scratch), then
compared against fk_expect_exact.expectation().  Nothing is imported from
gate1_placement_share.py.

From SPEC 1:  s[k] = resp[k] s[k-1] + F(s[k-1]) dlam + f[k] dlam,  s[-1]=0,
              kappa_a = sum_k w_k s_a[k].
Order-counting in F and in the drive's third cumulant gives, exactly,

  <kappa^(2)_A kappa^(1)_B>
     = sum_j dlam Wd_j F_Abc sum_{m,n<=j-1} U[j,m] U[j,n]
       sum_p dlam Wd_p  <f_b[m] f_c[n] f_B[p]>_c ,
  U[j,m] = dlam (D_m/D_{j-1})^2 ,   Wd_n = sum_{p>=n} w_p (D_n/D_p)^2 .

With f = z + 0.5 Q:(zz - A) and z a Gaussian AR(1) chain, the O(Q) connected
third cumulant of f is exactly

  <f_b[m] f_c[n] f_B[p]>_c = Q[m]_bij R[m,n]_ic R[m,p]_jB + (n) + (p) ,
  R[k,l] = rho^|k-l| A[min(k,l)] .

So the estimator's expectation IS the exact FK diagram of the colored field:
the FK KERNEL is sigma-independent and exact, and every bit of sigma-dependence
sits in the drive's third-cumulant DENSITY.  That reduces the sigma->0 question
to: does R -> Sigma2 delta strongly enough for the fold to converge, and at what
rate.
"""
from __future__ import annotations

import sys
import numpy as np

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/"
                   "driver_field_emulators/code/mc_fk_complete")
import fk_expect_exact as ex   # noqa: E402


def kernels_from_spec(N, lam_min=406.0):
    """Wd, Hd, U, resp built from SPEC 1 alone."""
    lam_f = float(ex._BG.LAM_SOURCE_BASELINE)
    bg = ex._BG.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    lam = np.linspace(lam_min, lam_f, N)
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), float)
    w = np.full(N, dlam); w[0] = w[-1] = 0.5 * dlam
    Wd = np.array([np.sum(w[n:] * (D[n] / D[n:]) ** 2) for n in range(N)])
    idx = np.arange(N)
    jm1 = np.maximum(idx - 1, 0)
    U = np.where((idx - 1)[:, None] >= idx[None, :],
                 dlam * (D[None, :] / D[jm1][:, None]) ** 2, 0.0)
    # Hd[m] = sum_j dlam Wd_j U[j,m]/dlam  (the F-leg response integral)
    Hd = (dlam * Wd)[:, None] * U / dlam
    Hd = np.einsum("jm->m", Hd) * 0 + np.array(
        [np.sum(dlam * Wd * U[:, m] / dlam * (D[m] / D[m]) ** 0) for m in range(N)])
    # Hd as used by the local reference: sum_j dlam Wd_j (D_m/D_{j-1})^4
    Hd4 = np.array([np.sum((dlam * Wd)[m + 1:] * (D[m] / D[jm1][m + 1:]) ** 4)
                    for m in range(N)])
    return lam, dlam, D, Wd, Hd4, U


def induced_cumulant(A, Q, rho, N):
    idx = np.arange(N)
    R = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None] \
        * A[np.minimum(idx[:, None], idx[None, :])]
    Xi = np.einsum("mbij,mnic,mpjB->mnpbcB", Q, R, R, optimize=True)
    Xi += np.einsum("ncij,nmib,npjB->mnpbcB", Q, R, R, optimize=True)
    Xi += np.einsum("pBij,pmib,pnjc->mnpbcB", Q, R, R, optimize=True)
    return Xi, R


def main():
    print(f"{'N':>4} {'gamma':>7} {'sigma':>7} {'s/dlam':>8} "
          f"{'estimator':>14} {'FK diagram':>14} {'max rel dev':>12} "
          f"{'local ref':>13} {'est/local':>10}")
    for N, gam, sig in [(20, 1.0, 60.0), (24, 1.0, 120.0), (24, 5.0, 60.0),
                        (28, 1.0, 30.0), (28, 1.0, 15.0), (28, 1.0, 8.0),
                        (32, 17.3, 40.0), (16, 0.5, 200.0)]:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        lam, dlam, D, Wd, Hd4, U = kernels_from_spec(N)
        grid = ex.build_grid(n_lambda=N)
        assert np.allclose(Wd, grid.Wd) and np.allclose(Hd4, grid.Hd), \
            "kernel re-derivation disagrees with fk_expect_exact"
        V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate="smeared")
        Xi, R = induced_cumulant(A, Q, rho, N)
        M = np.einsum("j,Abc,jm,jn,p,mnpbcB->AB", dlam * Wd, ex.F6, U, U,
                      dlam * Wd, Xi, optimize=True)
        diagram = M + M.T
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        est = (T1 + T2) + (T1 + T2).T
        dev = np.abs(est - diagram).max() / np.abs(diagram).max()
        B = ex.calib_matrix(grid, V, A, rho, sig, "smeared")
        loc = ex.fk_reference(grid, ex.zeta_injected(B, Q))[0, 3]
        print(f"{N:>4} {gam:>7.1f} {sig:>7.1f} {sig/dlam:>8.3f} "
              f"{est[0,3]:>14.7e} {diagram[0,3]:>14.7e} {dev:>12.2e} "
              f"{loc:>13.6e} {est[0,3]/loc:>10.6f}")


if __name__ == "__main__":
    main()
