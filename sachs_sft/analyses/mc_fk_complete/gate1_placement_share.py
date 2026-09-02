"""Gate 1 -- is the new estimator placement-complete?

The Monte-Carlo does not inject the three-point cumulant directly; it deforms a
Gaussian field, ``f = z + 0.5 Q (z z - A)``.  Wick-expanding that at O(Q) gives a
three-point function with THREE terms, one for each leg that carries the ``Q``:

    <f_b[m] f_c[n] f_B[p]> = Q[m]_bij R[m,n]_ic R[m,p]_jB
                           + Q[n]_cij R[n,m]_ib R[n,p]_jB
                           + Q[p]_Bij R[p,m]_ib R[p,n]_jc

This script builds that object exactly on a small grid, contracts it through the
discrete FK kernels to get the truth, and compares:

  * ``T1 + T1^T``      -- the ``simulate_fk_vr`` structure (Q on the observable
                          leg only).  Expected to be a fraction of the truth.
  * ``(T1+T2) + ...^T`` -- the new estimator.  Must equal the truth EXACTLY.

Nothing here is statistical: both sides are closed-form contractions, so the
verdict is a machine-precision identity, not a convergence statement.
"""
from __future__ import annotations

import numpy as np

import fk_expect_exact as ex


def discrete_truth(grid, A, Q, rho, F6):
    """Exact O(Q) FK on the discrete grid, from the induced three-point function."""
    N = grid.lam.size
    dlam, Wd, D = grid.dlam, grid.Wd, grid.D
    idx = np.arange(N)
    R = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None] \
        * A[np.minimum(idx[:, None], idx[None, :])]                 # (N,N,6,6)
    # Xi[m,n,p]_{bcB}
    Xi = np.einsum("mbij,mnic,mpjB->mnpbcB", Q, R, R, optimize=True)
    Xi = Xi + np.einsum("ncij,nmib,npjB->mnpbcB", Q, R, R, optimize=True)
    Xi = Xi + np.einsum("pBij,pmib,pnjc->mnpbcB", Q, R, R, optimize=True)
    # U[j,m] = dlam (D_m / D_{j-1})^2 for m <= j-1
    jm1 = np.maximum(idx - 1, 0)
    U = np.where((idx - 1)[:, None] >= idx[None, :],
                 dlam * (D[None, :] / D[jm1][:, None]) ** 2, 0.0)    # (N,N)
    # <kappa2_A kappa1_B> = sum_j dlam Wd_j F_Abc sum_mn U U sum_p dlam Wd_p Xi
    Mfull = np.einsum("j,Abc,jm,jn,p,mnpbcB->AB",
                      dlam * Wd, F6, U, U, dlam * Wd, Xi, optimize=True)
    return Mfull + Mfull.T, Xi, R, U


def main() -> int:
    N = 24
    sig = 60.0                     # coarse grid: keep sigma/dlam sane for the toy
    gam = 1.0
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=N, lam_min=406.0)
    print(f"toy grid: N={N}, dlam={grid.dlam:.1f} Mpc, sigma_lambda={sig}, "
          f"sigma/dlam={sig / grid.dlam:.2f}")
    V, A, Q, rho = ex.node_stats(grid, cosg, sig)
    truth, Xi, R, U = discrete_truth(grid, A, Q, rho, ex.F6)

    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    vr = T1 + T1.T
    new = (T1 + T2) + (T1 + T2).T

    def kk(M):
        return float(M[0, 3])          # cross-ray kappa-kappa

    print()
    print(f"{'quantity':<34} {'kappa-kappa':>14} {'/ truth':>10}")
    print(f"{'exact discrete truth':<34} {kk(truth):>14.6e} {1.0:>10.4f}")
    print(f"{'T1 + T1^T  (simulate_fk_vr)':<34} {kk(vr):>14.6e} "
          f"{kk(vr) / kk(truth):>10.4f}")
    print(f"{'T2 + T2^T  (F-leg placements)':<34} {kk(T2 + T2.T):>14.6e} "
          f"{kk(T2 + T2.T) / kk(truth):>10.4f}")
    print(f"{'(T1+T2) + transpose  (NEW)':<34} {kk(new):>14.6e} "
          f"{kk(new) / kk(truth):>10.4f}")

    err = np.abs(new - truth).max() / np.abs(truth).max()
    print(f"\nfull (6x6) max relative deviation, new vs truth: {err:.3e}")
    print("VERDICT:", "PLACEMENT-COMPLETE" if err < 1e-10 else "INCOMPLETE")
    return 0 if err < 1e-10 else 1


if __name__ == "__main__":
    raise SystemExit(main())
