"""Alternative calibration matrices for the FK deformation tensor Q.

ADVERSARIAL: the author picked ONE calibration matrix (B[k] = sum_l dlam R[k,l])
and got sigma->0 -> 1.  If the sigma->0 limit is genuinely a property of the
construction, every OTHER matrix that also tends to Sigma2 as sigma->0 must give
the same limit.  If instead the limit moves by more than a percent between
equally defensible choices, the calibration was tuned.

Nothing here is imported by anything outside review/.
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import numpy as np
import fk_expect_exact as ex


def base_stats(grid, cos_gamma, sigma_lambda):
    """V, A, Z, rho -- everything EXCEPT the calibration choice."""
    N = grid.lam.size
    sig = float(sigma_lambda)
    V = np.empty((N, 6, 6)); Z = np.empty((N, 6, 6, 6))
    for k, lk in enumerate(grid.lam):
        S2 = ex._CORE._assemble_6x6(grid.builder, float(cos_gamma), float(lk))
        V[k] = S2 / (2.0 * sig)
        Z[k] = ex._DS.zeta6(float(cos_gamma), float(lk))
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1.0 - rho * rho) * V[k]
    return V, A, Z, rho


def _R_rows(grid, A, rho, k, idx):
    """R[k, l] for all l, shape (N,6,6):  rho^|k-l| A[min(k,l)]."""
    wgt = rho ** np.abs(idx - k)
    return wgt[:, None, None] * A[np.minimum(idx, k)]


def calib_matrices(grid, V, A, rho, sig, which):
    """Return M[k] (N,6,6) for calibration name `which`.

    All of these satisfy M -> Sigma2 in the joint (dlam/sigma -> 0,
    sigma/L_Sigma2 -> 0) limit EXCEPT the ones flagged BAD-LIMIT below.
    """
    N = grid.lam.size
    dlam = grid.dlam
    idx = np.arange(N)
    S2 = V * (2.0 * sig)

    if which == "nominal":                      # M = Sigma2 exactly (published)
        return S2

    if which == "smeared":                      # AUTHOR'S CHOICE
        return ex.smeared_covariance(grid, A, rho)

    if which == "smeared_V":                    # B built from nominal V, not A
        M = np.empty_like(V)
        for k in range(N):
            wgt = rho ** np.abs(idx - k)
            M[k] = dlam * np.einsum("l,lij->ij", wgt, V[np.minimum(idx, k)],
                                    optimize=True)
        return M

    if which == "smeared_Vk":                   # stationary-at-k: rho^|k-l| V[k]
        w = np.array([dlam * (rho ** np.abs(idx - k)).sum() for k in range(N)])
        return w[:, None, None] * V

    if which == "kernel":                       # kernel-weighted: Psi[k]/Wd[k]
        wd = dlam * grid.Wd
        M = np.empty_like(V)
        for k in range(N):
            M[k] = np.einsum("l,lij->ij", wd, _R_rows(grid, A, rho, k, idx),
                             optimize=True) / grid.Wd[k]
        return M

    if which == "kernel_fk":                    # weighted by the FK kernel dlam*Wd*Hd
        wk = dlam * grid.Wd * grid.Hd
        den = np.where(grid.Hd > 0, grid.Wd * grid.Hd, np.inf)
        M = np.empty_like(V)
        for k in range(N):
            M[k] = np.einsum("l,lij->ij", wk, _R_rows(grid, A, rho, k, idx),
                             optimize=True) / den[k]
        return M

    if which == "causal":                       # one-sided, NOT renormalised: -> Sigma2/2 (BAD LIMIT)
        M = np.empty_like(V)
        for k in range(N):
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = dlam * R[: k + 1].sum(axis=0)
        return M

    if which == "causal2":                      # one-sided x2  -> Sigma2  (defensible)
        M = np.empty_like(V)
        for k in range(N):
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = 2.0 * dlam * R[: k + 1].sum(axis=0)
        return M

    if which == "anticausal2":                  # the other one-sided x2
        M = np.empty_like(V)
        for k in range(N):
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = 2.0 * dlam * R[k:].sum(axis=0)
        return M

    if which == "causal2t":     # one-sided x2, TRAPEZOID endpoint -> O(dlam^2) accurate
        M = np.empty_like(V)
        for k in range(N):
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = 2.0 * dlam * (R[: k + 1].sum(axis=0) - 0.5 * R[k])
        return M

    if which == "anticausal2t":
        M = np.empty_like(V)
        for k in range(N):
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = 2.0 * dlam * (R[k:].sum(axis=0) - 0.5 * R[k])
        return M

    if which.startswith("trunc"):               # symmetric truncation at n*sigma
        nsig = float(which[5:])
        half = max(1, int(round(nsig * sig / dlam)))
        M = np.empty_like(V)
        for k in range(N):
            lo, hi = max(0, k - half), min(N, k + half + 1)
            R = _R_rows(grid, A, rho, k, idx)
            M[k] = dlam * R[lo:hi].sum(axis=0)
        return M

    if which.startswith("rtrunc"):              # truncated then RENORMALISED to full scalar weight
        nsig = float(which[6:])
        half = max(1, int(round(nsig * sig / dlam)))
        M = np.empty_like(V)
        for k in range(N):
            lo, hi = max(0, k - half), min(N, k + half + 1)
            R = _R_rows(grid, A, rho, k, idx)
            full = (rho ** np.abs(idx - k)).sum()
            part = (rho ** np.abs(idx[lo:hi] - k)).sum()
            M[k] = dlam * R[lo:hi].sum(axis=0) * (full / part)
        return M

    if which == "scalar":                       # Sigma2[k] * scalar smear weight (NO rotation)
        w = np.array([dlam * (rho ** np.abs(idx - k)).sum() for k in range(N)])
        return (w / (2.0 * sig))[:, None, None] * S2

    if which == "scalar_stat":                  # Sigma2[k] * (r/2)coth(r/2), interior formula
        r = dlam / sig
        fac = (r / 2.0) / np.tanh(r / 2.0)
        return fac * S2

    if which == "geo":                          # geometric mean of smeared and nominal (metric-aware)
        B = ex.smeared_covariance(grid, A, rho)
        M = np.empty_like(V)
        for k in range(N):
            d, U = np.linalg.eigh(S2[k])
            # symmetric sqrt of |S2|
            sq = U @ np.diag(np.sqrt(np.abs(d))) @ U.T
            isq = U @ np.diag(1.0 / np.sqrt(np.abs(d))) @ U.T
            X = isq @ B[k] @ isq
            ev, W = np.linalg.eigh((X + X.T) / 2)
            Xs = W @ np.diag(np.sqrt(np.clip(ev, 0, None))) @ W.T
            M[k] = sq @ Xs @ sq
        return M

    raise ValueError(which)


def run(gamma_arcmin, N, sig, names, block=None, verbose=False):
    cosg = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    grid = ex.build_grid(n_lambda=N)
    V, A, Z, rho = base_stats(grid, cosg, sig)
    if block is None:
        block = 128 if N <= 2000 else 64
    out = {}
    for nm in names:
        M = calib_matrices(grid, V, A, rho, sig, nm)
        Q = np.array([ex._DS.solve_Q(M[k], Z[k]) for k in range(N)])
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho,
                                block=block)
        tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
        vr = (T1 + T1.T)[0, 3]
        out[nm] = dict(tot=float(tot), t1=float(vr), Mmid=M[N // 2].copy())
        if verbose:
            print(f"    {nm:>12s}  tot={tot:.6e}", flush=True)
    return grid, out
