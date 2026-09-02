"""Independent checks of cleanroom_exact.py.

(1) AR(1) covariance: build the joint covariance of the whole chain from the
    explicit linear map z[k] = sum_i c_{k,i} L[i] xi_i and compare with
    X[k,p] = rho^{|k-p|} Cov[min(k,p)].
(2) The O(Q) three-point formula (incl. the C = Cov tadpole subtraction) against
    a brute-force Monte-Carlo of the actual drive f in a small toy system.
(3) The fast (recursion-based) FK assembly against a literal O(N^3) triple sum
    built straight from R[j,p], Wd[j] and the explicit X tensor.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import cleanroom_exact as CR
from _env import ds


# ---------------------------------------------------------------------------
def check_ar1_covariance(seed=0, n=25, d=4, sigma=3.0, dlam=1.1):
    rng = np.random.default_rng(seed)
    V = np.empty((n, d, d))
    for k in range(n):
        M = rng.normal(size=(d, d))
        V[k] = M @ M.T / d + 0.3 * np.eye(d)
    rho = np.exp(-dlam / sigma)
    # explicit linear map
    c = np.zeros((n, n))
    c[0, 0] = 1.0
    for k in range(1, n):
        c[k, k] = np.sqrt(1.0 - rho**2)
    for k in range(1, n):
        c[k, :k] = rho * c[k - 1, :k]
    Xexp = np.einsum("ki,pi,ide->kpde", c, c, V, optimize=True)
    # formula
    Cov = np.empty_like(V)
    Cov[0] = V[0]
    for k in range(1, n):
        Cov[k] = rho**2 * Cov[k - 1] + (1 - rho**2) * V[k]
    idx = np.arange(n)
    Xfor = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None] * \
        Cov[np.minimum(idx[:, None], idx[None, :])]
    err = np.max(np.abs(Xexp - Xfor)) / np.max(np.abs(Xexp))
    print(f"[1] AR(1) covariance  max rel err = {err:.3e}")
    return err


# ---------------------------------------------------------------------------
def check_wick_3pt(seed=1, n=6, d=3, sigma=2.0, dlam=0.9, qscale=0.25,
                   nsamp=4_000_000, batch=250_000):
    """MC-check <f_p,a f_q,b f_l,c>|_O(Q) = sum of the 3 vertex placements."""
    rng = np.random.default_rng(seed)
    V = np.empty((n, d, d))
    for k in range(n):
        M = rng.normal(size=(d, d))
        V[k] = (M @ M.T) / d + 0.4 * np.eye(d)
    rho = float(np.exp(-dlam / sigma))
    Cov = np.empty_like(V)
    Cov[0] = V[0]
    for k in range(1, n):
        Cov[k] = rho**2 * Cov[k - 1] + (1 - rho**2) * V[k]
    idx = np.arange(n)
    X = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None] * \
        Cov[np.minimum(idx[:, None], idx[None, :])]
    Q = rng.normal(size=(n, d, d, d)) * qscale
    Q = 0.5 * (Q + Q.transpose(0, 1, 3, 2))
    L = np.linalg.cholesky(V)

    # analytic O(Q) prediction
    t3 = np.einsum("lcmn,lpma,lqnb->pqlabc", Q, X, X, optimize=True)
    pred = (np.einsum("pamn,pqmb,plnc->pqlabc", Q, X, X, optimize=True)
            + np.einsum("qbmn,qpma,qlnc->pqlabc", Q, X, X, optimize=True)
            + t3)

    acc = np.zeros((n, n, n, d, d, d))
    done = 0
    while done < nsamp:
        b = min(batch, nsamp - done)
        z = np.empty((b, n, d))
        z[:, 0] = rng.normal(size=(b, d)) @ L[0].T
        s = np.sqrt(1 - rho**2)
        for k in range(1, n):
            z[:, k] = rho * z[:, k - 1] + s * (rng.normal(size=(b, d)) @ L[k].T)
        f = (z + 0.5 * np.einsum("kamn,skm,skn->ska", Q, z, z, optimize=True)
             - 0.5 * np.einsum("kamn,kmn->ka", Q, Cov)[None])
        acc += np.einsum("spa,sqb,slc->pqlabc", f, f, f, optimize=True)
        done += b
    mc = acc / nsamp
    scale = np.max(np.abs(pred))
    err = np.max(np.abs(mc - pred)) / scale
    # MC noise floor estimate: compare a symmetric split
    print(f"[2] Wick O(Q) 3-pt   max |MC-pred|/max|pred| = {err:.3e}  "
          f"(nsamp={nsamp:,}, qscale={qscale})")
    return err


# ---------------------------------------------------------------------------
def check_fast_vs_brute(gamma_arcmin=1.0, N=40, sigma_lambda=8.0):
    cosg = float(np.cos(np.radians(gamma_arcmin / 60.0)))
    grid = CR.build_grid(N)
    S, Z = CR.driver_tables(cosg, grid["lam"])
    rho, V, Cov, rp = CR.ar1_stats(S, grid["dlam"], sigma_lambda)
    Q = CR.solve_Q_all(V, Z / (2.0 * sigma_lambda) ** 2)
    fast = CR.fk_exact(grid, Cov, rho, rp, Q)

    dlam, D2, Wd = grid["dlam"], grid["D2"], grid["Wd"]
    idx = np.arange(N)
    X = (rho ** np.abs(idx[:, None] - idx[None, :]))[:, :, None, None] * \
        Cov[np.minimum(idx[:, None], idx[None, :])]
    # R[j,p], j = 0..N-1 (j=0 row is zero: F(s[-1]) = 0)
    R = np.zeros((N, N))
    for j in range(1, N):
        R[j, : j] = dlam * D2[: j] / D2[j - 1]
    K = np.einsum("j,jp,jq->pq", Wd, R, R, optimize=True)

    T1A = T1B = T2A = T2B = 0.0
    for l in range(N):
        MA = np.einsum("pq,pmc,qnc->mn", K, X[l, :, :, 0:3], X[l, :, :, 0:3], optimize=True)
        MB = np.einsum("pq,pmc,qnc->mn", K, X[l, :, :, 3:6], X[l, :, :, 3:6], optimize=True)
        T1A -= Wd[l] * float(np.einsum("mn,mn->", Q[l, 3], MA))
        T1B -= Wd[l] * float(np.einsum("mn,mn->", Q[l, 0], MB))
    for p in range(N):
        Ep = np.einsum("q,qmc->mc", K[p], X[p], optimize=True)
        Yp = np.einsum("l,lnb->nb", Wd, X[p], optimize=True)
        T2A -= 2.0 * float(np.einsum("cmn,mc,n->", Q[p, 0:3], Ep[:, 0:3], Yp[:, 3]))
        T2B -= 2.0 * float(np.einsum("cmn,mc,n->", Q[p, 3:6], Ep[:, 3:6], Yp[:, 0]))
    brute = dict(T1_A=T1A, T1_B=T1B, T2_A=T2A, T2_B=T2B,
                 T1=T1A + T1B, T2=T2A + T2B, total=T1A + T1B + T2A + T2B)
    print(f"[3] fast vs brute (N={N}, sigma={sigma_lambda}):")
    for k in ("T1_A", "T1_B", "T2_A", "T2_B", "total"):
        r = fast[k] / brute[k] if brute[k] != 0 else np.nan
        print(f"     {k:5s} fast={fast[k]:+.10e} brute={brute[k]:+.10e} ratio={r:.12f}")
    return fast, brute


if __name__ == "__main__":
    # check_wick_3pt is a noisy MC; validate_wick_gh.py supersedes it with an
    # exact Gauss-Hermite evaluation.  Run it explicitly if you want it.
    check_ar1_covariance()
    check_fast_vs_brute()
