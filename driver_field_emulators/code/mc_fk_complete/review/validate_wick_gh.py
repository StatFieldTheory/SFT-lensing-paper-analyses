"""Wick-free verification of the O(Q) three-point formula and of C[k] = Cov[k].

Builds a small AR(1) chain (n nodes, d components), forms the exact joint
Gaussian covariance from the EXPLICIT linear map (no Wick), and evaluates
E[f_p,a f_q,b f_l,c] by tensor-product Gauss-Hermite quadrature, which is exact
for the degree-6 polynomial integrand.  <f f f> = alpha*t + beta*t^3 in the
vertex scale t, so the O(Q) coefficient is extracted exactly by

    alpha = (8 I(t) - I(2t)) / (6 t).

That alpha is compared with the analytic prediction

    Q[p]_{a m n} X[p,q]_{m b} X[p,l]_{n c} + (vertex on q) + (vertex on l).

Also verifies that using C[k] = V[k] (instead of Cov[k]) leaves a non-zero mean
<f> and a spurious tadpole, i.e. that C[k] = Cov[k] is forced.
"""
from __future__ import annotations
import itertools
import numpy as np
from numpy.polynomial.hermite import hermgauss


def build(n=3, d=2, sigma=2.0, dlam=0.9, seed=7):
    rng = np.random.default_rng(seed)
    V = np.empty((n, d, d))
    for k in range(n):
        M = rng.normal(size=(d, d))
        V[k] = (M @ M.T) / d + 0.5 * np.eye(d)
    rho = float(np.exp(-dlam / sigma))
    # explicit linear map z[k] = sum_i c[k,i] L[i] xi_i  (no Wick anywhere)
    c = np.zeros((n, n))
    c[0, 0] = 1.0
    for k in range(1, n):
        c[k, k] = np.sqrt(1.0 - rho**2)
    for k in range(1, n):
        c[k, :k] = rho * c[k - 1, :k]
    L = np.linalg.cholesky(V)
    Amap = np.zeros((n * d, n * d))
    for k in range(n):
        for i in range(k + 1):
            Amap[k*d:(k+1)*d, i*d:(i+1)*d] = c[k, i] * L[i]
    Sig = Amap @ Amap.T
    Cov = np.array([Sig[k*d:(k+1)*d, k*d:(k+1)*d] for k in range(n)])
    X = np.array([[Sig[k*d:(k+1)*d, p*d:(p+1)*d] for p in range(n)] for k in range(n)])
    Q0 = rng.normal(size=(n, d, d, d))
    Q0 = 0.5 * (Q0 + Q0.transpose(0, 1, 3, 2))
    return dict(n=n, d=d, V=V, Cov=Cov, X=X, Amap=Amap, Q0=Q0, rho=rho)


def gh_moments(B, t, C, nodes=6):
    """E[f_p,a f_q,b f_l,c] by exact tensor Gauss-Hermite. C is the subtraction."""
    n, d, Amap, Q0 = B["n"], B["d"], B["Amap"], B["Q0"]
    k = n * d
    x, w = hermgauss(nodes)
    grids = np.array(list(itertools.product(*([x] * k))))       # (nodes^k, k)
    wts = np.prod(np.array(list(itertools.product(*([w] * k)))), axis=1)
    wts = wts / wts.sum()
    z = (np.sqrt(2.0) * grids) @ Amap.T                          # (S, k)
    z = z.reshape(-1, n, d)
    Q = t * Q0
    f = (z + 0.5 * np.einsum("kamn,skm,skn->ska", Q, z, z, optimize=True)
         - 0.5 * np.einsum("kamn,kmn->ka", Q, C)[None])
    m3 = np.einsum("s,spa,sqb,slc->pqlabc", wts, f, f, f, optimize=True)
    m1 = np.einsum("s,ska->ka", wts, f, optimize=True)
    return m3, m1


def main():
    B = build()
    X, Cov, Q0 = B["X"], B["Cov"], B["Q0"]
    t = 0.05
    I1, mean1 = gh_moments(B, t, Cov)
    I2, _ = gh_moments(B, 2 * t, Cov)
    alpha = (8 * I1 - I2) / (6 * t)                # exact O(Q) coefficient / t-slope
    pred = (np.einsum("pamn,pqmb,plnc->pqlabc", Q0, X, X, optimize=True)
            + np.einsum("qbmn,qpma,qlnc->pqlabc", Q0, X, X, optimize=True)
            + np.einsum("lcmn,lpma,lqnb->pqlabc", Q0, X, X, optimize=True))
    rel = np.max(np.abs(alpha - pred)) / np.max(np.abs(pred))
    print(f"[GH] O(Q) 3-pt coefficient: max rel err vs Wick prediction = {rel:.3e}")
    print(f"[GH] |<f>| with C = Cov[k] : {np.max(np.abs(mean1)):.3e}"
          f"   (rel to |<z z>| scale {np.max(np.abs(Cov)):.3e})")
    _, meanV = gh_moments(B, t, B["V"])
    print(f"[GH] |<f>| with C = V[k]   : {np.max(np.abs(meanV)):.3e}  <-- non-zero: WRONG choice")
    # tadpole size if C = V
    I1v, _ = gh_moments(B, t, B["V"])
    I2v, _ = gh_moments(B, 2 * t, B["V"])
    alphav = (8 * I1v - I2v) / (6 * t)
    relv = np.max(np.abs(alphav - pred)) / np.max(np.abs(pred))
    print(f"[GH] with C = V[k]: max rel err vs Wick prediction = {relv:.3e}  <-- tadpole survives")


if __name__ == "__main__":
    main()
