"""Exact expectation of the placement-complete FK estimator (no Monte-Carlo).

Everything the Monte-Carlo of ``fk_complete_core`` averages is a polynomial in
the Gaussian base field ``z``, so its expectation can be evaluated in closed
form by Wick contraction on the *discrete* grid the simulation actually uses.
That gives a noise-free reference: it says what the estimator converges to, so
a disagreement with the analytic FK can be attributed to the estimator's
structure rather than to sampling.

Discrete conventions (identical to ``sachs_mc_core``)
----------------------------------------------------
Grid ``lam = linspace(lam_min, lam_source, N)``, step ``dlam``.  The linear
response telescopes, ``P[j->k] = prod_{i=j+1}^{k} resp[i] = (D_j / D_k)^2``, and
the driving field enters node ``k`` with weight ``dlam``:

    s1[k]  = resp[k] s1[k-1] + z[k] dlam            (linear arm)
    s2[k]  = resp[k] s2[k-1] + F(s1[k-1]) dlam      (one F vertex)
    kappa  = sum_k w_k s[k],   w_k = dlam (half at the two ends)

so ``kappa1_A = sum_j dlam Wd[j] z_A[j]`` with ``Wd[j] = sum_{k>=j} w_k (D_j/D_k)^2``.

The colored (AR(1)) base field has ``z[0] = L_0 xi``, ``z[k] = rho z[k-1] +
sqrt(1-rho^2) L_k xi_k``, hence the EXACT variance recursion
``A[k] = rho^2 A[k-1] + (1-rho^2) V[k]`` (not ``V[k]`` itself -- ``V`` varies
along the ray, so the chain is not stationary) and
``R[k,p] = rho^{|k-p|} A[min(k,p)]``.

The two placement classes
-------------------------
With the skew drive ``df[p] = 0.5 Q[p] : (z_p z_p - A[p])``,

    T1_AB = Cov(kappa2_A[z],  kappa_d1_B[df])   -- Q on the far observable leg
    T2_AB = Cov(kappa1_A[z],  kappa_d2_B)       -- Q on an F-vertex input leg

and ``FK = (T1 + T2) + (T1 + T2)^T``.  Wick gives

    T1_AB = sum_{j,p} dlam^2 Wd_j Wd_p  F_Abc Q[p]_Buv Phi[j,p]_bu Phi[j,p]_cv
    T2_AB = sum_{j,p} dlam Wd_j U[j,p]  Fs_Bbc Q[p]_cuv Psi[p]_Au Phi[j,p]_bv

with ``Phi[j,p] = <s1[j-1] z[p]^T>``, ``Psi[p]_Au = <kappa1_A z_u[p]>``,
``U[j,p] = dlam (D_p / D_{j-1})^2`` for ``p <= j-1`` and ``Fs = F + F^T`` on the
two input slots.  ``T1`` alone is what ``simulate_fk_vr`` measures.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import _bootstrap

_DS, _BG, _CORE = _bootstrap.wire()


# --- the paper F-vertex, and its linearisation's doubled form ---------------
F = _CORE._F.copy()                      # (3,3,3), F_abc s_b s_c
FS = F + F.transpose(0, 2, 1)            # symmetrised: DF(s).d = FS_abc s_b d_c


def _f6() -> NDArray[np.float64]:
    """Lift the per-ray (3,3,3) vertex to the two-ray (6,6,6) block-diagonal."""
    out = np.zeros((6, 6, 6))
    out[0:3, 0:3, 0:3] = F
    out[3:6, 3:6, 3:6] = F
    return out


def _fs6() -> NDArray[np.float64]:
    out = np.zeros((6, 6, 6))
    out[0:3, 0:3, 0:3] = FS
    out[3:6, 3:6, 3:6] = FS
    return out


F6 = _f6()
FS6 = _fs6()


@dataclass
class Grid:
    lam: NDArray[np.float64]
    dlam: float
    D: NDArray[np.float64]
    resp: NDArray[np.float64]
    Wd: NDArray[np.float64]              # (N,) discrete LOS window
    Hd: NDArray[np.float64]              # (N,) discrete F-leg response integral
    bg: object
    builder: object


def build_grid(n_lambda: int = 1000, lam_min: float = 406.0,
               lam_source: float | None = None, apply_anchor: bool = False) -> Grid:
    lam_f = float(lam_source if lam_source is not None else _BG.LAM_SOURCE_BASELINE)
    bg = _BG.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    lam = np.linspace(lam_min, lam_f, n_lambda)
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), float)
    resp = np.empty(n_lambda)
    resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    w = np.full(n_lambda, dlam)
    w[0] = w[-1] = 0.5 * dlam
    # Wd[j] = sum_{k>=j} w_k (D_j/D_k)^2
    tail = np.cumsum((w / D**2)[::-1])[::-1]
    Wd = D**2 * tail
    # Hd[m] = sum_{j>=m+1} dlam Wd[j] (D_m / D_{j-1})^4    (F-leg response)
    gj = dlam * Wd[1:] / D[:-1] ** 4                       # index j-1 = 0..N-2
    tail4 = np.zeros(n_lambda)
    tail4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D**4 * tail4
    builder = _DS.Sigma2Builder(background=bg, apply_c0=apply_anchor)
    return Grid(lam, dlam, D, resp, Wd, Hd, bg, builder)


def smeared_covariance(grid: Grid, A: NDArray, rho: float) -> NDArray[np.float64]:
    """``B[k] = sum_l dlam <z[k] z[l]^T>`` -- the covariance the drive actually sees.

    In the white-noise limit ``B -> Sigma2``; at finite ``sigma_lambda`` it is the
    AR(1)-smeared version, which is NOT proportional to ``Sigma2[k]``: the
    two-ray covariance changes direction along the ray, so smearing rotates it.
    """
    N = grid.lam.size
    idx = np.arange(N)
    B = np.empty_like(A)
    for k in range(N):
        wgt = rho ** np.abs(idx - k)
        B[k] = grid.dlam * np.einsum("l,lij->ij", wgt, A[np.minimum(idx, k)],
                                     optimize=True)
    return B


def node_stats(grid: Grid, cos_gamma: float, sigma_lambda: float,
               rho: float | None = None, calibrate: str = "smeared"):
    """Per-node ``V``, exact AR(1) variance ``A``, and the deformation ``Q``.

    ``calibrate``:

    * ``"nominal"`` -- the published choice: solve ``Q`` against ``V = Sigma2 /
      (2 sigma)`` with the target ``zeta / (2 sigma)^2``, i.e. assume the field
      is delta-correlated on the scale over which ``Sigma2`` varies.
    * ``"smeared"`` -- solve ``Q`` against ``B``, the covariance the AR(1) field
      actually presents to the drive.  Identical in the white-noise limit:
      ``solve_Q`` is homogeneous of degree ``-2`` in its matrix argument, since
      ``Z_abc = Q_amn M_mb M_nc + 2 perms`` is QUADRATIC in ``M``, so
      ``solve_Q(c M, Z) = c^-2 solve_Q(M, Z)`` and hence
      ``solve_Q(2 sigma V, zeta) == solve_Q(V, zeta / (2 sigma)^2)`` exactly.
      Since ``B -> Sigma2 = 2 sigma V`` as ``sigma_lambda -> 0``, the two
      calibrations share the white-noise limit; the smeared one simply does not
      assume ``Sigma2`` is constant over a correlation length.
    """
    N = grid.lam.size
    sig = float(sigma_lambda)
    V = np.empty((N, 6, 6))
    Z = np.empty((N, 6, 6, 6))
    for k, lk in enumerate(grid.lam):
        S2 = _CORE._assemble_6x6(grid.builder, float(cos_gamma), float(lk))
        V[k] = S2 / (2.0 * sig)
        Z[k] = _DS.zeta6(float(cos_gamma), float(lk))
    if rho is None:
        rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V)
    A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1.0 - rho * rho) * V[k]
    if calibrate == "nominal":
        M = V * (2.0 * sig)                       # = Sigma2, by construction
    elif calibrate == "smeared":
        M = smeared_covariance(grid, A, rho)
    else:
        raise ValueError(f"unknown calibrate={calibrate!r}")
    Q = np.array([_DS.solve_Q(M[k], Z[k]) for k in range(N)])
    return V, A, Q, rho


def calib_matrix(grid: Grid, V: NDArray, A: NDArray, rho: float,
                 sigma_lambda: float, calibrate: str) -> NDArray[np.float64]:
    """The matrix ``Q`` was calibrated against (for the local-limit reference)."""
    if calibrate == "nominal":
        return V * (2.0 * sigma_lambda)
    return smeared_covariance(grid, A, rho)


def zeta_injected(M: NDArray, Q: NDArray) -> NDArray:
    """The three-point cumulant density the deformation ACTUALLY injects.

    ``M`` is the matrix ``Q`` was calibrated against.  ``solve_Q`` fully
    symmetrises its target, so a tabulated ``zeta6`` that is not exactly
    symmetric (the vertex table is ~2.7% asymmetric under leg exchange, an
    interpolation artefact) is not reproduced exactly.  This is the honest
    target for the estimator; the residual against ``zeta6`` is the injection
    fidelity, reported separately.
    """
    out = np.empty_like(Q)
    for k in range(Q.shape[0]):
        out[k] = _DS.cum3_from_Q(Q[k], M[k])
    return out


def fk_reference(grid: Grid, zeta: NDArray[np.float64]) -> NDArray[np.float64]:
    """Local (white-noise) FK on this discrete grid: ``2 sum_m dlam Wd Hd F:zeta``."""
    kern = grid.dlam * grid.Wd * grid.Hd                      # (N,)
    M = np.einsum("m,Abc,mbcB->AB", kern, F6, zeta, optimize=True)
    return M + M.T


def _phi_block(grid: Grid, A: NDArray, rho: float, ps: NDArray[np.int64]):
    """``Phi[j, p] = <s1[j-1] z[p]^T>`` and ``R[k, p] = <z[k] z[p]^T>``.

    Vectorised over a block of source nodes ``ps`` (the recursion in ``k`` is
    inherently sequential, so the block axis is what keeps this out of a
    million-iteration Python loop).  Shapes ``(N, nb, 6, 6)``.
    """
    N = grid.lam.size
    nb = ps.size
    idx = np.arange(N)
    lag = np.abs(idx[:, None] - ps[None, :])                  # (N, nb)
    R = (rho ** lag)[:, :, None, None] * A[np.minimum(idx[:, None], ps[None, :])]
    phi = np.empty((N, nb, 6, 6))
    cur = np.zeros((nb, 6, 6))
    for k in range(N):
        cur = grid.resp[k] * cur + grid.dlam * R[k]
        phi[k] = cur
    Phi = np.empty_like(phi)
    Phi[0] = 0.0
    Phi[1:] = phi[:-1]
    return Phi, R


def expectation(grid: Grid, cos_gamma: float, sigma_lambda: float,
                V=None, A=None, Q=None, rho=None, block: int = 128,
                verbose: bool = False):
    """Exact ``T1``, ``T2`` (6x6 each, before the A<->B symmetrisation)."""
    if V is None:
        V, A, Q, rho = node_stats(grid, cos_gamma, sigma_lambda)
    N = grid.lam.size
    dlam, Wd, D = grid.dlam, grid.Wd, grid.D
    T1 = np.zeros((6, 6))
    T2 = np.zeros((6, 6))
    cj1 = dlam * Wd                                           # weight of the F node j
    jm1 = np.maximum(np.arange(N) - 1, 0)
    Dj1 = D[jm1]                                              # D_{j-1}
    jminus = np.arange(N) - 1
    for lo in range(0, N, block):
        ps = np.arange(lo, min(lo + block, N))
        Phi, R = _phi_block(grid, A, rho, ps)                 # (N, nb, 6, 6)
        # --- T1: Q on the far observable leg -------------------------------
        Ph = np.moveaxis(Phi, 1, 0).reshape(ps.size, N, 36)   # (nb, N, 36)
        G = np.einsum("pjb,j,pjc->pbc", Ph, cj1, Ph, optimize=True)
        G = G.reshape(ps.size, 6, 6, 6, 6)                    # [p, b, u, c, v]
        T1 += np.einsum("p,Abc,pBuv,pbucv->AB",
                        dlam * Wd[ps], F6, Q[ps], G, optimize=True)
        # --- T2: Q on an F-vertex input leg --------------------------------
        # U[j,p] = dlam (D_p / D_{j-1})^2 for p <= j-1
        U = np.where(jminus[:, None] >= ps[None, :],
                     dlam * (D[ps][None, :] / Dj1[:, None]) ** 2, 0.0)   # (N, nb)
        S = np.einsum("j,jp,jpbv->pbv", cj1, U, Phi, optimize=True)
        Psi = np.einsum("m,mpAu->pAu", dlam * Wd, R, optimize=True)
        T2 += np.einsum("Bbc,pcuv,pAu,pbv->AB",
                        FS6, Q[ps], Psi, S, optimize=True)
        if verbose:
            print(f"  [exact] p<{ps[-1] + 1}/{N}", flush=True)
    return T1, T2
