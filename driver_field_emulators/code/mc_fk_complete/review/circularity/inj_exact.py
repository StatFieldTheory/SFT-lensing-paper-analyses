"""Exact (Wick) expectation of the FK estimator, with DELIBERATE defects.

Mirrors ``fk_expect_exact.expectation`` term for term, but every structural
ingredient that is NOT read from the shared vertex table is switchable:

    response exponent, the causal window on the F leg, the lambda measure,
    the trapezoid end weights, the F-vertex index contraction, the A<->B
    symmetrisation, and the amplitude of the injected cumulant.

Nothing in the production folder is modified; this file only imports it.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import fk_expect_exact as ex

_DS, _BG, _CORE = ex._DS, ex._BG, ex._CORE


def build_grid(n_lambda=1000, lam_min=406.0, lam_source=None,
               resp_pow=2.0, end_weights=True):
    lam_f = float(lam_source if lam_source is not None else _BG.LAM_SOURCE_BASELINE)
    bg = _BG.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    lam = np.linspace(lam_min, lam_f, n_lambda)
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), float)
    resp = np.empty(n_lambda); resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** resp_pow
    w = np.full(n_lambda, dlam)
    if end_weights:
        w[0] = w[-1] = 0.5 * dlam
    # Wd[j] = sum_{k>=j} w_k prod resp  -> (D_j/D_k)^resp_pow
    Dp = D ** resp_pow
    tail = np.cumsum((w / Dp)[::-1])[::-1]
    Wd = Dp * tail
    gj = dlam * Wd[1:] / (D[:-1] ** (2 * resp_pow))
    tail4 = np.zeros(n_lambda); tail4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D ** (2 * resp_pow) * tail4
    builder = _DS.Sigma2Builder(background=bg, apply_c0=False)
    return ex.Grid(lam, dlam, D, resp, Wd, Hd, bg, builder)


def _phi_block(grid, A, rho, ps, causal_shift=0):
    N = grid.lam.size; nb = ps.size
    idx = np.arange(N)
    lag = np.abs(idx[:, None] - ps[None, :])
    R = (rho ** lag)[:, :, None, None] * A[np.minimum(idx[:, None], ps[None, :])]
    phi = np.empty((N, nb, 6, 6))
    cur = np.zeros((nb, 6, 6))
    for k in range(N):
        cur = grid.resp[k] * cur + grid.dlam * R[k]
        phi[k] = cur
    Phi = np.empty_like(phi)
    if causal_shift == 0:            # Phi[j] = <s1[j-1] z^T>   (correct)
        Phi[0] = 0.0; Phi[1:] = phi[:-1]
    else:                            # Phi[j] = <s1[j]   z^T>   (acausal)
        Phi[:] = phi
    return Phi, R


def expectation(grid, cos_gamma, sigma_lambda, V, A, Q, rho, block=128,
                defect="none"):
    """T1, T2 (6x6) with one structural defect injected."""
    N = grid.lam.size
    dlam, Wd, D = grid.dlam, grid.Wd, grid.D
    F6, FS6 = ex.F6.copy(), ex.FS6.copy()
    Qs = Q
    u_pow = 2.0
    causal_shift = 0
    u_offset = 1          # U nonzero for p <= j-u_offset
    t2_dlam = 1.0

    if defect == "none":
        pass
    elif defect == "acausal_U":            # p <= j   instead of p <= j-1
        u_offset = 0
    elif defect == "acausal_Phi":          # F sees the NEW state
        causal_shift = 1
    elif defect == "acausal_both":
        u_offset = 0; causal_shift = 1
    elif defect == "U_pow1":               # one propagator on the F leg
        u_pow = 1.0
    elif defect == "F_no_leg_sym":         # FS -> F in T2
        FS6 = ex.F6.copy()
    elif defect == "F_coef_half":          # F_101 = F_202 : -2 -> -1
        F6 = ex.F6.copy(); FS6 = ex.FS6.copy()
        for base in (0, 3):
            F6[base + 1, base + 0, base + 1] = -1.0
            F6[base + 2, base + 0, base + 2] = -1.0
        FS6 = F6 + F6.transpose(0, 2, 1)
    elif defect == "F_legswap":            # transpose the two F input legs
        F6 = ex.F6.transpose(0, 2, 1).copy()
        FS6 = F6 + F6.transpose(0, 2, 1)
    elif defect == "zeta_x1.3":            # amplitude, estimator side only
        Qs = 1.3 * Q
    elif defect == "drop_dlam_T2":         # one lambda measure missing in T2
        t2_dlam = 1.0 / dlam
    elif defect == "no_transpose":
        pass                                # handled by the caller
    else:
        raise ValueError(defect)

    T1 = np.zeros((6, 6)); T2 = np.zeros((6, 6))
    cj1 = dlam * Wd
    jm1 = np.maximum(np.arange(N) - 1, 0)
    Dj1 = D[jm1]
    jminus = np.arange(N) - u_offset
    for lo in range(0, N, block):
        ps = np.arange(lo, min(lo + block, N))
        Phi, R = _phi_block(grid, A, rho, ps, causal_shift)
        Ph = np.moveaxis(Phi, 1, 0).reshape(ps.size, N, 36)
        G = np.einsum("pjb,j,pjc->pbc", Ph, cj1, Ph, optimize=True)
        G = G.reshape(ps.size, 6, 6, 6, 6)
        T1 += np.einsum("p,Abc,pBuv,pbucv->AB", dlam * Wd[ps], F6, Qs[ps], G,
                        optimize=True)
        U = np.where(jminus[:, None] >= ps[None, :],
                     dlam * (D[ps][None, :] / Dj1[:, None]) ** u_pow, 0.0)
        S = np.einsum("j,jp,jpbv->pbv", cj1, U, Phi, optimize=True)
        Psi = np.einsum("m,mpAu->pAu", t2_dlam * dlam * Wd, R, optimize=True)
        T2 += np.einsum("Bbc,pcuv,pAu,pbv->AB", FS6, Qs[ps], Psi, S,
                        optimize=True)
    return T1, T2


def fk_kk(T1, T2, defect="none"):
    tot = T1 + T2
    if defect == "no_transpose":
        return float(tot[0, 3])
    return float((tot + tot.T)[0, 3])


def build_grid_warped(n_lambda=1000, lam_min=406.0, lam_source=None, eps=0.0,
                      mode="linear"):
    """Grid whose background D(lambda) is deliberately warped by (1 + eps u).

    ``mode='linear'``  u = (lam-lam_min)/(lam_f-lam_min)   (a tilt)
    ``mode='const'``   u = 1                               (pure rescaling)
    D is a SHARED input: the MC's background.py is a thin wrapper over exactly
    the D_callable.py that the sft-wick fold uses as its R_time_module.  So a
    common error in D cancels in the comparison; this measures how big that
    blind spot is.
    """
    lam_f = float(lam_source if lam_source is not None else _BG.LAM_SOURCE_BASELINE)
    bg = _BG.Background(lam_min=120.0, lam_max=lam_f + 5.0)
    lam = np.linspace(lam_min, lam_f, n_lambda)
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), float)
    u = np.ones_like(lam) if mode == "const" else \
        (lam - lam[0]) / (lam[-1] - lam[0])
    D = D * (1.0 + eps * u)
    resp = np.empty(n_lambda); resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    w = np.full(n_lambda, dlam); w[0] = w[-1] = 0.5 * dlam
    tail = np.cumsum((w / D ** 2)[::-1])[::-1]
    Wd = D ** 2 * tail
    gj = dlam * Wd[1:] / D[:-1] ** 4
    t4 = np.zeros(n_lambda); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    Hd = D ** 4 * t4
    builder = _DS.Sigma2Builder(background=bg, apply_c0=False)
    return ex.Grid(lam, dlam, D, resp, Wd, Hd, bg, builder)
