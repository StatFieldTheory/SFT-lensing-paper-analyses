"""FF as an explicit QUADRATIC functional of the state two-point function C.

The FF moment measured by ``sachs_mc_core.simulate_ff_crn`` is
``<kk>_Fon - <kk>_Foff``.  For a Gaussian driving field this is, through O(F^2),

    FF_AB = <k2_A k2_B> + <k1_A k3_B> + <k3_A k1_B>

with (continuum, R(t,la) = [D(la)/D(t)]^2, W(la) = int_la^{lf} R dt)

    k1_A(n) = int W(la) f_A(n,la) dla                       -> <k1 s1> = int C dt
    k2_A(n) = int W(la) Fs_Abc s1_b(n,la) s1_c(n,la) dla
    k3_A(n) = 2 int W(la) Fs_Abc s1_b(n,la) s2_c(n,la) dla
    s2_c(n,la) = int_0^la R(la,v) Fs_cfg s1_f(n,v) s1_g(n,v) dv

EVERY term pairs exactly FOUR s1's, i.e. it is exactly quadratic in the state
two-point function C_ij(n,la; n',mu) = <s1_i(n,la) s1_j(n',mu)>.  So FF is a
fixed quadratic functional and the only question is WHICH C the two sides feed
it.  This module evaluates that one functional for any C.

Three kernels, all built from the SAME corr_op table:
  full : C_ij(n,la; n',mu) straight from corr_op         (the formalism)
  col  : [D(la)D(mu)]^-2 A(min(la,mu)),  A = D^4 C(t,t)  (equal-time collapse;
         matches the true C EXACTLY on the diagonal, differs off it)
  mc   : col with A -> A - A(lam_min)                    (what the MC realises:
         it starts the ray at lam_min = 406 with s = 0)
"""
from __future__ import annotations

import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "..")
import _bootstrap
_DS, _BGM, _CORE = _bootstrap.wire()
import corr_op_C_callable as _co

N1 = np.array([0.0, 0.0, 1.0])
LO = 406.0
LF = _BGM.LAM_SOURCE_BASELINE
BG = _BGM.Background(lam_min=120.0, lam_max=LF + 5.0)

# symmetric F vertex: F(s)_a = Fs_abc s_b s_c
_F = _CORE._F
FS = 0.5 * (_F + _F.transpose(0, 2, 1))


def n2_of(g_arcmin):
    g = float(g_arcmin) / 60.0 * np.pi / 180.0
    return np.array([np.sin(g), 0.0, np.cos(g)])


def grid(n):
    lam = np.linspace(LO, LF, n)
    h = lam[1] - lam[0]
    w = np.full(n, h); w[0] = w[-1] = 0.5 * h
    D = np.asarray(BG.D(lam), float)
    # H(la) = int_la^{LF} D^-2 dt  (high-order GL, then W = D^2 H)
    x, gw = leggauss(120)
    H = np.empty(n)
    for i, l in enumerate(lam):
        if l >= LF:
            H[i] = 0.0; continue
        t = 0.5 * (LF - l) * x + 0.5 * (LF + l)
        H[i] = float(np.sum(0.5 * (LF - l) * gw / np.asarray(BG.D(t), float) ** 2))
    W = D ** 2 * H
    # lower-triangular trapezoid weights: Wtri[i, j] integrates v over [lam_0, lam_i]
    Wtri = np.zeros((n, n))
    for i in range(n):
        if i == 0:
            continue
        Wtri[i, :i + 1] = h
        Wtri[i, 0] = 0.5 * h
        Wtri[i, i] = 0.5 * h
    R = np.where(np.arange(n)[None, :] <= np.arange(n)[:, None],
                 (D[None, :] / D[:, None]) ** 2, 0.0)      # R[mu, v]
    return lam, h, w, D, W, Wtri, R


def _C_pairs(na, nb, lam):
    """C_ij(na, lam_p; nb, lam_q) as (n, n, 3, 3)."""
    n = lam.size
    P = np.repeat(lam, n)
    Q = np.tile(lam, n)
    return _co.C_fn_batch(na, P, nb, Q).reshape(n, n, 3, 3)


def kernels(g_arcmin, lam, D, kind):
    """Return (Cx, Cxx, Cyy) for the requested kernel kind."""
    n2 = n2_of(g_arcmin)
    if kind == "full":
        Cx = _C_pairs(N1, n2, lam)
        Cxx = _C_pairs(N1, N1, lam)
        Cyy = _C_pairs(n2, n2, lam)
        return Cx, Cxx, Cyy
    # collapsed surrogates from the SAME table's equal-time diagonal
    D4 = D ** 4
    Ax = D4[:, None, None] * _co.C_fn_batch(N1, lam, n2, lam)      # (n,3,3)
    Aw1 = D4[:, None, None] * _co.C_fn_batch(N1, lam, N1, lam)
    Aw2 = D4[:, None, None] * _co.C_fn_batch(n2, lam, n2, lam)
    if kind == "mc":
        Ax = Ax - Ax[0]; Aw1 = Aw1 - Aw1[0]; Aw2 = Aw2 - Aw2[0]
    elif kind != "col":
        raise ValueError(kind)
    n = lam.size
    mn = np.minimum(np.arange(n)[:, None], np.arange(n)[None, :])
    pref = 1.0 / (D[:, None] ** 2 * D[None, :] ** 2)
    Cx = pref[:, :, None, None] * Ax[mn]
    Cxx = pref[:, :, None, None] * Aw1[mn]
    Cyy = pref[:, :, None, None] * Aw2[mn]
    return Cx, Cxx, Cyy


def ff_terms(g_arcmin, n=241, kind="full"):
    lam, h, w, D, W, Wtri, R = grid(n)
    Cx, Cxx, Cyy = kernels(g_arcmin, lam, D, kind)
    Cw1 = np.einsum("ppij->pij", Cxx)          # within-ray-1 equal time
    Cw2 = np.einsum("ppij->pij", Cyy)
    wW = w * W

    # --- <k2_A k2_B> : connected + disconnected -------------------------
    t1_conn = 2.0 * np.einsum("p,q,Abc,Bde,pqbd,pqce->AB",
                              wW, wW, FS, FS, Cx, Cx, optimize=True)
    m2A = np.einsum("p,Abc,pbc->A", wW, FS, Cw1, optimize=True)   # <k2(n1)>
    m2B = np.einsum("p,Bde,pde->B", wW, FS, Cw2, optimize=True)   # <k2(n2)>
    t1_disc = np.outer(m2A, m2B)

    # --- Psi: <k1 s1> ----------------------------------------------------
    PsiA = np.einsum("t,tqAd->qAd", w, Cx, optimize=True)    # <k1_A(n1) s1_d(n2,q)>
    PsiB = np.einsum("t,qtdB->qdB", w, Cx, optimize=True)    # <k1_B(n2) s1_d(n1,q)>

    # --- <k1_A(n1) k3_B(n2)> --------------------------------------------
    RW = Wtri * R                                             # [mu, v]
    mean_s2_ray2 = np.einsum("mv,efg,vfg->me", RW, FS, Cw2, optimize=True)
    p1 = 2.0 * np.einsum("m,mAd,Bde,me->AB", wW, PsiA, FS, mean_s2_ray2, optimize=True)
    p2 = 4.0 * np.einsum("m,mv,Bde,efg,vAf,mvdg->AB",
                         wW, RW, FS, FS, PsiA, Cyy, optimize=True)
    term2 = p1 + p2

    # --- <k3_A(n1) k1_B(n2)> --------------------------------------------
    mean_s2_ray1 = np.einsum("mv,efg,vfg->me", RW, FS, Cw1, optimize=True)
    q1 = 2.0 * np.einsum("m,mdB,Ade,me->AB", wW, PsiB, FS, mean_s2_ray1, optimize=True)
    q2 = 4.0 * np.einsum("m,mv,Ade,efg,vfB,mvdg->AB",
                         wW, RW, FS, FS, PsiB, Cxx, optimize=True)
    term3 = q1 + q2

    return dict(conn22=t1_conn, disc=t1_disc, k1k3=term2, k3k1=term3,
                moment=t1_conn + t1_disc + term2 + term3,
                connected=t1_conn + term2 + term3)


def o0_from_kernel(g_arcmin, n=241, kind="full"):
    lam, h, w, D, W, Wtri, R = grid(n)
    Cx, _, _ = kernels(g_arcmin, lam, D, kind)
    return float(np.einsum("p,q,pq->", w, w, Cx[:, :, 0, 0], optimize=True))
