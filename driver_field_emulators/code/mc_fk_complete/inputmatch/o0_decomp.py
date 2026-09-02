"""Decompose the Order-0 deficit into (i) off-diagonal structure and
(ii) the lam < 406 truncation the MC imposes by starting at s(406)=0.

Everything below reads the SAME corr_op table.

Definitions
-----------
A(la)      = D(la)^4 C_00(n1,la; n2,la)          (smooth; no spline derivative)
H(la)      = int_la^{lam_f} D(t)^-2 dt            (so G(la) = D(la)^2 H(la))
Sigma2(la) = A'(la) / D(la)^4                     (driver_stats eq. 4)

FULL        = int_{lo}^{f} int_{lo}^{f} C_00(t1,t2) dt1 dt2      (true kernel)
COL_CONS    = same double integral of the SEPARABLE CAUSAL surrogate
              C_col(t1,t2) = [D(t1)D(t2)]^-2 A(min(t1,t2))
            = 2 int_{lo}^{f} A(t) H(t) D(t)^-2 dt
              (matches the true C EXACTLY on the diagonal, differs off it)
COL_MC      = int_{lo}^{f} Sigma2 G^2 dla = COL_CONS - A(lo) H(lo)^2
              (what driver_stats.order0_mc computes and what the MC realises,
               because the MC starts the ray at lam=406 with s=0)
"""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "..")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co

N1 = np.array([0.0, 0.0, 1.0])
LO = 406.0
LF = bgm.LAM_SOURCE_BASELINE

BG = bgm.Background()


def n2_of(g_arcmin):
    g = float(g_arcmin) / 60.0 * np.pi / 180.0
    return np.array([np.sin(g), 0.0, np.cos(g)])


def _H(la, ng=200):
    """H(la) = int_la^{LF} D(t)^-2 dt, vectorised, high order."""
    x, w = leggauss(ng)
    la = np.atleast_1d(np.asarray(la, float))
    out = np.empty_like(la)
    for i, l in enumerate(la):
        if l >= LF:
            out[i] = 0.0; continue
        t = 0.5 * (LF - l) * x + 0.5 * (LF + l)
        jw = 0.5 * (LF - l) * w
        out[i] = float(np.sum(jw / np.asarray(BG.D(t), float) ** 2))
    return out


def A_of(g_arcmin, la):
    n2 = n2_of(g_arcmin)
    la = np.atleast_1d(np.asarray(la, float))
    C = co.C_fn_batch(N1, la, n2, la)[:, 0, 0]
    return np.asarray(BG.D(la), float) ** 4 * C


def full_double(g_arcmin, ng=96, lo=LO):
    n2 = n2_of(g_arcmin)
    x, w = leggauss(ng)
    t2 = 0.5 * (LF - lo) * x + 0.5 * (LF + lo)
    w2 = 0.5 * (LF - lo) * w
    T1, T2, W = [], [], []
    for j, tb in enumerate(t2):
        t1 = 0.5 * (tb - lo) * x + 0.5 * (tb + lo)
        w1 = 0.5 * (tb - lo) * w
        T1.append(t1); T2.append(np.full_like(t1, tb)); W.append(w1 * w2[j])
    T1 = np.concatenate(T1); T2 = np.concatenate(T2); W = np.concatenate(W)
    C = co.C_fn_batch(N1, T1, n2, T2)[:, 0, 0]
    return 2.0 * float(np.sum(W * C))


def col_consistent(g_arcmin, ng=200, lo=LO):
    x, w = leggauss(ng)
    t = 0.5 * (LF - lo) * x + 0.5 * (LF + lo)
    jw = 0.5 * (LF - lo) * w
    A = A_of(g_arcmin, t)
    H = _H(t)
    D2 = np.asarray(BG.D(t), float) ** 2
    return 2.0 * float(np.sum(jw * A * H / D2))


def col_mc(g_arcmin, ng=200, lo=LO):
    """COL_CONS - A(lo) H(lo)^2 ; identical to int Sigma2 G^2 dla."""
    A0 = float(A_of(g_arcmin, [lo])[0])
    H0 = float(_H([lo])[0])
    return col_consistent(g_arcmin, ng=ng, lo=lo) - A0 * H0 ** 2, A0 * H0 ** 2


if __name__ == "__main__":
    print("# cross-check: by-parts COL_MC vs driver_stats.order0_mc (apply_c0=False)")
    b = ds.Sigma2Builder(apply_c0=False, lam_lo=LO, lam_hi=2328.0)
    for g in (1.0, 5.0, 17.3):
        cm, bnd = col_mc(g)
        ref = [ds.order0_mc(np.cos(g / 60 * np.pi / 180), builder=b, n_gauss=n)
               for n in (48, 96, 160, 240)]
        print(f"  g={g:6.2f}'  by-parts={cm:.6e}   order0_mc(ng=48,96,160,240)="
              + " ".join(f"{r:.5e}" for r in ref))

    print("\n# convergence of COL_CONS and FULL at 1'")
    for ng in (64, 128, 200, 320):
        print(f"  ng={ng:4d} COL_CONS={col_consistent(1.0, ng=ng):.8e}")
    for ng in (64, 96, 128, 192):
        print(f"  ng={ng:4d} FULL    ={full_double(1.0, ng=ng):.8e}")

    print("\n# decomposition")
    hdr = ("gamma", "FULL", "COL_CONS", "COL_MC", "offdiag", "trunc", "total")
    print("  {:>7} {:>13} {:>13} {:>13} {:>9} {:>9} {:>9}".format(*hdr))
    gam = [0.5, 1.0, 2.0, 2.6, 5.0, 6.7, 10.0, 17.3, 30.0, 44.4, 60.0, 90.0, 150.0]
    rows = []
    for g in gam:
        F = full_double(g, ng=128)
        Cc = col_consistent(g, ng=240)
        Cm, bnd = col_mc(g, ng=240)
        rows.append((g, F, Cc, Cm, Cc / F, Cm / Cc, Cm / F))
        print(f"  {g:7.2f} {F:13.6e} {Cc:13.6e} {Cm:13.6e} "
              f"{Cc/F:9.5f} {Cm/Cc:9.5f} {Cm/F:9.5f}")
    np.savez("o0_decomp.npz", rows=np.array(rows))
