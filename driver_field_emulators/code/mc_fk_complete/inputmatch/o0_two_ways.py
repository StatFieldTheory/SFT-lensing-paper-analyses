"""Order-0 <kappa kappa> computed TWO ways from the SAME corr_op table.

(a) full  : int int C_00(n1,t1; n2,t2) dt1 dt2   -- genuine off-diagonal kernel
(b) collapsed : int Sigma2_00(la) G(la)^2 dla    -- what the MC's white-noise
                driver realises, with Sigma2 = d/dlam[D^4 C(lam,lam)]/D^4.

(b) is EXACTLY the double integral of the *separable causal* surrogate
    C_col(t1,t2) = D(t1)^-2 D(t2)^-2 int_0^{min} D^4 Sigma2 dla ,
which reproduces the true C on the diagonal by construction but not off it.
"""
import sys, time
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, "..")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co

N1 = np.array([0.0, 0.0, 1.0])
LAM_LO = 406.0
LAM_F = bgm.LAM_SOURCE_BASELINE   # 2313.029


def n2_of(gamma_arcmin):
    g = float(gamma_arcmin) / 60.0 * np.pi / 180.0
    return np.array([np.sin(g), 0.0, np.cos(g)])


def o0_full(gamma_arcmin, ng=48, lam_lo=LAM_LO, lam_f=LAM_F, comp=(0, 0)):
    """int int C_ab dt1 dt2 via 2x lower-triangle (C_00 is t1<->t2 symmetric)."""
    n2 = n2_of(gamma_arcmin)
    x, w = leggauss(int(ng))
    # outer t2 over [lam_lo, lam_f]
    t2 = 0.5 * (lam_f - lam_lo) * x + 0.5 * (lam_f + lam_lo)
    w2 = 0.5 * (lam_f - lam_lo) * w
    T1, T2, W = [], [], []
    for j, tb in enumerate(t2):
        t1 = 0.5 * (tb - lam_lo) * x + 0.5 * (tb + lam_lo)
        w1 = 0.5 * (tb - lam_lo) * w
        T1.append(t1); T2.append(np.full_like(t1, tb)); W.append(w1 * w2[j])
    T1 = np.concatenate(T1); T2 = np.concatenate(T2); W = np.concatenate(W)
    C = co.C_fn_batch(N1, T1, n2, T2)[:, comp[0], comp[1]]
    return 2.0 * float(np.sum(W * C))


def o0_collapsed(gamma_arcmin, ng=48, lam_lo=LAM_LO, lam_f=LAM_F, builder=None):
    b = builder
    x, w = leggauss(int(ng))
    la = 0.5 * (lam_f - lam_lo) * x + 0.5 * (lam_f + lam_lo)
    jw = 0.5 * (lam_f - lam_lo) * w
    sig = np.array([b.matrix(np.cos(float(gamma_arcmin) / 60.0 * np.pi / 180.0),
                             float(l))[0, 0] for l in la])
    G = ds.order0_window(b.bg, la, lam_f, n_gauss=ng)
    return float(np.sum(jw * sig * G ** 2))


if __name__ == "__main__":
    b = ds.Sigma2Builder(apply_c0=False, lam_lo=LAM_LO, lam_hi=2328.0)
    print("# convergence of the double integral at gamma=1'")
    for ng in (24, 32, 48, 64, 96, 128):
        t0 = time.time()
        v = o0_full(1.0, ng=ng)
        print(f"  ng={ng:4d}  O0_full = {v:.8e}   ({time.time()-t0:.1f}s)")
    print("\n# convergence of the collapsed integral at gamma=1'")
    for ng in (24, 32, 48, 64, 96):
        v = o0_collapsed(1.0, ng=ng, builder=b)
        print(f"  ng={ng:4d}  O0_col  = {v:.8e}")

    print("\n# gamma dependence  (ng=96 full, ng=64 collapsed)")
    gam = [0.5, 1.0, 2.0, 2.6, 5.0, 6.7, 10.0, 17.3, 30.0, 44.4, 90.0, 150.0]
    rows = []
    for g in gam:
        a = o0_full(g, ng=96)
        c = o0_collapsed(g, ng=64, builder=b)
        rows.append((g, a, c, c / a))
        print(f"  gamma={g:7.2f}'  full={a:+.6e}  collapsed={c:+.6e}  "
              f"col/full={c/a:.6f}")
    np.savez("o0_two_ways.npz", gamma=np.array([r[0] for r in rows]),
             full=np.array([r[1] for r in rows]),
             collapsed=np.array([r[2] for r in rows]),
             ratio=np.array([r[3] for r in rows]))
