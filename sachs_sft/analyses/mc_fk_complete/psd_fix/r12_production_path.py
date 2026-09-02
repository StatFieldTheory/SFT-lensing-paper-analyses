"""Does the deformation actually inject the intended cumulant, node by node,
on the PRODUCTION FK path?

The claim under attack is that `fk_reference` is blind to the builder because
    cum3_from_Q(solve_Q(M, Z), M) == sym(Z)   for any M.
That identity requires M to be positive definite.  Here it is checked node by
node on the real production configuration -- build_grid(n_lambda=1000),
node_stats(sigma_lambda=8), both the "smeared" (production) and "nominal"
calibrations -- together with the PSD status of V, A and the calibration matrix.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R
import fk_expect_exact as ex

def sym6(Z):
    return (Z + Z.transpose(0,2,1) + Z.transpose(1,0,2) + Z.transpose(1,2,0)
            + Z.transpose(2,0,1) + Z.transpose(2,1,0))/6.0

def npsd(X):
    n = 0; worst = np.inf
    for M in X:
        w = np.linalg.eigvalsh(0.5*(M+M.T))
        worst = min(worst, w.min()/max(abs(w).max(), 1e-300))
        if w.min() < -1e-14*abs(w).max(): n += 1
    return n, worst

GAM = 1.0
SIG = 8.0
cos = float(np.cos(np.deg2rad(GAM/60.0)))
for calib in ("smeared", "nominal"):
    print(f"\n########  calibrate={calib!r}, gamma={GAM}', sigma_lambda={SIG}, "
          f"n_lambda=1000  ########")
    for bname, cls in (("stock", None), ("R1", R.R1Knots)):
        g = ex.build_grid(n_lambda=1000)
        if cls is not None:
            g.builder = cls(background=g.bg, apply_c0=False)
        V, A, Q, rho = ex.node_stats(g, cos, SIG, calibrate=calib)
        M = ex.calib_matrix(g, V, A, rho, SIG, calib)
        Z = np.array([ds.zeta6(cos, float(l)) for l in g.lam])
        inj = np.array([ds.cum3_from_Q(Q[k], M[k]) for k in range(g.lam.size)])
        Zs = np.array([sym6(z) for z in Z])
        den = np.maximum(np.abs(Zs).max(axis=(1,2,3)), 1e-300)
        err = np.abs(inj - Zs).max(axis=(1,2,3))/den
        nV, wV = npsd(V); nA, wA = npsd(A); nM, wM = npsd(M)
        bad = np.where(err > 1e-6)[0]
        print(f"  [{bname}] non-PSD: V {nV}/1000 (min/max {wV:+.2e}) | "
              f"A {nA}/1000 ({wA:+.2e}) | calib M {nM}/1000 ({wM:+.2e})")
        print(f"          identity cum3_from_Q(Q,M) vs sym(zeta): "
              f"median {np.median(err):.2e}  max {err.max():.2e}  "
              f"nodes with err>1e-6: {bad.size}")
        for k in bad[:6]:
            print(f"            node {k:>4} lam={g.lam[k]:8.2f} err={err[k]:.3e} "
                  f"|Q|max={np.abs(Q[k]).max():.3e} "
                  f"M_00={M[k,0,0]:+.4e} mineig/max="
                  f"{np.linalg.eigvalsh(0.5*(M[k]+M[k].T)).min()/abs(np.linalg.eigvalsh(0.5*(M[k]+M[k].T))).max():+.3e}")
        fkr = ex.fk_reference(g, inj)
        print(f"          fk_reference[against injected cumulant] kk = {fkr[0,0]:.10e}")
