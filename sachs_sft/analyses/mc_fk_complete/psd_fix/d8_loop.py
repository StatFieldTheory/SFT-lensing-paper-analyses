"""Are the Sigma2 ringing windows the SAME lambda that made the nominal
calibration blow up?  If so, the PSD defect and the 'lambda-grid graininess'
that forced the smeared calibration are one and the same problem."""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import fk_expect_exact as ex
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/psd_fix")
import sigma2_repaired as R

gam, sig, N = 17.3, 8.0, 1000
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
grid21 = np.load(__import__("corr_op_C_callable").TABLE_PATH,
                 allow_pickle=True)["lambda_grid"]

for name, cls in (("stock Sigma2", ds.Sigma2Builder), ("R1-repaired Sigma2", R.R1Knots)):
    g = ex.build_grid(n_lambda=N)
    g.builder = cls(background=g.bg, apply_c0=False)
    print(f"=== {name} ===")
    for calib in ("nominal", "smeared"):
        V, A, Q, rho = ex.node_stats(g, cosg, sig, calibrate=calib)
        qmax = np.abs(Q).reshape(N, -1).max(axis=1)
        spikes = np.where(qmax > 10 * np.median(qmax))[0]
        near = [float(np.min(np.abs(grid21 - g.lam[i]))) for i in spikes]
        M = ex.calib_matrix(g, V, A, rho, sig, calib)
        ref = ex.fk_reference(g, ex.zeta_injected(M, Q))[0, 3]
        T1, T2 = ex.expectation(g, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
        print(f"  {calib:>8}: max|Q| = {qmax.max():.3e}   spikes(>10x med) = "
              f"{spikes.size:3d}   median dist to a table node = "
              f"{np.median(near) if near else float('nan'):6.1f} Mpc   "
              f"FK/ref = {tot/ref:+.4f}")
    e = np.array([np.linalg.eigvalsh(core._assemble_6x6(g.builder, cosg, float(l))).min()
                  for l in g.lam])
    print(f"  non-PSD lambda: {(e < 0).sum()}/{N}")
