"""Downstream effect of the Sigma2 repair, part 2: the FK numbers themselves.

For each gamma, with the stock builder (baseline) and with ``R1Knots`` swapped
into ``grid.builder`` (apply_c0 matched):

  ref_true   fk_reference against the tabulated zeta6           -- builder-free,
             reported as a control that the grid/kernel side did not move
  ref_inj    fk_reference against the cumulant the deformation actually injects,
             ``zeta_injected(calib_matrix, Q)``                 -- builder-dependent
  exact      the exact expectation ``(T1+T2)+(T1+T2)^T`` of the estimator
             (only where --exact is given; it is the expensive one)

and the ratio of each to the paper's converged fold
``sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz``.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np

HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, HERE)
sys.path.insert(0, HERE + "/psd_fix")

import fk_expect_exact as ex
import sigma2_repaired as R
from run_gate2 import analytic_fk


def one(gam, sig, n_lambda, apply_c0, calibrate, builder_cls, want_exact):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=n_lambda, apply_anchor=apply_c0)
    if builder_cls is not None:
        grid.builder = builder_cls(background=grid.bg, apply_c0=apply_c0)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calibrate)
    M = ex.calib_matrix(grid, V, A, rho, sig, calibrate)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    ref_true = float(ex.fk_reference(grid, zt)[0, 3])
    ref_inj = float(ex.fk_reference(grid, ex.zeta_injected(M, Q))[0, 3])
    exact = np.nan
    if want_exact:
        t0 = time.time()
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        exact = float(((T1 + T2) + (T1 + T2).T)[0, 3])
        print(f"      [expectation {time.time()-t0:.0f}s]", flush=True)
    return ref_true, ref_inj, exact


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 5.0, 17.3])
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--apply-c0", type=int, default=0)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--exact-gammas", type=float, nargs="*", default=[])
    ap.add_argument("--out", default=HERE + "/psd_fix/_e2_downstream.npz")
    a = ap.parse_args()

    print(f"sigma={a.sigma} N={a.n_lambda} apply_c0={bool(a.apply_c0)} "
          f"calib={a.calibrate}")
    hdr = (f"{'gamma':>7} {'fold':>12} | {'ref_true':>12} {'/fold':>8} | "
           f"{'ref_inj base':>13} {'/fold':>8} | {'ref_inj R1':>13} {'/fold':>8} | "
           f"{'R1/base':>8}")
    print(hdr); print("-" * len(hdr))
    rows = []
    for gam in a.gammas:
        fold = float(analytic_fk([gam])[0])
        we = gam in a.exact_gammas
        rt_b, ri_b, ex_b = one(gam, a.sigma, a.n_lambda, bool(a.apply_c0),
                               a.calibrate, None, we)
        rt_r, ri_r, ex_r = one(gam, a.sigma, a.n_lambda, bool(a.apply_c0),
                               a.calibrate, R.R1Knots, we)
        assert abs(rt_b / rt_r - 1) < 1e-12, (rt_b, rt_r)
        print(f"{gam:>7.2f} {fold:>12.6e} | {rt_b:>12.6e} {rt_b/fold:>8.4f} | "
              f"{ri_b:>13.6e} {ri_b/fold:>8.4f} | {ri_r:>13.6e} {ri_r/fold:>8.4f} | "
              f"{ri_r/ri_b:>8.5f}", flush=True)
        if we:
            print(f"{'  exact':>7} {'':>12} | {'':>12} {'':>8} | "
                  f"{ex_b:>13.6e} {ex_b/fold:>8.4f} | {ex_r:>13.6e} {ex_r/fold:>8.4f} | "
                  f"{ex_r/ex_b:>8.5f}", flush=True)
        rows.append((gam, fold, rt_b, ri_b, ri_r, ex_b, ex_r))
    np.savez(a.out, rows=np.array(rows, float),
             cols=np.array(["gamma", "fold", "ref_true", "ref_inj_base",
                            "ref_inj_R1", "exact_base", "exact_R1"]))
    print("wrote", a.out)
