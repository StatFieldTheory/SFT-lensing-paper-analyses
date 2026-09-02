"""Gate 2 -- one point.  gamma = 1 arcmin, production settings, many seeds.

Reports the Monte-Carlo against three references, so the residual can be
attributed:

  exact       the closed-form expectation of the SAME estimator (MC error only)
  ref_inj     the local-limit FK built from the cumulant actually injected
              (adds the finite-sigma_lambda error of the estimator)
  analytic    the paper's converged fold, products/table_..._permfix_xi.npz
              (adds the vertex table's leg asymmetry and the discretisation)
"""
from __future__ import annotations
import argparse, math, sys, time
from pathlib import Path
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc
import _bootstrap


def analytic_fk(gammas):
    """kappa-kappa order-2 entries of the converged permutation-closed fold."""
    d = np.load(_bootstrap.FOLD, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return np.interp(gammas, g[o], v[o])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--n-real", type=int, default=24_000)
    ap.add_argument("--batch", type=int, default=2_000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
    grid = ex.build_grid(n_lambda=a.n_lambda)
    nodes = ex.node_stats(grid, cosg, a.sigma, calibrate=a.calibrate)
    V, A, Q, rho = nodes
    M = ex.calib_matrix(grid, V, A, rho, a.sigma, a.calibrate)
    zi = ex.zeta_injected(M, Q)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    ref_inj = ex.fk_reference(grid, zi)[0, 3]
    ref_true = ex.fk_reference(grid, zt)[0, 3]
    T1x, T2x = ex.expectation(grid, cosg, a.sigma, V=V, A=A, Q=Q, rho=rho)
    exact = ((T1x + T2x) + (T1x + T2x).T)[0, 3]
    ana = float(analytic_fk([a.gamma])[0])

    print(f"gamma={a.gamma}' sigma={a.sigma} N={a.n_lambda} dlam={grid.dlam:.3f} "
          f"sigma/dlam={a.sigma/grid.dlam:.2f} calib={a.calibrate} "
          f"n_real={a.n_real} seeds={a.seeds}")
    print(f"  analytic FK (paper fold)      {ana:.6e}")
    print(f"  local ref, tabulated zeta     {ref_true:.6e}  "
          f"({ref_true/ana:.4f} x analytic)")
    print(f"  local ref, injected  zeta     {ref_inj:.6e}  "
          f"({ref_inj/ref_true:.4f} injection fidelity)")
    print(f"  exact expectation of estimator{exact:.6e}  "
          f"({exact/ref_inj:.4f} finite-sigma)")
    print()

    vals, vals_vr, vals_full = [], [], []
    for s in range(a.seeds):
        t0 = time.time()
        r = fc.simulate_fk_complete(a.gamma, sigma_lambda=a.sigma,
                                    n_lambda=a.n_lambda, n_real=a.n_real,
                                    batch_size=a.batch, seed=20260827 + 991 * s,
                                    grid=grid, nodes=nodes)
        vals.append(r.kk); vals_vr.append(float(r.fk_vr_like[0, 0]))
        vals_full.append(float(r.fk_full[0, 0]))
        print(f"  seed {s}: MC={r.kk:.6e}  /exact={r.kk/exact:.4f} "
              f"/analytic={r.kk/ana:.4f}  batchSE={r.fk_batch_se[0,0]:.2e} "
              f"n_good={r.n_good} [{time.time()-t0:.0f}s]", flush=True)
    v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(v.size)
    print()
    print(f"  MC mean          {v.mean():.6e} +/- {sem:.2e} (seed SEM, {a.seeds} seeds)")
    print(f"  seed scatter     {v.std(ddof=1)/abs(v.mean()):.2%}")
    print(f"  MC / exact       {v.mean()/exact:.4f} +/- {sem/abs(exact):.4f}")
    print(f"  MC / ref_inj     {v.mean()/ref_inj:.4f} +/- {sem/abs(ref_inj):.4f}")
    print(f"  MC / analytic    {v.mean()/ana:.4f} +/- {sem/abs(ana):.4f}")
    vf = np.array(vals_full)
    print(f"  nonlinear-arm flavour / perturbative = {vf.mean()/v.mean():.4f} "
          f"(F^2+ contamination)")
    vv = np.array(vals_vr)
    print(f"  T1-only (fk_vr structure) / analytic = {vv.mean()/ana:.4f}")
    if a.out:
        np.savez(a.out, gamma=a.gamma, sigma=a.sigma, n_lambda=a.n_lambda,
                 n_real=a.n_real, mc=v, mc_vr=vv, mc_full=vf, exact=exact,
                 ref_inj=ref_inj, ref_true=ref_true, analytic=ana)
        print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
