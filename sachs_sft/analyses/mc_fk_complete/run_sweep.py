"""Gate 3 / Gate 4 sweeps: sigma_lambda dependence and the gamma sweep.

For every point the exact expectation of the same estimator is computed too, so
the Monte-Carlo residual (sampling) is separated from the estimator's
finite-sigma_lambda residual.
"""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc
from run_gate2 import analytic_fk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+", default=[1.0])
    ap.add_argument("--sigmas", type=float, nargs="+", default=[4.0, 8.0, 16.0])
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--n-real", type=int, default=24_000)
    ap.add_argument("--batch", type=int, default=2_000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--seed-base", type=int, default=20260827)
    ap.add_argument("--no-exact", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    grid = ex.build_grid(n_lambda=a.n_lambda)
    rows = []
    print(f"N={a.n_lambda} dlam={grid.dlam:.3f} n_real={a.n_real} "
          f"seeds={a.seeds} calib={a.calibrate}")
    print(f"{'gamma':>7} {'sigma':>6} {'MC':>12} {'seedSEM':>10} {'scat%':>6} "
          f"{'exact':>12} {'MC/exact':>9} {'ref_inj':>12} {'analytic':>12} "
          f"{'MC/ana':>8} {'exact/ref':>9} {'vr/ana':>7} {'full/pert':>9}")
    for gam in a.gammas:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        ana = float(analytic_fk([gam])[0])
        zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
        ref_true = ex.fk_reference(grid, zt)[0, 3]
        for sig in a.sigmas:
            t0 = time.time()
            nodes = ex.node_stats(grid, cosg, sig, calibrate=a.calibrate)
            V, A, Q, rho = nodes
            M = ex.calib_matrix(grid, V, A, rho, sig, a.calibrate)
            ref_inj = ex.fk_reference(grid, ex.zeta_injected(M, Q))[0, 3]
            if a.no_exact:
                exact = np.nan
            else:
                T1x, T2x = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
                exact = ((T1x + T2x) + (T1x + T2x).T)[0, 3]
            vals, vv, vf = [], [], []
            for s in range(a.seeds):
                r = fc.simulate_fk_complete(gam, sigma_lambda=sig,
                                            n_lambda=a.n_lambda, n_real=a.n_real,
                                            batch_size=a.batch,
                                            seed=a.seed_base + 991 * s,
                                            grid=grid, nodes=nodes)
                vals.append(r.kk); vv.append(float(r.fk_vr_like[0, 0]))
                vf.append(float(r.fk_full[0, 0]))
            v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(v.size)
            rows.append(dict(gamma=gam, sigma=sig, mc=v.mean(), sem=sem,
                             scatter=v.std(ddof=1) / abs(v.mean()),
                             exact=exact, ref_inj=ref_inj, ref_true=ref_true,
                             analytic=ana, vr=np.mean(vv), full=np.mean(vf),
                             seeds=v))
            print(f"{gam:>7.2f} {sig:>6.1f} {v.mean():>12.5e} {sem:>10.2e} "
                  f"{100*v.std(ddof=1)/abs(v.mean()):>6.2f} {exact:>12.5e} "
                  f"{v.mean()/exact:>9.4f} {ref_inj:>12.5e} {ana:>12.5e} "
                  f"{v.mean()/ana:>8.4f} {exact/ref_inj:>9.4f} "
                  f"{np.mean(vv)/ana:>7.4f} {np.mean(vf)/v.mean():>9.4f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    np.savez(a.out, **{k: np.array([r[k] for r in rows]) for k in
                       ("gamma", "sigma", "mc", "sem", "scatter", "exact",
                        "ref_inj", "ref_true", "analytic", "vr", "full")},
             seeds=np.array([r["seeds"] for r in rows]),
             n_lambda=a.n_lambda, n_real=a.n_real, calibrate=a.calibrate)
    print(f"-> {a.out}")


if __name__ == "__main__":
    raise SystemExit(main())
