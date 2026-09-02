"""Adversarial runner: independent seed blocks, recording per-seed value AND
the within-seed batch SE, so between-seed and within-seed variance can be
compared directly.  Writes to review/stats/ only.
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))   # mc_fk_complete/
import fk_expect_exact as ex          # noqa: E402
import fk_complete_core as fc         # noqa: E402
from run_gate2 import analytic_fk     # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--n-real", type=int, default=24_000)
    ap.add_argument("--batch", type=int, default=2_000)
    ap.add_argument("--seeds", type=int, default=32)
    ap.add_argument("--seed-base", type=int, required=True)
    ap.add_argument("--seed-step", type=int, default=991)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
    grid = ex.build_grid(n_lambda=a.n_lambda)
    nodes = ex.node_stats(grid, cosg, a.sigma, calibrate=a.calibrate)
    V, A, Q, rho = nodes
    T1x, T2x = ex.expectation(grid, cosg, a.sigma, V=V, A=A, Q=Q, rho=rho)
    exact = float(((T1x + T2x) + (T1x + T2x).T)[0, 3])
    ana = float(analytic_fk([a.gamma])[0])

    vals, bse, t1s, t2s, vrs = [], [], [], [], []
    t00 = time.time()
    for s in range(a.seeds):
        r = fc.simulate_fk_complete(a.gamma, sigma_lambda=a.sigma,
                                    n_lambda=a.n_lambda, n_real=a.n_real,
                                    batch_size=a.batch,
                                    seed=a.seed_base + a.seed_step * s,
                                    grid=grid, nodes=nodes,
                                    calibrate=a.calibrate)
        vals.append(r.kk); bse.append(float(r.fk_batch_se[0, 0]))
        t1s.append(float(r.t1[0, 3])); t2s.append(float(r.t2[0, 3]))
        vrs.append(float(r.fk_vr_like[0, 0]))
    v = np.array(vals); b = np.array(bse)
    sem = v.std(ddof=1) / np.sqrt(v.size)
    print(f"gamma={a.gamma} sigma={a.sigma} N={a.n_lambda} n_real={a.n_real} "
          f"batch={a.batch} seeds={a.seeds} base={a.seed_base} "
          f"[{time.time()-t00:.0f}s]")
    print(f"  exact={exact:.6e}  analytic={ana:.6e}")
    print(f"  MC   ={v.mean():.6e} +/- {sem:.3e}   MC/exact={v.mean()/exact:.4f}"
          f" +/- {sem/exact:.4f}   MC/analytic={v.mean()/ana:.4f} +/- {sem/ana:.4f}")
    print(f"  seed sd={v.std(ddof=1):.4e}   mean within-seed batchSE="
          f"{b.mean():.4e}   ratio(seed sd / batchSE)={v.std(ddof=1)/b.mean():.3f}")
    np.savez(a.out, gamma=a.gamma, sigma=a.sigma, n_lambda=a.n_lambda,
             n_real=a.n_real, batch=a.batch, seed_base=a.seed_base,
             seed_step=a.seed_step, seeds=v, batch_se=b, t1=np.array(t1s),
             t2=np.array(t2s), vr=np.array(vrs), exact=exact, analytic=ana)
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
