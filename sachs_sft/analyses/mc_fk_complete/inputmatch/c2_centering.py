"""C2. The in-sample centering, measured directly by pairing.

fk_expect_exact assumes the drive is  df[p] = 0.5 Q[p]:(z_p z_p - A[p])  with
A the EXACT AR(1) variance.  The production Monte-Carlo subtracts the batch
SAMPLE second moment instead.  Those are two different estimators, so
MC/exact != 1 by construction.  Run both on IDENTICAL fields (same rng seed,
same batch structure) and read the difference off directly.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np
HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/inputmatch")
import fk_expect_exact as ex      # noqa: E402
import mc_variant as mv           # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--gamma", type=float, default=1.0)
ap.add_argument("--sigma", type=float, default=8.0)
ap.add_argument("--batches", type=int, nargs="+",
                default=[250, 500, 1000, 2000, 4000, 6000])
ap.add_argument("--n-real", type=int, default=24_000)
ap.add_argument("--seeds", type=int, default=8)
ap.add_argument("--pilot-at", type=int, nargs="*", default=[500, 2000])
ap.add_argument("--out", required=True)
a = ap.parse_args()

cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
grid = ex.build_grid(n_lambda=1000)
nodes = ex.node_stats(grid, cosg, a.sigma, calibrate="smeared")
V, A, Q, rho = nodes
T1x, T2x = ex.expectation(grid, cosg, a.sigma, V=V, A=A, Q=Q, rho=rho)
exact = float(((T1x + T2x) + (T1x + T2x).T)[0, 3])
print(f"gamma={a.gamma}' sigma={a.sigma} n_real={a.n_real} exact={exact:.7e}")
print(f"{'batch':>6} {'n':>3} {'sample/exact':>13} {'A-cent/exact':>13} "
      f"{'delta=samp-A':>13} {'+/-':>10} {'delta [%]':>10} {'pilot-A [%]':>12}")
rows = []
t00 = time.time()
for mb in a.batches:
    ds, dA, dp = [], [], []
    for s in range(a.seeds):
        sd = 900_000_011 + 7919 * s + 13 * mb
        r1 = mv.run(a.gamma, sigma_lambda=a.sigma, n_real=a.n_real,
                    batch_size=mb, seed=sd, grid=grid, nodes=nodes,
                    center="sample")
        r2 = mv.run(a.gamma, sigma_lambda=a.sigma, n_real=a.n_real,
                    batch_size=mb, seed=sd, grid=grid, nodes=nodes,
                    center="exact")
        ds.append(r1["kk"]); dA.append(r2["kk"])
        if mb in a.pilot_at:
            r3 = mv.run(a.gamma, sigma_lambda=a.sigma, n_real=a.n_real,
                        batch_size=mb, seed=sd, grid=grid, nodes=nodes,
                        center="pilot")
            dp.append(r3["kk"])
    ds = np.array(ds); dA = np.array(dA)
    d = ds - dA
    se = d.std(ddof=1) / np.sqrt(d.size)
    ps = ""
    if dp:
        dp = np.array(dp); dd = dp - dA
        ps = (f"{100*dd.mean()/exact:>7.3f}+-"
              f"{100*dd.std(ddof=1)/np.sqrt(dd.size)/exact:.3f}")
    rows.append((mb, ds.mean()/exact, dA.mean()/exact, d.mean()/exact,
                 se/exact, (np.mean(dp)-dA.mean())/exact if len(dp) else np.nan))
    print(f"{mb:>6d} {a.seeds:>3d} {ds.mean()/exact:>13.5f} {dA.mean()/exact:>13.5f} "
          f"{d.mean()/exact:>13.5f} {se/exact:>10.5f} "
          f"{100*d.mean()/exact:>9.3f}% {ps:>12} [{time.time()-t00:.0f}s]",
          flush=True)
r = np.array(rows)
np.savez(a.out, batch=r[:, 0], samp=r[:, 1], acent=r[:, 2], delta=r[:, 3],
         dse=r[:, 4], pilot=r[:, 5], exact=exact, gamma=a.gamma, sigma=a.sigma,
         n_real=a.n_real, seeds=a.seeds)
# 1/batch fit
x = 1.0 / r[:, 0]
w = 1.0 / r[:, 4] ** 2
A1 = np.vstack([np.ones_like(x), x]).T
C = np.linalg.inv(A1.T @ (w[:, None] * A1))
b = C @ (A1.T @ (w * r[:, 3]))
print(f"\nweighted fit  delta/exact = c0 + c1/batch :")
print(f"   c0 = {b[0]:+.5f} +/- {np.sqrt(C[0,0]):.5f}   (must be 0 if the whole "
      f"offset is the in-sample centering)")
print(f"   c1 = {b[1]:+.4f} +/- {np.sqrt(C[1,1]):.4f}  -> at batch=2000: "
      f"{b[1]/2000:+.5f}")
print(f"-> {a.out}  [{time.time()-t00:.0f}s]")
