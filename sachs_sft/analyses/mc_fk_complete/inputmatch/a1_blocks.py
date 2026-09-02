"""A1. The decisive test of the 'per-block common offset'.

Identical configuration everywhere; the ONLY thing that differs between blocks
is the seed base.  Two seeding schemes are run on the same block bases:

  int    seeds = base + 991*s          (production: arithmetic progression)
  spawn  SeedSequence(base).spawn()[s] (guaranteed independent streams)

Per seed we keep the 12 per-batch values, so between-block / between-seed /
within-seed variance can be separated with the right dof.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np

HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/inputmatch")
import fk_expect_exact as ex          # noqa: E402
import mc_variant as mv               # noqa: E402
from run_gate2 import analytic_fk     # noqa: E402

BASES = [11_000_003, 22_000_013, 33_000_037, 44_000_051,
         55_000_069, 66_000_083, 77_000_097, 88_000_121]

ap = argparse.ArgumentParser()
ap.add_argument("--scheme", required=True, choices=["int", "spawn"])
ap.add_argument("--gamma", type=float, default=1.0)
ap.add_argument("--sigmas", type=float, nargs="+", default=[8.0, 4.0])
ap.add_argument("--n-lambda", type=int, default=1000)
ap.add_argument("--n-real", type=int, default=24_000)
ap.add_argument("--batch", type=int, default=2_000)
ap.add_argument("--blocks", type=int, default=8)
ap.add_argument("--seeds", type=int, default=16)
ap.add_argument("--out", required=True)
a = ap.parse_args()

cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
grid = ex.build_grid(n_lambda=a.n_lambda)
ana = float(analytic_fk([a.gamma])[0])
res = {}
t00 = time.time()
for sig in a.sigmas:
    nodes = ex.node_stats(grid, cosg, sig, calibrate="smeared")
    V, A, Q, rho = nodes
    T1x, T2x = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    exact = float(((T1x + T2x) + (T1x + T2x).T)[0, 3])
    vals = np.empty((a.blocks, a.seeds))
    bat = np.empty((a.blocks, a.seeds, a.n_real // a.batch))
    ngood = np.empty((a.blocks, a.seeds), int)
    for b, base in enumerate(BASES[:a.blocks]):
        for s in range(a.seeds):
            rng = mv.make_rng(seed=base + 991 * s, base=base, index=s,
                              mode=a.scheme, nmax=a.seeds)
            r = mv.run(a.gamma, sigma_lambda=sig, n_lambda=a.n_lambda,
                       n_real=a.n_real, batch_size=a.batch, grid=grid,
                       nodes=nodes, rng=rng)
            vals[b, s] = r["kk"]; bat[b, s] = r["kk_batches"]
            ngood[b, s] = r["n_good"].sum()
        print(f"  [{a.scheme} sig={sig}] block {b} base={base} "
              f"mean/exact={vals[b].mean()/exact:.4f} "
              f"sd={vals[b].std(ddof=1)/exact:.4f} [{time.time()-t00:.0f}s]",
              flush=True)
    res[f"vals_{sig}"] = vals; res[f"bat_{sig}"] = bat
    res[f"ngood_{sig}"] = ngood; res[f"exact_{sig}"] = exact
np.savez(a.out, analytic=ana, bases=np.array(BASES[:a.blocks]),
         scheme=a.scheme, gamma=a.gamma, n_real=a.n_real, batch=a.batch,
         n_lambda=a.n_lambda, **res)
print(f"-> {a.out}  [{time.time()-t00:.0f}s]")
