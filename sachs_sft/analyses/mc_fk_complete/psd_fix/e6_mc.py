"""Paired-seed Monte-Carlo: does the repair move the published marker value?

Same seeds, same grid, same n_real; only ``grid.builder`` differs.  Pairing
cancels almost all of the sampling noise, so the seed scatter of the DIFFERENCE
is a far tighter bound on the repair's effect than either run's own error bar.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np
HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, HERE); sys.path.insert(0, HERE + "/psd_fix")
import fk_expect_exact as ex
import fk_complete_core as fc
import sigma2_repaired as R
from run_gate2 import analytic_fk

ap = argparse.ArgumentParser()
ap.add_argument("--gamma", type=float, default=1.0)
ap.add_argument("--sigmas", type=float, nargs="+", default=[8.0, 4.0])
ap.add_argument("--seeds", type=int, default=12)
ap.add_argument("--n-lambda", type=int, default=1000)
ap.add_argument("--n-real", type=int, default=24_000)
ap.add_argument("--seed-base", type=int, default=20260828)
ap.add_argument("--out", default=HERE + "/psd_fix/_e6_mc.npz")
a = ap.parse_args()

cosg = float(np.cos(np.deg2rad(a.gamma / 60.0)))
fold = float(analytic_fk([a.gamma])[0])
store = {}
print(f"gamma={a.gamma}' N={a.n_lambda} n_real={a.n_real} seeds={a.seeds} "
      f"fold={fold:.6e}")
for sig in a.sigmas:
    vals = {}
    for tag, cls in (("base", None), ("R1", R.R1Knots)):
        grid = ex.build_grid(n_lambda=a.n_lambda)
        if cls is not None:
            grid.builder = cls(background=grid.bg, apply_c0=False)
        nodes = ex.node_stats(grid, cosg, sig, calibrate="smeared")
        v = []
        t0 = time.time()
        for s in range(a.seeds):
            r = fc.simulate_fk_complete(a.gamma, sigma_lambda=sig,
                                        n_lambda=a.n_lambda, n_real=a.n_real,
                                        batch_size=2000,
                                        seed=a.seed_base + 991 * s,
                                        grid=grid, nodes=nodes)
            v.append(r.kk)
        vals[tag] = np.array(v)
        sem = float(vals[tag].std(ddof=1) / np.sqrt(a.seeds))
        print(f"  sigma={sig:>4.1f} {tag:>5}: MC = {vals[tag].mean():.6e} "
              f"+/- {sem:.2e}   /fold = {vals[tag].mean()/fold:.4f} "
              f"+/- {sem/fold:.4f}   [{time.time()-t0:.0f}s]")
        store[f"{sig}_{tag}"] = vals[tag]
    d = vals["R1"] - vals["base"]
    dsem = d.std(ddof=1) / np.sqrt(a.seeds)
    print(f"  sigma={sig:>4.1f} PAIRED DIFF R1-base = {d.mean():+.4e} "
          f"+/- {dsem:.2e}   relative = {d.mean()/vals['base'].mean():+.3e} "
          f"+/- {dsem/vals['base'].mean():.3e}")
# paired sigma -> 0 extrapolation, per seed
if set(a.sigmas) >= {8.0, 4.0}:
    print()
    for tag in ("base", "R1"):
        e = 2 * store[f"4.0_{tag}"] - store[f"8.0_{tag}"]
        print(f"  extrapolated {tag:>5}: {e.mean():.6e} +/- "
              f"{e.std(ddof=1)/np.sqrt(a.seeds):.2e}   /fold = {e.mean()/fold:.4f}"
              f" +/- {e.std(ddof=1)/np.sqrt(a.seeds)/fold:.4f}")
    eb = 2 * store["4.0_base"] - store["8.0_base"]
    er = 2 * store["4.0_R1"] - store["8.0_R1"]
    d = er - eb
    print(f"  PAIRED extrapolated diff R1-base = {d.mean()/fold:+.5f} +/- "
          f"{d.std(ddof=1)/np.sqrt(a.seeds)/fold:.5f}  (in units of the fold)")
np.savez(a.out, fold=fold, **store)
print("wrote", a.out)
