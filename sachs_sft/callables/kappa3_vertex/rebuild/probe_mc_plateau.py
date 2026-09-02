"""Is the large-separation FK Monte-Carlo plateau a signal or a seed fluctuation?

The paper figure quotes, at every gamma above ~300 arcmin, FK_mc ~ -2.0e-6 with
a standard error of 7.3e-7, while the analytic FK has fallen to ~1e-8.  That
error bar is the scatter of EIGHT seeds divided by sqrt(8), and those same eight
seeds (1..8) are reused at every gamma, so the nine plotted points are not nine
independent measurements: they are one eight-seed mean, re-evaluated on a
covariance that barely changes once the rays have decorrelated.

A single seed at gamma = 957' returns +1.2e-6, i.e. the seed-to-seed spread is
about 2.1e-6, the size of the plateau itself.  This probe therefore draws
FRESH seeds and asks whether the mean moves to zero.

  * mean -> 0 with the expected 1/sqrt(n_seed) shrinkage  => the plateau is a
    fluctuation of seeds 1..8, and the figure needs more seeds, not a new term.
  * mean stays at -2e-6                                   => reproducible, and
    the estimator or the diagram needs a real explanation.

Why the FK estimator can scatter this much while the FF one does not: FK is
``T = sum_k W(k) Q[k] . Cov_sample(kappa_g, z z)`` and Q, the Levy vertex, is
~1e10.  Sampling noise in that covariance is amplified by ten orders of
magnitude.  The FF arm has no Q.  "Same machinery" does not mean "same error
scale", which is why FF sitting cleanly at 6e-7 is no evidence against a large
FK scatter.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[5]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
_FIX = Path(__file__).resolve().parents[1] / "equal_time_limber_cut15360_permaware"

def _worker(job):
    """One (gamma, seed) Monte-Carlo. Imports live here so each process is clean.

    Every setting travels inside ``job``: macOS spawns worker processes, which
    re-import this module rather than inheriting the parent's globals, so a
    module-level config would arrive unset.
    """
    gamma, seed, table, n_real, n_lambda, sigma = job
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(var, "2")
    sys.path.insert(0, str(_MC))
    sys.path.insert(0, str(_FIX))
    import driver_stats as ds
    import perm_aware_kappa3_callable as pa
    pa.TABLE_PATH = Path(table)
    pa._CACHE = None
    ds._k3 = pa
    import sachs_mc_core as mc
    cfg = mc.MCConfig(n_real=n_real, batch_size=3_000, n_lambda=n_lambda,
                      use_f_vertex=True, skew_scale=1.0,
                      sigma_lambda=sigma, seed=seed)
    t0 = time.time()
    r = mc.simulate_fk_vr(cfg, gamma)
    return gamma, seed, float(r.fk[0, 0]), float(r.fk_err[0, 0]), time.time() - t0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[293.9, 957.2, 1944.1])
    ap.add_argument("--seeds", type=int, nargs="+", default=list(range(9, 33)),
                    help="FRESH seeds; the figure already used 1..8")
    ap.add_argument("--n-real", type=int, default=24_000)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    table = str(a.table.resolve())
    jobs = [(g, s, table, a.n_real, a.n_lambda, a.sigma)
            for g in a.gammas for s in a.seeds]
    print(f"[probe] {len(jobs)} runs = {len(a.gammas)} gammas x {len(a.seeds)} "
          f"fresh seeds, n_real={a.n_real}, n_lambda={a.n_lambda}, "
          f"sigma={a.sigma}, {a.workers} workers", flush=True)

    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, (g, s, v, se, dt) in enumerate(ex.map(_worker, jobs), 1):
            rows.append((g, s, v, se))
            print(f"  [{i:3d}/{len(jobs)}] gamma={g:8.2f} seed={s:3d} "
                  f"fk={v:+.4e} (batch se {se:.1e})  {dt:.0f}s", flush=True)
    print(f"[probe] done in {(time.time()-t0)/60:.1f} min", flush=True)

    arr = np.array(rows)
    np.savez(a.out, gamma=arr[:, 0], seed=arr[:, 1], fk=arr[:, 2], se=arr[:, 3],
             n_real=a.n_real, n_lambda=a.n_lambda, sigma=a.sigma)

    print("\ngamma      n   mean          sem          spread(1sig)   mean/sem")
    for g in a.gammas:
        v = arr[arr[:, 0] == g][:, 2]
        m, sd = v.mean(), v.std(ddof=1)
        sem = sd / np.sqrt(len(v))
        print(f"{g:8.1f} {len(v):3d}  {m:+.4e}  {sem:.3e}  {sd:.3e}   "
              f"{m/sem:+.2f}")
    print(f"\n[probe] -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
