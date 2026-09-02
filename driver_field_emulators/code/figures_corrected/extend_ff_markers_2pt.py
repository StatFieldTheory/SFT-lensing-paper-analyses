"""Restore the FF Monte-Carlo markers at the original figure's full gamma grid.

Restricting the validation figure's marker list to small separations is a
statement about the FK channel: its estimator depends on the AR(1) regulator
sigma_lambda, so beyond a few arcminutes it measures the regulator rather
than the diagram.  The FF channel has no such limitation -- it is Gaussian,
carries no deformation tensor, and its large-separation plateau agreeing with
the analytic fold is one of the published figure's checks -- so the FF
markers keep the original figure's full gamma grid.

This script runs ``simulate_ff_crn`` (identical seed and settings to the
generator) at the gammas the restricted run dropped, merges them into the
generator's cache with the FK entries set to NaN (matplotlib draws nothing
there), and replots through the generator's own ``_plot``.
"""

from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"

# the original generator grid beyond the restricted FK range
GAMMAS_LARGE = (1.0, 2.6)


def _ff_one(job):
    gamma, ff_real, ff_nl = job
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(var, "2")
    sys.path.insert(0, str(_MC))
    import sachs_mc_core as mc
    cfg = mc.MCConfig(n_real=ff_real, batch_size=3_000, n_lambda=ff_nl,
                      use_f_vertex=True, apply_anchor=False, seed=11)
    r = mc.simulate_ff_crn(cfg, gamma)
    return gamma, float(r.ff_moment[0, 0]), float(r.ff_moment_err[0, 0])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ff-real", type=int, default=48_000)
    ap.add_argument("--ff-nl", type=int, default=4000)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    cache = _MC / "outputs" / "appendix_mc_curve.npz"
    D = dict(np.load(cache, allow_pickle=True))
    print(f"[ff-ext] cache has markers at {np.asarray(D['g_mc']).round(1)}",
          flush=True)

    jobs = [(g, a.ff_real, a.ff_nl) for g in GAMMAS_LARGE]
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for g, v, se in ex.map(_ff_one, jobs):
            rows.append((g, v, se))
            print(f"  gamma={g:8.1f}  FF_mc={v:+.3e} +/-{se:.1e}", flush=True)

    g_new = np.array([r[0] for r in rows])
    ff_new = np.array([r[1] for r in rows])
    se_new = np.array([r[2] for r in rows])
    g_all = np.concatenate([np.asarray(D["g_mc"], float), g_new])
    order = np.argsort(g_all)
    nan = np.full(g_new.size, np.nan)
    D["g_mc"] = g_all[order]
    D["fk_mc"] = np.concatenate([np.asarray(D["fk_mc"], float), nan])[order]
    D["fk_se"] = np.concatenate([np.asarray(D["fk_se"], float), nan])[order]
    D["ffc_mc"] = np.concatenate([np.asarray(D["ffc_mc"], float), ff_new])[order]
    D["ffc_se"] = np.concatenate([np.asarray(D["ffc_se"], float), se_new])[order]
    np.savez(cache, **D)
    print(f"[ff-ext] cache -> {cache}", flush=True)

    print("[ff-ext] done (no replot)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
