"""One real FF sweep at a single source distance, in its own process.

The FF order-2 sweep is expensive and deadlocks under loky, so it is run with
sweep.n_jobs = 1. Serially that is one core for five source distances; instead
each distance gets its own single-threaded PROCESS, which is the same
parallelism without loky. Each process also gets PRIVATE expansion and
propagator caches, so concurrent runs cannot race on the shared cache
directory inside the production run folder (the production output npz is never
touched: run_vertex_sweep redirects output.path to a private scratch file).

Usage::

    <sft-wick python> run_ff_one_lambda.py --lam 2094.894 --tag z1p7
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys
import time
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO / "talk" / "scripts"))
import run_multiz_components as R  # noqa: E402

PRODUCTS = _REPO / "driver_field_emulators" / "products"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--tag", type=str, required=True)
    a = ap.parse_args()

    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "2"

    cache_root = PRODUCTS / "_ff_caches" / a.tag
    cache_root.mkdir(parents=True, exist_ok=True)

    # give this process private caches by wrapping the sweep entry point
    real_sweep = R.run_vertex_sweep
    from sft_wick.workflow.config import load_workflow_config  # noqa: E402

    original_load = load_workflow_config

    def patched_load(path):
        cfg = original_load(path)
        exp = dataclasses.replace(
            cfg.expand, cache_path=str(cache_root / "expand"))
        props = dataclasses.replace(
            cfg.propagators, cache_path=str(cache_root / "propagators"))
        return dataclasses.replace(cfg, expand=exp, propagators=props)

    import sft_wick.workflow.config as _cfgmod
    _cfgmod.load_workflow_config = patched_load

    t0 = time.time()
    print(f"[{a.tag}] lam_source = {a.lam}, private caches at {cache_root}",
          flush=True)
    gamma, xi = real_sweep(R.CONFIGS["ff"], np.array([a.lam]), f"ff_{a.tag}",
                           n_jobs=1)
    ff = np.asarray(xi[0], float)
    out = PRODUCTS / f"multiz_ff_real_{a.tag}.npz"
    np.savez(out, lam=a.lam, gamma=gamma, ff=ff)
    print(f"[{a.tag}] done in {(time.time()-t0)/60:.1f} min; "
          f"xi(0.5')={ff[0]:+.4e}  xi(5000')={ff[-1]:+.4e}  -> {out}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
