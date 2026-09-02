"""Regenerate the Monte-Carlo validation figure with the corrected vertex.

Both sides of that figure read the same three-point vertex: the analytic FK
curve through the folded sweep, and the Monte-Carlo through
``driver_stats``, which imports the vertex callable directly. Correcting
only one side would manufacture a disagreement that is not there, so this
driver redirects both:

* ``fig_xi_channels._load_kk`` is wrapped so the FK run reads the corrected
  folded result instead of the deployed one;
* ``driver_stats._k3`` is replaced by the permutation-aware callable bound
  to the corrected table, so the Monte-Carlo is seeded with the same
  statistics.

Everything else, including the realisation count, the colored-noise
correlation length and the plotting, is the paper's own.

Usage::

    <PyCCL python> regenerate_mc.py \
        --table ../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/table_permclosed.npz \
        --n-real 24000
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[4]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--fk-xi", type=Path, required=True,
                    help="the corrected folded FK sweep")
    ap.add_argument("--n-real", type=int, default=24_000)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--ff-real", type=int, default=48_000)
    ap.add_argument("--ff-nl", type=int, default=4000)
    ap.add_argument("--gammas", type=float, nargs="+", default=None,
                    help="override the generator's GAMMA_MC. The Monte-Carlo "
                         "has discriminating power only at small separation: "
                         "its AR(1) regulator sigma_lambda subtends "
                         "sigma_lambda/chi ~ 14', and the FK amplitude is "
                         "proportional to sigma_lambda, so beyond a few "
                         "arcminutes the markers measure the regulator, not "
                         "the diagram. Restricting the marker list keeps the "
                         "figure honest about where the check applies.")
    args = ap.parse_args()

    sys.path.insert(0, str(_MC))
    sys.path.insert(0, str(_HERE.parent / "callable_fixed"))

    # 1. the Monte-Carlo side
    import driver_stats as ds  # noqa: E402
    import perm_aware_kappa3_callable as pa  # noqa: E402

    pa.TABLE_PATH = args.table.resolve()
    pa._CACHE = None
    ds._k3 = pa
    print(f"[mc] Monte-Carlo vertex -> {pa.TABLE_PATH.name} "
          f"via the permutation-aware callable")

    # 2. the analytic side
    fig = _load(_MC / "fig_xi_channels.py", "fig_xi_channels")
    original = fig._load_kk

    def patched(name: str):
        if name != "C_corr_op_K_limber_FK":
            return original(name)
        d = np.load(args.fk_xi, allow_pickle=True)
        m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
        import math
        gammas = []
        for xi, yi in zip(d["x"][m], d["y"][m]):
            xi = np.asarray(xi, float); yi = np.asarray(yi, float)
            gammas.append(math.degrees(math.acos(float(np.clip(
                np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)),
                -1, 1)))) * 60)
        order = np.argsort(gammas)
        return np.asarray(gammas)[order], np.asarray(d["value"], float)[m][order]

    fig._load_kk = patched
    if args.gammas is not None:
        fig.GAMMA_MC = tuple(args.gammas)
        print(f"[mc] GAMMA_MC restricted to {fig.GAMMA_MC}")
    print(f"[mc] analytic FK -> {args.fk_xi.name}")
    print(f"[mc] {args.n_real} realisations, n_lambda {args.n_lambda}, "
          f"sigma_lambda {args.sigma}")
    fig.main(args.n_real, args.n_lambda, args.sigma, args.seeds,
             args.ff_real, args.ff_nl, False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
