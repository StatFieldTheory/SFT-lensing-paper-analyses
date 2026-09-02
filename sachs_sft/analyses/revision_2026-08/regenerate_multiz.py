"""Regenerate the multi-redshift figure's data with the corrected vertex.

Only the FK slice needs recomputing. Order-0 carries no three-point vertex,
and FF is the nonlinear-propagation term, built on the F tensor rather than
on the three-point vertex K, so neither is touched by the corrections. This
driver therefore materialises one corrected config for FK, redirects the
paper's own config map at it, and calls the paper's own compute script, so
the redshift grid, the sweep geometry and the output schema are unchanged.

Usage::

    <sft-wick python> regenerate_multiz.py \
        --table ../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/table_permclosed.npz \
        --callable ../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[3]
_A3 = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "analysis3")
_TALK = _A3  # multiz_sweep.py lives beside the compute script since 2026-09-02


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--callable", dest="callable_src", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path,
                    default=_HERE / "outputs" / "_variant_multiz_conv")
    ap.add_argument("--reuse", type=Path, default=None,
                    help="an existing multi-z npz whose o0 and ff slices are "
                         "reused instead of recomputed")
    args = ap.parse_args()

    sys.path.insert(0, str(_TALK))
    sys.path.insert(0, str(_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables" / "kappa3_vertex" / "rebuild"))
    import run_fk_variant as rfv  # noqa: E402

    # One corrected FK config, written where every other variant goes so the
    # deployed run directory is never touched.
    _, config = rfv.materialise_variant(
        args.table.resolve(), args.work_dir.resolve(),
        n_jobs=6, callable_src=args.callable_src.resolve())
    print(f"[multiz] corrected FK config: {config}")

    compute = _load(_A3 / "compute_multiz_kappa_2pcf.py", "compute_multiz")
    compute.R.CONFIGS = dict(compute.R.CONFIGS)
    compute.R.CONFIGS["fk"] = config
    print(f"[multiz] redshifts: {compute.Z.tolist()}")

    if args.reuse is not None:
        # The paper's compute script sweeps all three channels.  Only FK reads
        # the three-point vertex, so Order-0 and FF are served from the existing
        # result and the sweep runs once.  Besides the hours saved this keeps
        # the two untouched channels bit-identical to the published figure, so
        # any difference in the new figure is attributable to the correction.
        import numpy as np
        cached = np.load(args.reuse, allow_pickle=True)
        gamma_cached = np.asarray(cached["gamma"], float)
        real_sweep = compute.R.run_vertex_sweep

        def sweep(cfg, lam, label):
            if label == "fk":
                return real_sweep(cfg, lam, label)
            print(f"[{label}] reused from {args.reuse.name} "
                  f"(no three-point vertex; unaffected by the correction)",
                  flush=True)
            return gamma_cached, np.asarray(cached[label], float)

        compute.R.run_vertex_sweep = sweep
        if not np.allclose(np.asarray(cached["z"], float), compute.Z):
            raise SystemExit(f"redshift grid of {args.reuse} does not match "
                             f"the compute script's: {cached['z']} vs {compute.Z}")
    return compute.main()


if __name__ == "__main__":
    raise SystemExit(main())
