"""Independent local FK contraction using existing discrete-grid machinery."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from runtime import (DEFAULT_FOLD, DEFAULT_RESPONSE, DEFAULT_TABLE, SOURCE,
                     check_response, configure, fold_values, sha256)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--response", type=Path, default=DEFAULT_RESPONSE)
    parser.add_argument("--fold", type=Path, default=DEFAULT_FOLD)
    parser.add_argument("--allow-missing-fold", action="store_true")
    parser.add_argument("--lambda-source", type=float, default=SOURCE)
    parser.add_argument("--gammas", type=float, nargs="+", default=[1.0, 2.6, 6.7, 17.3])
    parser.add_argument("--n-lambda", type=int, nargs="+", default=[500, 1000, 2000, 4000])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if min(args.n_lambda) < 2:
        parser.error("Each grid must have at least two nodes")
    if args.out.exists():
        raise FileExistsError(args.out)
    runtime = configure(args.table, args.response, args.lambda_source)
    result = dict(runtime.provenance)
    result["runner_sha256"] = sha256(Path(__file__))
    if args.fold.is_file():
        values, provenance = fold_values(args.fold, args.gammas, args.lambda_source,
                                         args.table, args.response)
        result.update(provenance)
        result["fold_values"] = values.tolist()
    elif args.allow_missing_fold:
        result["fold_status"] = "Unavailable, no fold comparison established"
    else:
        raise FileNotFoundError(args.fold)
    result["gamma"] = args.gammas
    result["n_lambda"] = args.n_lambda
    result["scope"] = "Existing fk_expect_exact.fk_reference, target cumulant contraction on discrete source grids. Grid refinement checks the continuum limit separately from finite-correlation MC."
    rows = []
    for size in args.n_lambda:
        started = time.monotonic()
        grid = runtime.ex.build_grid(n_lambda=size, lam_source=args.lambda_source, apply_anchor=False)
        check_response(runtime, grid.lam, grid.resp)
        row = []
        for gamma in args.gammas:
            cosine = float(np.cos(np.deg2rad(gamma / 60)))
            target = np.array([runtime.ex._DS.zeta6(cosine, float(lam)) for lam in grid.lam])
            row.append(float(runtime.ex.fk_reference(grid, target)[0, 3]))
        rows.append(row)
        print(json.dumps({"n_lambda": size, "local_fk": row,
                          "elapsed_seconds": time.monotonic() - started}), flush=True)
    result["local_fk"] = rows
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
