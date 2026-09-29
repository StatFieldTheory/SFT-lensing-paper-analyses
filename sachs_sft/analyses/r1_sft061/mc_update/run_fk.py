"""One fresh-process, checkpointed corrected-input FK Monte Carlo block."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from runtime import (DEFAULT_FOLD, DEFAULT_RESPONSE, DEFAULT_TABLE, SOURCE,
                     check_response, configure, fold_values, sha256)


def sampled_nodes(runtime, grid, cosine, sigma):
    """Calibrate against the covariance actually sampled after PSD clipping."""
    ex = runtime.ex
    V, _, _, rho = ex.node_stats(grid, cosine, sigma, calibrate="smeared")
    factors = np.array([runtime.core._cholesky_psd(value) for value in V])
    sampled = factors @ np.swapaxes(factors, 1, 2)
    A = np.empty_like(sampled)
    A[0] = sampled[0]
    for index in range(1, len(A)):
        A[index] = rho**2 * A[index - 1] + (1 - rho**2) * sampled[index]
    matrix = ex.smeared_covariance(grid, A, rho)
    target = np.array([ex._DS.zeta6(cosine, float(lam)) for lam in grid.lam])
    Q = np.array([ex._DS.solve_Q(mat, cumulant) for mat, cumulant in zip(matrix, target)])
    injected = ex.zeta_injected(matrix, Q)
    scale = max(float(np.max(np.abs(V))), np.finfo(float).tiny)
    eigenvalues = np.linalg.eigvalsh(V)
    stats = {
        "covariance_psd_max_abs_change_over_peak": float(np.max(np.abs(sampled - V)) / scale),
        "covariance_negative_node_count": int(np.count_nonzero(eigenvalues[:, 0] < -1e-14 * scale)),
        "covariance_min_eigenvalue_over_peak": float(eigenvalues.min() / scale),
        "injection_max_abs_error_over_peak": float(np.max(np.abs(injected - target)) / max(np.max(np.abs(target)), np.finfo(float).tiny)),
    }
    # Keep original V so the sampler constructs exactly the factors used above.
    return (V, A, Q, rho), target, injected, stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--response", type=Path, default=DEFAULT_RESPONSE)
    parser.add_argument("--fold", type=Path, default=DEFAULT_FOLD)
    parser.add_argument("--allow-missing-fold", action="store_true", help="Probe only, records missing comparison")
    parser.add_argument("--lambda-source", type=float, default=SOURCE)
    parser.add_argument("--gammas", type=float, nargs="+", default=[1.0, 2.6, 6.7, 17.3])
    parser.add_argument("--sigmas", type=float, nargs="+", default=[8.0, 4.0])
    parser.add_argument("--n-lambda", type=int, default=1000)
    parser.add_argument("--n-real", type=int, default=24000)
    parser.add_argument("--batch", type=int, default=2000)
    parser.add_argument("--seeds", type=int, default=24)
    parser.add_argument("--seed-base", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-only", action="store_true")
    args = parser.parse_args()
    if args.n_lambda < 2 or min(args.n_real, args.batch, args.seeds) < 1 or args.n_real % args.batch:
        parser.error("Require n_lambda >= 2 and positive complete MC batches")
    if min(args.sigmas) <= 0 or len(set(args.gammas)) != len(args.gammas) or len(set(args.sigmas)) != len(args.sigmas):
        parser.error("Require positive sigma and distinct gamma/sigma values")
    args.out.mkdir(parents=True, exist_ok=False)
    runtime = configure(args.table, args.response, args.lambda_source)
    manifest = dict(runtime.provenance)
    manifest["arguments"] = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    manifest["runner_sha256"] = sha256(Path(__file__))
    if args.fold.is_file():
        analytic, fold_info = fold_values(args.fold, args.gammas, args.lambda_source,
                                          args.table, args.response)
        manifest.update(fold_info)
    elif args.allow_missing_fold:
        analytic = np.full(len(args.gammas), np.nan)
        manifest["fold_status"] = "Missing at probe time, no paper-fold comparison established"
    else:
        raise FileNotFoundError(args.fold)
    grid = runtime.ex.build_grid(n_lambda=args.n_lambda, lam_source=args.lambda_source, apply_anchor=False)
    manifest["response_max_abs_error"] = check_response(runtime, grid.lam, grid.resp)
    manifest["expectation_scope"] = "Discrete polynomial expectation with population centering. MC uses finite-batch sample centering and covariance normalization, which can add O(1/batch) bias."
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for gamma_index, gamma in enumerate(args.gammas):
        cosine = float(np.cos(np.deg2rad(gamma / 60)))
        for sigma_index, sigma in enumerate(args.sigmas):
            started = time.monotonic()
            nodes, target, injected, diagnostics = sampled_nodes(runtime, grid, cosine, sigma)
            first, second = runtime.ex.expectation(grid, cosine, sigma, V=nodes[0], A=nodes[1], Q=nodes[2], rho=nodes[3])
            exact = float((first + second + (first + second).T)[0, 3])
            row = dict(gamma=gamma, sigma=sigma, exact=exact,
                       ref_true=float(runtime.ex.fk_reference(grid, target)[0, 3]),
                       ref_inj=float(runtime.ex.fk_reference(grid, injected)[0, 3]),
                       analytic=float(analytic[gamma_index]), **diagnostics)
            seed_values, full_values, retained = [], [], []
            row_path = args.out / f"row_g{gamma_index:02d}_s{sigma_index:02d}.npz"
            for seed_index in range(0 if args.reference_only else args.seeds):
                result = runtime.fc.simulate_fk_complete(
                    gamma, sigma_lambda=sigma, n_lambda=args.n_lambda,
                    n_real=args.n_real, batch_size=args.batch,
                    seed=args.seed_base + 991 * seed_index, grid=grid, nodes=nodes,
                )
                if not np.isfinite(result.kk) or result.n_good != args.n_real:
                    raise RuntimeError("Nonfinite FK estimate or rejected realizations")
                seed_values.append(result.kk)
                full_values.append(float(result.fk_full[0, 0]))
                retained.append(result.n_good)
                np.savez(row_path, **row, seeds=seed_values, full_seeds=full_values,
                         retained=retained, complete=False)
            elapsed = time.monotonic() - started
            np.savez(row_path, **row, seeds=seed_values, full_seeds=full_values,
                     retained=retained, complete=True, elapsed_seconds=elapsed)
            print(json.dumps(dict(row, mc=float(np.mean(seed_values)) if seed_values else None,
                                  elapsed_seconds=elapsed)), flush=True)
    (args.out / "complete.json").write_text(json.dumps({"status": "complete", "reference_only": args.reference_only}) + "\n")


if __name__ == "__main__":
    main()
