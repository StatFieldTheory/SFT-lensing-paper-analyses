"""FF Monte Carlo and continuum contractions using its actual local covariance.

The state covariance is reconstructed from the sampler's L_white factors.
The FF contraction is the existing ff_functional.ff_terms implementation.
It is a continuum evaluation with interpolated finite-step covariance, not an
exact Euler prediction. Compare both --n-lambda and --reference-n values.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from runtime import (DEFAULT_RESPONSE, DEFAULT_TABLE, SOURCE, check_response,
                     configure, sha256)


def sampled_covariance(runtime, cfg, gamma):
    cosine = float(np.cos(np.deg2rad(gamma / 60)))
    bg, lam, dlam, resp, factors, _, _, _ = runtime.core._precompute(cfg, cosine)
    check_response(runtime, lam, resp)
    D = np.asarray(bg.D(lam))
    noise_cov = factors @ np.swapaxes(factors, 1, 2)
    state = np.empty_like(noise_cov)
    accumulated = np.zeros((6, 6))
    for index in range(len(lam)):
        accumulated = resp[index] ** 2 * accumulated + noise_cov[index]
        state[index] = accumulated
    latent = D[:, None, None] ** 4 * state
    # Independent telescoping reconstruction of the same sampled covariance.
    cumulative = np.cumsum(D[:, None, None] ** 4 * noise_cov, axis=0)
    # Cross components can cross zero. Scale each component by its own peak
    # along lambda so cancellation does not turn roundoff into a false failure.
    component_error = np.max(np.abs(latent - cumulative), axis=0)
    component_peak = np.max(np.abs(cumulative), axis=0)
    component_tolerance = 2e-12 * component_peak + 1e-16
    if (not np.all(np.isfinite(component_error))
            or not np.all(np.isfinite(component_tolerance))
            or np.any(component_error > component_tolerance)):
        raise AssertionError("Sampled covariance reconstructions exceed componentwise peak-scaled tolerance")
    component_tolerance_ratio = float(np.max(component_error / component_tolerance))

    def kernels(_gamma, requested, requested_D, kind):
        if kind != "sampled" or _gamma != gamma:
            raise ValueError("This closure is specific to the requested sampled covariance")
        if requested[0] < lam[0] or requested[-1] > lam[-1]:
            raise ValueError("Cannot extrapolate sampled covariance")
        interpolated = np.empty((len(requested), 6, 6))
        for a in range(6):
            for b in range(6):
                interpolated[:, a, b] = np.interp(requested, lam, latent[:, a, b])
        minimum = np.minimum(np.arange(len(requested))[:, None], np.arange(len(requested))[None, :])
        prefactor = 1 / (requested_D[:, None] ** 2 * requested_D[None, :] ** 2)
        covariance = prefactor[:, :, None, None] * interpolated[minimum]
        return covariance[:, :, :3, 3:], covariance[:, :, :3, :3], covariance[:, :, 3:, 3:]

    # Reuse r6_exact_ff_mean.py's preceding-node convention independently.
    weights = np.full(len(lam), dlam)
    weights[[0, -1]] *= 0.5
    window = D**2 * np.cumsum((weights / D**2)[::-1])[::-1]
    previous = np.zeros_like(state)
    previous[1:] = state[:-1]
    means = np.empty(6)
    for offset in (0, 3):
        means[offset:offset + 3] = np.einsum(
            "j,abc,jbc->a", dlam * window, runtime.ff.FS,
            previous[:, offset:offset + 3, offset:offset + 3], optimize=True,
        )
    diagnostics = {
        "euler_order_F_mean_ray1": means[:3].tolist(),
        "euler_order_F_mean_ray2": means[3:].tolist(),
        "euler_disconnected_kk": float(means[0] * means[3]),
        "sampled_noise_sha256": hashlib.sha256(noise_cov.tobytes()).hexdigest(),
        "state_covariance_recursion_max_relative": float(np.max(np.abs(latent - cumulative)) / np.max(np.abs(cumulative))),
        "state_covariance_recursion_max_component_tolerance_ratio": component_tolerance_ratio,
    }
    return kernels, diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--response", type=Path, default=DEFAULT_RESPONSE)
    parser.add_argument("--lambda-source", type=float, default=SOURCE)
    parser.add_argument("--gammas", type=float, nargs="+", default=[1.0, 2.6, 6.7, 17.3])
    parser.add_argument("--n-lambda", type=int, default=4000)
    parser.add_argument("--reference-n", type=int, nargs="+", default=[241, 481])
    parser.add_argument("--n-real", type=int, default=96000)
    parser.add_argument("--batch", type=int, default=3000)
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--seed-base", type=int, default=7010001)
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if min(args.n_real, args.batch, args.seeds) < 1 or args.n_lambda < 2 or min(args.reference_n) < 3:
        parser.error("Require positive MC sizes and at least three reference nodes")
    if args.batch % 2 or args.n_real % args.batch:
        parser.error("Require an even batch size and complete batches")
    if len(set(args.gammas)) != len(args.gammas) or len(set(args.reference_n)) != len(args.reference_n):
        parser.error("Require distinct angles and reference grid sizes")
    args.out.mkdir(parents=True, exist_ok=False)
    runtime = configure(args.table, args.response, args.lambda_source)
    manifest = dict(runtime.provenance)
    manifest["arguments"] = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    manifest["runner_sha256"] = sha256(Path(__file__))
    manifest["reference_scope"] = "Continuum Order-F-squared contraction with the actual finite-step local sampled covariance. Not the nonlocal cosmological FF and not an exact Euler reference."
    manifest["mc_scope"] = "Antithetic Gaussian F-on minus F-off moment, including higher F orders"
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    cfg = runtime.core.MCConfig(lam_source=args.lambda_source, n_lambda=args.n_lambda,
                               n_real=args.n_real, batch_size=args.batch,
                               apply_anchor=False, seed=args.seed_base)
    original_kernels = runtime.ff.kernels
    for gamma_index, gamma in enumerate(args.gammas):
        started = time.monotonic()
        kernels, diagnostics = sampled_covariance(runtime, cfg, gamma)
        references = []
        try:
            runtime.ff.kernels = kernels
            for size in args.reference_n:
                terms = runtime.ff.ff_terms(gamma, n=size, kind="sampled")
                references.append(terms)
        finally:
            runtime.ff.kernels = original_kernels
        means, errors, pairs = [], [], []
        for seed_index in range(0 if args.reference_only else args.seeds):
            seeded_cfg = replace(cfg, seed=args.seed_base + 991 * seed_index)
            mean, error, count, blown = runtime.anti.ff_moment_anti(runtime.core, seeded_cfg, gamma)
            if blown or count != args.n_real // 2 or not np.isfinite(mean):
                raise RuntimeError("Nonfinite FF estimate or rejected realizations")
            means.append(mean)
            errors.append(error)
            pairs.append(count)
            np.savez(args.out / f"seeds_g{gamma_index:02d}.npz", gamma=gamma,
                     means=means, errors=errors, pairs=pairs, complete=False)
        np.savez(args.out / f"seeds_g{gamma_index:02d}.npz", gamma=gamma,
                 means=means, errors=errors, pairs=pairs, complete=True)
        payload = dict(gamma=gamma, **diagnostics,
                       reference_n=args.reference_n,
                       references=[{key: value.tolist() for key, value in ref.items()} for ref in references],
                       mc=float(np.mean(means)) if means else None,
                       mc_se=float(np.sqrt(np.sum(np.square(errors))) / len(means)) if means else None,
                       elapsed_seconds=time.monotonic() - started)
        (args.out / f"reference_g{gamma_index:02d}.json").write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps({key: value for key, value in payload.items() if key != "references"}), flush=True)
    (args.out / "complete.json").write_text(json.dumps({"status": "complete", "reference_only": args.reference_only}) + "\n")


if __name__ == "__main__":
    main()
