"""Selected-point higher-F-order diagnostic using the existing FF sampler.

Only the coefficient multiplying the existing nonlinear vertex is changed.
The same Gaussian draws are used as in a completed baseline run. The existing
Order-F-squared reference motivates division by the squared coupling. A finite
coupling comparison measures coupling dependence, not an exact F-squared limit.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from runtime import configure, sha256


def scaled_core(core, strength: float):
    """Minimal sampler interface, without mutating the imported core module."""
    def vertex(state):
        return strength * core._f_vertex(state)

    return SimpleNamespace(_precompute=core._precompute, _f_vertex=vertex)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--strength", type=float, default=0.5)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.strength) or not 0 < args.strength < 1:
        parser.error("Require a finite coupling strength strictly between zero and one")
    complete = json.loads((args.baseline / "complete.json").read_text())
    if complete["status"] != "complete" or complete["reference_only"]:
        raise ValueError("A completed Monte Carlo baseline is required")
    baseline_path = args.baseline / "manifest.json"
    baseline = json.loads(baseline_path.read_text())
    settings = baseline["arguments"]
    if args.gamma not in settings["gammas"] or settings["seeds"] < 2:
        raise ValueError("Require the selected angle and at least two baseline seeds")
    gamma_index = settings["gammas"].index(args.gamma)
    seed_path = args.baseline / f"seeds_g{gamma_index:02d}.npz"
    with np.load(seed_path) as product:
        baseline_means = np.asarray(product["means"], dtype=float)
        if not product["complete"] or len(baseline_means) != settings["seeds"]:
            raise ValueError("Incomplete baseline seed sample")
        if not np.all(np.isfinite(baseline_means)):
            raise ValueError("Nonfinite baseline seed sample")
        if not np.all(product["pairs"] == settings["n_real"] // 2):
            raise ValueError("Baseline rejected realizations or used incomplete batches")
    runtime = configure(Path(baseline["table_path"]), Path(baseline["response_path"]),
                        baseline["lambda_source_mpc"])
    for key in ("table_sha256", "response_sha256", "covariance_sha256", "source_sha256"):
        if runtime.provenance[key] != baseline[key]:
            raise ValueError(f"Baseline and coupling probe provenance differ: {key}")
    args.out.mkdir(parents=True, exist_ok=False)
    manifest = dict(runtime.provenance)
    manifest.update({
        "baseline_manifest": str(baseline_path.resolve()),
        "baseline_manifest_sha256": sha256(baseline_path),
        "baseline_seeds": str(seed_path.resolve()),
        "baseline_seeds_sha256": sha256(seed_path),
        "runner_sha256": sha256(Path(__file__)),
        "arguments": {"gamma": args.gamma, "strength": args.strength,
                      "baseline": str(args.baseline.resolve())},
        "sampling_settings": settings,
        "scope": "Same sampler and random seed streams with a scaled existing F vertex. Difference of normalized finite-coupling estimates, not an exact perturbative-order decomposition.",
    })
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    cfg = runtime.core.MCConfig(
        lam_source=baseline["lambda_source_mpc"], n_lambda=settings["n_lambda"],
        n_real=settings["n_real"], batch_size=settings["batch"],
        apply_anchor=False, seed=settings["seed_base"],
    )
    sampler = scaled_core(runtime.core, args.strength)
    means, errors = [], []
    started = time.monotonic()
    for index in range(settings["seeds"]):
        seeded = replace(cfg, seed=settings["seed_base"] + 991 * index)
        mean, error, count, blown = runtime.anti.ff_moment_anti(sampler, seeded, args.gamma)
        if blown or count != settings["n_real"] // 2 or not np.isfinite(mean):
            raise RuntimeError("Nonfinite estimate or rejected realizations")
        means.append(mean / args.strength**2)
        errors.append(error / args.strength**2)
        np.savez(args.out / "seeds.npz", normalized_means=means,
                 normalized_errors=errors, baseline_means=baseline_means,
                 strength=args.strength, complete=False)
    normalized = np.asarray(means)
    differences = baseline_means - normalized
    baseline_mean = float(baseline_means.mean())
    difference = float(differences.mean())
    difference_se = float(differences.std(ddof=1) / np.sqrt(len(differences)))
    np.savez(args.out / "seeds.npz", normalized_means=means,
             normalized_errors=errors, baseline_means=baseline_means,
             paired_differences=differences, strength=args.strength, complete=True)
    summary = {
        "gamma": args.gamma, "strength": args.strength,
        "baseline_mean": baseline_mean, "normalized_reduced_coupling_mean": float(normalized.mean()),
        "paired_mean_difference": difference, "paired_difference_se": difference_se,
        "difference_over_baseline": difference / baseline_mean,
        "difference_se_over_baseline": difference_se / abs(baseline_mean),
        "elapsed_seconds": time.monotonic() - started,
        "scope": "Selected-angle finite-coupling sensitivity. Both estimates include higher F orders and use the same finite Euler step.",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.out / "complete.json").write_text(json.dumps({"status": "complete"}) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
