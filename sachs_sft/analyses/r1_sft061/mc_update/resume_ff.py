"""Recover a stopped, single-angle FF run without repeating completed seeds.

The caller must first verify the old process has stopped. This helper checks
the original input and implementation hashes, archives the existing checkpoint,
and runs only the missing seed suffix. The completion sentinel is written last.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import shutil
import time

import numpy as np

from run_ff import sampled_covariance
from runtime import configure, sha256


def atomic_npz(path, **values):
    temporary = path.with_name(path.stem + ".recovery-tmp.npz")
    if temporary.exists():
        raise FileExistsError(temporary)
    np.savez(temporary, **values)
    temporary.replace(path)


def atomic_json(path, values):
    temporary = path.with_name(path.name + ".recovery-tmp")
    if temporary.exists():
        raise FileExistsError(temporary)
    temporary.write_text(json.dumps(values, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    folder = args.run.resolve(strict=True)
    manifest_path = folder / "manifest.json"
    seed_path = folder / "seeds_g00.npz"
    if (folder / "complete.json").exists() or (folder / "reference_g00.json").exists():
        raise ValueError("Recovery only supports a stopped run without final reference products")
    manifest = json.loads(manifest_path.read_text())
    settings = manifest["arguments"]
    if settings["reference_only"] or len(settings["gammas"]) != 1:
        raise ValueError("Recovery is restricted to one Monte Carlo angle")
    if manifest["runner_sha256"] != sha256(Path(__file__).with_name("run_ff.py")):
        raise ValueError("The original FF runner changed after this run started")
    for source, digest in manifest["source_sha256"].items():
        if sha256(Path(source)) != digest:
            raise ValueError(f"Original implementation hash changed: {source}")
    with np.load(seed_path) as product:
        gamma = float(product["gamma"])
        means = np.asarray(product["means"], dtype=float).tolist()
        errors = np.asarray(product["errors"], dtype=float).tolist()
        pairs = np.asarray(product["pairs"], dtype=int).tolist()
        if bool(product["complete"]):
            raise ValueError("Checkpoint is already marked complete")
    completed = len(means)
    if gamma != settings["gammas"][0] or not 0 < completed < settings["seeds"]:
        raise ValueError("Expected a nonempty incomplete prefix of the configured seed sample")
    if len(errors) != completed or len(pairs) != completed:
        raise ValueError("Checkpoint column lengths differ")
    if not np.all(np.isfinite(means)) or not np.all(np.isfinite(errors)) or min(errors) < 0:
        raise ValueError("Invalid checkpoint values")
    if pairs != [settings["n_real"] // 2] * completed:
        raise ValueError("Prior samples rejected realizations or used incomplete batches")
    runtime = configure(Path(manifest["table_path"]), Path(manifest["response_path"]),
                        manifest["lambda_source_mpc"])
    for key in ("table_sha256", "response_sha256", "covariance_sha256", "source_sha256"):
        if runtime.provenance[key] != manifest[key]:
            raise ValueError(f"Recovered runtime differs from original: {key}")
    recovery = {
        "original_manifest_sha256": sha256(manifest_path),
        "original_checkpoint_sha256": sha256(seed_path),
        "resume_helper_sha256": sha256(Path(__file__)),
        "completed_seed_count": completed,
        "missing_seed_indices": list(range(completed, settings["seeds"])),
        "missing_rng_seeds": [settings["seed_base"] + 991 * index
                              for index in range(completed, settings["seeds"])],
        "status": "validated, no computation performed" if args.check_only else "running",
    }
    print(json.dumps(recovery, indent=2), flush=True)
    if args.check_only:
        return
    archive = folder / f"seeds_g00.partial{completed}.npz"
    recovery_path = folder / "recovery.json"
    if archive.exists() or recovery_path.exists():
        raise FileExistsError("Existing recovery evidence must be inspected before another attempt")
    shutil.copy2(seed_path, archive)
    recovery_path.write_text(json.dumps(recovery, indent=2) + "\n")
    cfg = runtime.core.MCConfig(
        lam_source=manifest["lambda_source_mpc"], n_lambda=settings["n_lambda"],
        n_real=settings["n_real"], batch_size=settings["batch"],
        apply_anchor=False, seed=settings["seed_base"],
    )
    started = time.monotonic()
    for index in range(completed, settings["seeds"]):
        seeded = replace(cfg, seed=settings["seed_base"] + 991 * index)
        mean, error, count, blown = runtime.anti.ff_moment_anti(runtime.core, seeded, gamma)
        if (blown or count != settings["n_real"] // 2 or not np.isfinite(mean)
                or not np.isfinite(error) or error < 0):
            raise RuntimeError("Recovered sample rejected realizations or produced a nonfinite mean")
        means.append(mean)
        errors.append(error)
        pairs.append(count)
        atomic_npz(seed_path, gamma=gamma, means=means, errors=errors, pairs=pairs, complete=False)
        print(f"Recovered seed index {index}", flush=True)
    kernels, diagnostics = sampled_covariance(runtime, cfg, gamma)
    original_kernels = runtime.ff.kernels
    references = []
    try:
        runtime.ff.kernels = kernels
        for size in settings["reference_n"]:
            references.append(runtime.ff.ff_terms(gamma, n=size, kind="sampled"))
    finally:
        runtime.ff.kernels = original_kernels
    payload = dict(
        gamma=gamma, **diagnostics, reference_n=settings["reference_n"],
        references=[{key: value.tolist() for key, value in reference.items()} for reference in references],
        mc=float(np.mean(means)), mc_se=float(np.sqrt(np.sum(np.square(errors))) / len(means)),
        elapsed_seconds=time.monotonic() - started,
        elapsed_seconds_scope="Recovery continuation only, earlier seed durations unavailable",
    )
    atomic_json(folder / "reference_g00.json", payload)
    atomic_npz(seed_path, gamma=gamma, means=means, errors=errors, pairs=pairs, complete=True)
    recovery["status"] = "complete"
    recovery["completed_checkpoint_sha256"] = sha256(seed_path)
    atomic_json(recovery_path, recovery)
    atomic_json(folder / "complete.json", {"status": "complete", "reference_only": False})
    print(json.dumps({"mc": payload["mc"], "mc_se": payload["mc_se"],
                      "elapsed_seconds": payload["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
