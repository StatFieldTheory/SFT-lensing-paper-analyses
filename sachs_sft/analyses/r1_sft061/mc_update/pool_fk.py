"""Pool independent FK blocks using the historical linear sigma extrapolation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from runtime import DEFAULT_FOLD, fold_values, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", nargs="+", type=Path, required=True)
    parser.add_argument("--fold", type=Path, default=DEFAULT_FOLD)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.suffix != ".npz":
        parser.error("Output must have an explicit .npz suffix")
    if len(args.blocks) < 2 or len({path.resolve() for path in args.blocks}) != len(args.blocks):
        parser.error("At least two distinct independent blocks are required")
    if args.out.exists() or args.out.with_suffix(".json").exists():
        raise FileExistsError(args.out)
    manifests = []
    for path in args.blocks:
        complete = json.loads((path / "complete.json").read_text())
        if complete["status"] != "complete" or complete["reference_only"]:
            raise ValueError(f"Incomplete or reference-only block: {path}")
        manifests.append(json.loads((path / "manifest.json").read_text()))
    first = manifests[0]
    settings = first["arguments"]
    if settings["sigmas"] != [8.0, 4.0]:
        raise ValueError("This pooling implementation requires paired sigma=8 and sigma=4 rows")
    signature_keys = ["lambda_source_mpc", "table_sha256", "response_sha256", "covariance_sha256", "source_sha256", "runner_sha256"]
    argument_keys = ["gammas", "sigmas", "n_lambda", "n_real", "batch", "seeds"]
    seed_sets = []
    for manifest in manifests:
        if any(manifest[key] != first[key] for key in signature_keys):
            raise ValueError("Input or implementation provenance differs across blocks")
        if any(manifest["arguments"][key] != settings[key] for key in argument_keys):
            raise ValueError("Sampling settings differ across blocks")
        base = manifest["arguments"]["seed_base"]
        seeds = {base + 991 * index for index in range(settings["seeds"])}
        if any(seeds & previous for previous in seed_sets):
            raise ValueError("Seed streams overlap between purportedly independent blocks")
        seed_sets.append(seeds)
    gammas = np.asarray(settings["gammas"])
    analytic, provenance = fold_values(args.fold, gammas, first["lambda_source_mpc"],
                                        Path(first["table_path"]), Path(first["response_path"]))
    block_values = np.empty((len(args.blocks), len(gammas)))
    exact = np.empty_like(block_values)
    injected = np.empty_like(block_values)
    targets = np.empty_like(block_values)
    row_hashes = {}
    for block_index, folder in enumerate(args.blocks):
        for gamma_index, gamma in enumerate(gammas):
            rows = []
            for sigma_index, sigma in enumerate([8.0, 4.0]):
                path = folder / f"row_g{gamma_index:02d}_s{sigma_index:02d}.npz"
                with np.load(path) as data:
                    row = {key: np.array(data[key]) for key in data.files}
                if not row["complete"] or row["gamma"] != gamma or row["sigma"] != sigma:
                    raise ValueError(f"Incomplete or mismatched row: {path}")
                if row["seeds"].size != settings["seeds"] or not np.all(np.isfinite(row["seeds"])):
                    raise ValueError(f"Missing or nonfinite seed estimates: {path}")
                rows.append(row)
                row_hashes[str(path.resolve())] = sha256(path)
            block_values[block_index, gamma_index] = np.mean(2 * rows[1]["seeds"] - rows[0]["seeds"])
            exact[block_index, gamma_index] = 2 * rows[1]["exact"] - rows[0]["exact"]
            injected[block_index, gamma_index] = rows[1]["ref_inj"]
            targets[block_index, gamma_index] = rows[1]["ref_true"]
    for values in (exact, injected, targets):
        np.testing.assert_allclose(values, np.broadcast_to(values[0], values.shape), rtol=1e-13, atol=0)
    mean = block_values.mean(axis=0)
    error = block_values.std(axis=0, ddof=1) / np.sqrt(len(args.blocks))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, gamma=gammas, fk=mean, err=error, analytic=analytic,
             blocks=block_values, extrapolated_expectation=exact[0],
             ref_inj=injected[0], ref_true=targets[0],
             n_blocks=len(args.blocks), seeds_per_block=settings["seeds"])
    provenance.update({
        "input_manifest_paths": [str((path / "manifest.json").resolve()) for path in args.blocks],
        "input_manifest_sha256": {str((path / "manifest.json").resolve()): sha256(path / "manifest.json") for path in args.blocks},
        "row_sha256": row_hashes, "runner_sha256": sha256(Path(__file__)),
        "error_scope": "Standard error from independent block means after paired linear sigma extrapolation",
        "systematics": "Finite-lambda-step and finite-sigma errors are recorded separately from sampling error. The expectation uses population centering, unlike finite-batch MC.",
        "ratio": (mean / analytic).tolist(),
        "ratio_error": (error / np.abs(analytic)).tolist(),
        "expectation_over_local_injected": (exact[0] / injected[0]).tolist(),
        "injection_over_target": (injected[0] / targets[0]).tolist(),
        "local_target_over_fold": (targets[0] / analytic).tolist(),
    })
    args.out.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({key: provenance[key] for key in ("ratio", "ratio_error", "expectation_over_local_injected", "injection_over_target", "local_target_over_fold")}, indent=2))


if __name__ == "__main__":
    main()
