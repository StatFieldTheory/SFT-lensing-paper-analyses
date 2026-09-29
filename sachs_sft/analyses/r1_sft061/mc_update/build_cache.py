"""Assemble the figure-6 cache from explicit, provenance-matched new products."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from runtime import DEFAULT_FOLD, SOURCE, fold_values, sha256

ROUNDOFF_OLD_RUNNER = "b43c7cf59e3aeb0930d600bf239548e79d9309d43317d8d98c9121a9bc6cd953"
ROUNDOFF_NEW_RUNNER = "301c56173e434ac41af189782f007945529570bad3dd5af356237a0b83ac0efe"
ROUNDOFF_AUDIT_SHA256 = "f18468d850704a133e8b49b5985226d4d6bf1605a2f04079176145cbbb33a8f6"


def validate_ff_runner_compatibility(reference: dict, mc: dict) -> dict | None:
    """Allow only the explicitly audited diagnostic-only runner transition."""
    if reference["runner_sha256"] == mc["runner_sha256"]:
        return None
    if (reference["runner_sha256"] != ROUNDOFF_NEW_RUNNER
            or mc["runner_sha256"] != ROUNDOFF_OLD_RUNNER):
        raise ValueError("FF runner mismatch is not the approved roundoff-diagnostic transition")
    audit_path = Path(__file__).with_name("ff_covariance_roundoff_audit.json")
    if sha256(audit_path) != ROUNDOFF_AUDIT_SHA256:
        raise ValueError("FF roundoff audit differs from the approved immutable evidence")
    audit = json.loads(audit_path.read_text())
    if (audit["status"] != "passed" or audit["old_runner_sha256"] != ROUNDOFF_OLD_RUNNER
            or audit["new_runner_sha256"] != ROUNDOFF_NEW_RUNNER):
        raise ValueError("FF roundoff audit does not approve these runner versions")
    for manifest in (reference, mc):
        if any(manifest[key] != value for key, value in audit["physical_provenance"].items()):
            raise ValueError("FF physical provenance differs from the audited calculation")
        if manifest["arguments"]["n_lambda"] != audit["n_lambda"]:
            raise ValueError("FF covariance grid differs from the audited calculation")
    if reference["arguments"]["reference_n"] != audit["reference_n"]:
        raise ValueError("FF contraction grids differ from the audited calculation")
    requested = reference["arguments"]["gammas"]
    audited = [row["gamma"] for row in audit["rows"]]
    if requested != audited:
        raise ValueError("FF angular grid differs from the audited calculation")
    return audit


def completed_manifest(folder: Path, allow_reference_only: bool) -> dict:
    complete = json.loads((folder / "complete.json").read_text())
    if complete["status"] != "complete" or (complete["reference_only"] and not allow_reference_only):
        raise ValueError(f"Required completed Monte Carlo product missing: {folder}")
    return json.loads((folder / "manifest.json").read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fk-fold", type=Path, default=DEFAULT_FOLD)
    parser.add_argument("--order0", type=Path,
                        default=Path(__file__).resolve().parents[1] / "figure_products" / "order0.npz")
    parser.add_argument("--fk-markers", type=Path, required=True)
    parser.add_argument("--ff-reference-dir", type=Path, required=True)
    parser.add_argument("--ff-mc-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.suffix != ".npz":
        parser.error("Output must have an explicit .npz suffix")
    if args.out.exists() or args.out.with_suffix(".json").exists():
        raise FileExistsError(args.out)
    ref_manifest = completed_manifest(args.ff_reference_dir, True)
    mc_manifest = completed_manifest(args.ff_mc_dir, False)
    keys = ["lambda_source_mpc", "response_sha256", "covariance_sha256", "source_sha256"]
    if any(ref_manifest[key] != mc_manifest[key] for key in keys):
        raise ValueError("FF reference and MC provenance differ")
    runner_audit = validate_ff_runner_compatibility(ref_manifest, mc_manifest)
    if ref_manifest["arguments"]["n_lambda"] != mc_manifest["arguments"]["n_lambda"]:
        raise ValueError("FF reference must use the same finite-step sampled covariance as MC")
    if ref_manifest["lambda_source_mpc"] != SOURCE:
        raise ValueError("FF product has the wrong source plane")
    provenance = {"source_redshift": 5.0, "lambda_source_mpc": SOURCE,
                  "ff_scope": "Continuum Order-F-squared reference for actual local sampled covariance",
                  "ff_mc_scope": "F-on minus F-off local Gaussian Monte Carlo, including higher F orders",
                  "runner_sha256": sha256(Path(__file__)), "input_sha256": {}}

    def record(path):
        provenance["input_sha256"][str(path.resolve())] = sha256(path)

    if runner_audit is not None:
        audit_path = Path(__file__).with_name("ff_covariance_roundoff_audit.json")
        record(audit_path)
        provenance["ff_runner_compatibility"] = {
            "scope": "Pinned diagnostic-only roundoff transition, physical provenance unchanged",
            "audit_path": str(audit_path.resolve()), "audit_sha256": ROUNDOFF_AUDIT_SHA256,
            "reference_runner_sha256": ROUNDOFF_NEW_RUNNER,
            "mc_runner_sha256": ROUNDOFF_OLD_RUNNER,
        }

    for path in [args.order0, args.ff_reference_dir / "manifest.json", args.ff_mc_dir / "manifest.json", args.fk_markers,
                 args.fk_markers.with_suffix(".json")]:
        record(path)
    reference_rows = []
    for index, gamma in enumerate(ref_manifest["arguments"]["gammas"]):
        path = args.ff_reference_dir / f"reference_g{index:02d}.json"
        record(path)
        row = json.loads(path.read_text())
        if row["gamma"] != gamma:
            raise ValueError("FF reference row has a mismatched angle")
        if (runner_audit is not None
                and row["sampled_noise_sha256"] != runner_audit["rows"][index]["sampled_noise_sha256"]):
            raise ValueError("FF reference sampled covariance differs from the audited calculation")
        chosen = int(np.argmax(row["reference_n"]))
        reference_rows.append((gamma, row["references"][chosen]["moment"][0][0]))
    reference_rows = np.asarray(sorted(reference_rows))
    mc_rows = []
    for index, gamma in enumerate(mc_manifest["arguments"]["gammas"]):
        path = args.ff_mc_dir / f"reference_g{index:02d}.json"
        record(path)
        row = json.loads(path.read_text())
        if row["gamma"] != gamma or row["mc"] is None:
            raise ValueError("FF Monte Carlo row is absent or mismatched")
        seed_path = args.ff_mc_dir / f"seeds_g{index:02d}.npz"
        record(seed_path)
        with np.load(seed_path) as data:
            if not data["complete"] or len(data["means"]) != mc_manifest["arguments"]["seeds"]:
                raise ValueError("FF Monte Carlo row lacks its completed seed sample")
        mc_rows.append((gamma, row["mc"], row["mc_se"]))
    mc_rows = np.asarray(sorted(mc_rows))
    if mc_rows[0, 0] < reference_rows[0, 0] or mc_rows[-1, 0] > reference_rows[-1, 0]:
        raise ValueError("FF reference does not cover all marker angles")

    with np.load(args.order0, allow_pickle=True) as data:
        mask = (data["a"] == 0) & (data["b"] == 0) & (data["order"] == 0)
        if not np.any(mask) or not np.allclose(data["t_final"][mask], SOURCE, atol=1e-10, rtol=0):
            raise ValueError("Order-0 guide product must contain matching source Order-0 rows")
        x = np.asarray(list(data["x"][mask]), dtype=float)
        y = np.asarray(list(data["y"][mask]), dtype=float)
        cosine = np.sum(x * y, axis=1) / (np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1))
        angles = np.rad2deg(np.arccos(np.clip(cosine, -1, 1))) * 60
        ordering = np.argsort(angles)
        angles = angles[ordering]
        o0 = np.asarray(data["value"], dtype=float)[mask][ordering]
    fk, fold_provenance = fold_values(args.fk_fold, angles, SOURCE,
                                      Path(mc_manifest["table_path"]), Path(mc_manifest["response_path"]))
    provenance.update(fold_provenance)
    marker_provenance = json.loads(args.fk_markers.with_suffix(".json").read_text())
    if marker_provenance["fold_sha256"] != fold_provenance["fold_sha256"]:
        raise ValueError("Pooled FK markers were compared against a different fold")
    with np.load(args.fk_markers) as markers:
        columns = {"g": angles, "o0": o0, "fk3": fk,
                   "g_ff": reference_rows[:, 0], "ff_mom": reference_rows[:, 1],
                   "g_mc": mc_rows[:, 0], "ffc_mc": mc_rows[:, 1], "ffc_se": mc_rows[:, 2],
                   "g_fk": np.array(markers["gamma"]), "fk_mc": np.array(markers["fk"]),
                   "fk_se": np.array(markers["err"]), "fk_analytic": np.array(markers["analytic"])}
    if any(not np.all(np.isfinite(value)) for value in columns.values()):
        raise ValueError("Nonfinite figure data")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **columns, source_redshift=5.0, lambda_source=SOURCE,
             ff_reference_label="FF (local covariance)")
    args.out.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"output": str(args.out), "columns": {key: list(value.shape) for key, value in columns.items()}}, indent=2))


if __name__ == "__main__":
    main()
