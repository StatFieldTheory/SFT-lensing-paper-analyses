"""Assemble completed, aligned convergence products into fresh plot inputs.

No numerical models, interpolation or missing-component synthesis are performed.
All inputs are validated before the output directory is created. The resulting
assembly manifest is evidence for the caller, not the active product manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sys

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ALIGNED = HERE.parent / "analysis1" / "r1_aligned"
sys.path.insert(0, str(HERE.parents[2] / "reproduce"))
from archived_inputs import locate  # noqa: E402
REDSHIFTS = np.asarray([1.0, 1.7, 2.5, 3.2, 4.0])
CUTS = (960, 1920, 3840, 7680, 15360)
FRAMEWORK_COMMIT = "6945815ad9e6efa2b9e768891d6bf92b1214670d"


class Inputs:
    """Collect and verify every consumed input, including transitive records."""

    def __init__(self):
        self.records: dict[str, dict] = {}

    def record(self, path: Path, expected: dict | None = None) -> dict:
        path = path.resolve()
        stored = locate(path)  # an archived input is hashed in the archive
        with stored.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        record = {"path": str(path), "sha256": digest, "bytes": stored.stat().st_size}
        if expected is not None and (expected["sha256"] != digest or
                                    expected.get("bytes", record["bytes"]) != record["bytes"]):
            raise ValueError(f"Input hash or size mismatch: {path}")
        previous = self.records.get(str(path))
        if previous is not None and previous != record:
            raise ValueError(f"Input changed during assembly: {path}")
        self.records[str(path)] = record
        return record

    def read_json(self, path: Path, expected: dict | None = None) -> dict:
        self.record(path, expected)
        return json.loads(path.read_text())


def finite_array(value, shape: tuple, label: str) -> np.ndarray:
    array = np.asarray(value)
    if array.shape != shape or not np.isrealobj(array) or not np.all(np.isfinite(array)):
        raise ValueError(f"Invalid shape or nonfinite/complex values for {label}: {array.shape}")
    return array.astype(float)


def source_map(path: Path, inputs: Inputs) -> dict:
    mapping = inputs.read_json(path)
    z = np.asarray(mapping["z"], float)
    lam = finite_array(mapping["lambda_mpc"], z.shape, "source-map lambda")
    if z.ndim != 1 or not np.all(np.isfinite(z)) or not np.all(np.diff(z) > 0) or not np.all(np.diff(lam) > 0):
        raise ValueError("Source map must have finite strictly increasing redshift and affine arrays")
    for filename, digest in mapping["mapping_provenance"]["source_sha256"].items():
        inputs.record(Path(filename), {"sha256": digest})
    provenance = mapping["mapping_provenance"]
    inputs.record(Path(provenance["C_table_path"]), {"sha256": provenance["C_table_sha256"]})
    return mapping


def endpoint(mapping: dict, redshift: float) -> float:
    matches = np.flatnonzero(np.asarray(mapping["z"]) == redshift)
    if len(matches) != 1:
        raise ValueError(f"No unique source-map entry for z={redshift}")
    return float(mapping["lambda_mpc"][int(matches[0])])


def read_batch(path: Path, mapping: dict, inputs: Inputs, kind: str) -> dict:
    plan = inputs.read_json(path / "prepared.json")
    if plan["source_mode"] != "true-redshift":
        raise ValueError(f"Historical-endpoint batch cannot supply revised figures: {path}")
    if plan["framework"]["version"] != "0.6.1" or plan["framework"]["commit"] != FRAMEWORK_COMMIT:
        raise ValueError(f"Unrecognized framework provenance: {path}")
    for key in ("z", "lambda_mpc", "response_module", "response_attr", "mapping_provenance"):
        if plan["source_mapping"][key] != mapping[key]:
            raise ValueError(f"Batch source mapping differs in {key}: {path}")
    if kind == "fk" and plan.get("input_mode") != "corrected-nphi512":
        raise ValueError("FK figure inputs require the explicitly corrected vertex tables")
    if kind == "ff" and (plan.get("suite") != "companion" or plan.get("profile") != "full"):
        raise ValueError("FF figure inputs require the full companion profile")
    for record in plan["inputs"] + plan["generated_inputs"]:
        inputs.record(Path(record["path"]), record)
    return plan


def read_rows(target: dict, plan: dict, inputs: Inputs, *, main_pairs: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read completed order-two convergence rows in their exact config order."""
    result = Path(target["result"])
    completed = inputs.read_json(result.parent / "completed.json")
    if completed["framework"] != plan["framework"] or completed["source_mode"] != "true-redshift":
        raise ValueError(f"Completion provenance differs from the prepared batch: {result}")
    if Path(completed["result"]["path"]).resolve() != result.resolve():
        raise ValueError("Completion record points to a different result")
    inputs.record(result, completed["result"])
    config_path = Path(target["prepared_config"])
    inputs.record(config_path)
    config = yaml.safe_load(config_path.read_text())
    if target["n_gauss"] != config["sweep"]["n_gauss"]:
        raise ValueError("Runtime target quadrature differs from its pinned prepared config")
    expected_vertices = ["F"] if target.get("channel") == "ff" else ["FK"]
    if config["sweep"]["vertex_types"] != expected_vertices:
        raise ValueError("Prepared diagram selection does not match the requested FF or FK product")
    expected_pairs = [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]] if main_pairs else [[0, 0]]
    if config["system"]["field"]["n_components"] != 3 or config["sweep"]["component_pairs"] != expected_pairs:
        raise ValueError("Expected all three internal fields and the requested observable pairs")
    positions = config["sweep"]["positions_grid"]
    x = finite_array(positions["x"], (1, 3), "reference direction")[0]
    y = np.asarray(positions["y"], float)
    y = finite_array(y, (len(y), 3), "angular directions")
    if np.linalg.norm(x) == 0 or np.any(np.linalg.norm(y, axis=1) == 0):
        raise ValueError("Zero-length external direction")
    gamma = np.degrees(np.arccos(np.clip((y / np.linalg.norm(y, axis=1)[:, None]) @
                                        (x / np.linalg.norm(x)), -1, 1))) * 60
    if not np.all(np.diff(gamma) > 0):
        raise ValueError("Prepared angular grid must be strictly increasing")
    lambdas = np.asarray(target["target_lambdas_mpc"], float)
    if not np.array_equal(config["sweep"]["t_final_grid"], lambdas):
        raise ValueError("Prepared target and config source endpoints disagree")
    count = len(y) * len(lambdas)
    # Framework output positions use object arrays. They are read only after
    # matching this local product to its completed-result cryptographic hash.
    with np.load(result, allow_pickle=True) as data:
        total = count * len(expected_pairs)
        all_values = finite_array(data["value"], (total,), "fold values")
        for key in ("a", "b", "order"):
            if data[key].shape != (total,):
                raise ValueError(f"Unexpected {key} row count: {result}")
        if not np.all(data["order"] == 2):
            raise ValueError(f"Unexpected perturbative order: {result}")
        for a, b in expected_pairs:
            if np.count_nonzero((data["a"] == a) & (data["b"] == b)) != count:
                raise ValueError(f"Incomplete observable pair {(a, b)}: {result}")
        selected = (data["a"] == 0) & (data["b"] == 0)
        values = all_values[selected]
        row_x = finite_array(np.asarray(list(data["x"]), float), (total, 3), "result x")[selected]
        row_y = finite_array(np.asarray(list(data["y"]), float), (total, 3), "result y")[selected]
        row_lam = finite_array(data["t_final"], (total,), "result endpoints")[selected]
    if not np.all(row_x == x):
        raise ValueError("Result reference directions differ from the prepared config")
    lookup = {(float(lam), tuple(direction)): (i, j)
              for i, lam in enumerate(lambdas) for j, direction in enumerate(y)}
    if len(lookup) != count:
        raise ValueError("Prepared source/angle coordinates are not unique")
    matrix = np.empty((len(lambdas), len(y)))
    seen = set()
    for lam, direction, value in zip(row_lam, row_y, values, strict=True):
        row_key = (float(lam), tuple(direction))
        if row_key not in lookup or row_key in seen:
            raise ValueError("Duplicate or unrequested source/angle row in fold product")
        seen.add(row_key)
        matrix[lookup[row_key]] = value
    if len(seen) != count:
        raise ValueError("Incomplete fold product")
    return x, y, gamma, lambdas, matrix


def read_order0(path: Path, mapping: dict, gamma: np.ndarray, inputs: Inputs) -> np.ndarray:
    acceptance = inputs.read_json(path)
    if acceptance["status"] != "accepted_for_multiz_order0":
        raise ValueError("Order-0 acceptance manifest has not passed its publication gate")
    accepted_mapping = inputs.read_json(Path(acceptance["source_map"]["path"]), acceptance["source_map"])
    if accepted_mapping != mapping:
        raise ValueError("Order-0 acceptance pins a different source map")
    validation_record = acceptance["validation"]
    validation = inputs.read_json(path.parent / validation_record["path"], validation_record)
    if validation["status"] != "passed":
        raise ValueError("Order-0 validation did not pass")
    expected_keys = {str(float(z)) for z in REDSHIFTS}
    if set(acceptance["products"]) != expected_keys:
        raise ValueError("Order-0 acceptance must contain exactly the five required source planes")
    values = []
    for z in REDSHIFTS:
        entry = acceptance["products"][str(float(z))]
        if entry["z_source"] != z:
            raise ValueError("Order-0 acceptance redshift key and entry disagree")
        product = path.parent / entry["path"]
        inputs.record(product, entry)
        metadata = inputs.read_json(path.parent / entry["metadata_path"], {"sha256": entry["metadata_sha256"]})
        with np.load(product, allow_pickle=False) as data:
            if json.loads(str(data["metadata_json"])) != metadata:
                raise ValueError("Order-0 embedded metadata and accepted sidecar disagree")
            if not np.array_equal(data["theta_arcmin"], gamma):
                raise ValueError("Order-0 and framework angle grids differ; no interpolation is allowed")
            values.append(finite_array(data["kk"], gamma.shape, "Order-0 convergence"))
        if metadata["observable_operator"] != "full" or metadata["source_redshift"] != z:
            raise ValueError("Order-0 operator or source redshift mismatch")
        if metadata["source_lambda_mpc"] != endpoint(mapping, float(z)):
            raise ValueError("Order-0 source affine endpoint mismatch")
        if metadata["source"]["mapping_provenance"] != mapping["mapping_provenance"]:
            raise ValueError("Order-0 background/source mapping differs from the framework")
        snapshot = path.parent / entry["source_snapshot_path"]
        if entry["source_snapshot_sha256"] != metadata["script_sha256"]:
            raise ValueError("Accepted source snapshot and product script hash disagree")
        inputs.record(snapshot, {"sha256": entry["source_snapshot_sha256"]})
        inputs.record(Path(metadata["power_table"]), {"sha256": metadata["power_table_sha256"]})
    return np.asarray(values)


def write_products(directory: Path, payloads: dict[str, dict], inputs: Inputs, details: dict) -> dict:
    if directory.exists():
        raise FileExistsError(f"Output directory must be fresh: {directory}")
    for name in payloads:
        if Path(name).name != name or Path(name).suffix != ".npz":
            raise ValueError("Product filenames must be local basenames with an explicit .npz suffix")
    inputs.record(Path(__file__))
    # Recheck after loading, before any writes, so provenance cannot silently
    # combine inputs that changed while this assembly was reading them.
    for record in list(inputs.records.values()):
        inputs.record(Path(record["path"]), record)
    provenance = {"kind": "validated convergence figure product assembly", "details": details,
                  "inputs": sorted(inputs.records.values(), key=lambda item: item["path"]),
                  "policy": "Exact existing grids and values, no interpolation or synthesized components"}
    directory.mkdir(parents=True, exist_ok=False)
    outputs = []
    for name, payload in payloads.items():
        path = directory / name
        with path.open("xb") as stream:
            np.savez(stream, **payload, metadata_json=json.dumps(provenance, sort_keys=True))
        outputs.append(inputs.record(path))
    manifest = {**provenance, "outputs": outputs}
    with (directory / "assembly_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    return manifest


def multiz(args: argparse.Namespace) -> dict:
    inputs = Inputs()
    mapping = source_map(args.source_map, inputs)
    fk_plan = read_batch(args.fk_batch, mapping, inputs, "fk")
    fk_targets = [target for target in fk_plan["targets"] if target["id"] == "multiz"]
    if len(fk_targets) != 1 or not np.array_equal(fk_targets[0]["nominal_z"], REDSHIFTS):
        raise ValueError("FK batch must contain one complete five-plane multiz target")
    x, y, gamma, lambdas, fk = read_rows(fk_targets[0], fk_plan, inputs)
    if len(gamma) != 40 or not np.array_equal(lambdas, [endpoint(mapping, z) for z in REDSHIFTS]):
        raise ValueError("FK source or 40-angle grid is incomplete")
    ff_plan = read_batch(args.ff_batch, mapping, inputs, "ff")
    ff = []
    for z in REDSHIFTS:
        targets = [target for target in ff_plan["targets"] if target.get("channel") == "ff" and target["z"] == z]
        if len(targets) != 1:
            raise ValueError(f"Expected exactly one genuine FF target at z={z}")
        fx, fy, fg, fl, fv = read_rows(targets[0], ff_plan, inputs)
        if not all(np.array_equal(a, b) for a, b in ((x, fx), (y, fy), (gamma, fg),
                                                     (np.asarray([endpoint(mapping, z)]), fl))):
            raise ValueError("FF and FK external coordinates differ; no interpolation is allowed")
        ff.append(fv[0])
    ff = np.asarray(ff)
    o0 = read_order0(args.o0_manifest, mapping, gamma, inputs)
    common = {"z": REDSHIFTS, "gamma": gamma, "lam": lambdas}
    return write_products(args.out_dir, {
        "multiz.npz": {**common, "o0": o0, "ff": ff, "fk": fk},
        "multiz_ff.npz": {**common, "ff": ff},
    }, inputs, {"sources": REDSHIFTS.tolist(), "lambda_mpc": lambdas.tolist(), "angle_count": 40})


def read_guard(path: Path, inputs: Inputs) -> dict:
    """Require the exact-query guard test and all of its pinned inputs."""
    validation = inputs.read_json(path)
    if (validation.get("status") != "passed" or validation.get("array_equal") is not True or
            validation.get("scalar_tensor_comparisons") != 2728 or
            validation.get("batch_tensor_comparisons") != 2728 or
            set(validation.get("rejections", [])) != {"one_ulp", "mixed_angles", "nonfinite_time"}):
        raise ValueError("Permutation-aware guard validation is incomplete or failed")
    for filename, digest in validation["sha256"].items():
        inputs.record(Path(filename), {"sha256": digest})
    contract_path = path.parent / "contract.json"
    required = [contract_path, path.parent / "guard.py",
                path.parent / "table_tree_cut15360.npz",
                path.parent / "folds_tree/tree_15360/callable_tree_15360.py"]
    if not all(str(required_path.resolve()) in validation["sha256"] for required_path in required):
        raise ValueError("Guard validation does not pin all required contract, guard and tested vertex files")
    contract = inputs.read_json(contract_path)
    for key in ("source", "loader", "source_map", "original_config", "framework_expansion_cache"):
        inputs.record(Path(contract[key]["path"]), contract[key])
    return contract


def physics_config(path: Path, inputs: Inputs) -> dict:
    """Remove only output locations, worker counts and the selected sky pairs."""
    inputs.record(path)
    config = yaml.safe_load(path.read_text())
    if len(config["system"]["nonlocal_vertices"]) != 1:
        raise ValueError("Expected exactly one nonlocal K vertex")
    config["system"]["nonlocal_vertices"][0].pop("coupling_module")
    for key in ("expand", "propagators", "sweep"):
        config[key].pop("cache_path", None)
        config[key].pop("n_jobs", None)
    config["sweep"].pop("positions_grid")
    config["sweep"].pop("component_pairs")
    return {key: config[key] for key in ("system", "expand", "propagators", "sweep")}


def main_fold_gate(args: argparse.Namespace, mapping: dict, inputs: Inputs,
                   contract: dict, reference: tuple, restricted: np.ndarray) -> dict:
    """Compare the 11 exact directions to the full accepted main FK fold."""
    plan = read_batch(args.main_batch, mapping, inputs, "fk")
    targets = [target for target in plan["targets"] if target["id"] == "main"]
    if len(targets) != 1 or targets[0]["nominal_z"] != [5.0]:
        raise ValueError("Expected one accepted full main target at true z=5")
    target = targets[0]
    if Path(target["table"]).resolve() != Path(contract["source"]["path"]).resolve():
        raise ValueError("Exact-query contract does not reference the accepted main vertex")
    x, y, gamma, lambdas, values = read_rows(target, plan, inputs, main_pairs=True)
    indices = np.asarray(contract["selected_indices"])
    if (indices.shape != (11,) or not np.issubdtype(indices.dtype, np.integer) or
            len(np.unique(indices)) != 11 or np.any(indices < 0) or np.any(indices >= len(y))):
        raise ValueError("Invalid main-fold direction indices in the query contract")
    if len(y) != 40 or not np.array_equal(lambdas, [endpoint(mapping, 5.0)]):
        raise ValueError("Main FK fold has a different source or angular grid")
    if not np.array_equal(contract["x"], x) or not np.array_equal(contract["y"], y[indices]):
        raise ValueError("Query contract directions differ from the exact main-fold subset")
    for wanted, actual in zip((x, y[indices], gamma[indices]), reference, strict=True):
        if not np.array_equal(wanted, actual):
            raise ValueError("Cutoff ladder and accepted main fold use different directions")
    selected = values[0, indices]
    if not np.array_equal(selected.view(np.uint64), restricted.view(np.uint64)):
        delta = float(np.max(np.abs(selected - restricted)))
        raise ValueError(f"Highest tree cutoff is not bitwise equal to the main fold: max absolute difference {delta:.17g}")
    return {"status": "passed", "comparison": "bitwise", "angle_count": 11,
            "main_result": str(Path(target["result"]).resolve()),
            "selected_main_indices": indices.tolist(), "max_absolute_difference": 0.0}


def ladder(args: argparse.Namespace) -> dict:
    inputs = Inputs()
    mapping = source_map(args.source_map, inputs)
    lam = endpoint(mapping, 5.0)
    contract = read_guard(args.guard_validation, inputs)
    if inputs.read_json(Path(contract["source_map"]["path"]), contract["source_map"]) != mapping:
        raise ValueError("Query contract pins a different source map")
    main_plan = read_batch(args.main_batch, mapping, inputs, "fk")
    main_targets = [target for target in main_plan["targets"] if target["id"] == "main"]
    if len(main_targets) != 1:
        raise ValueError("Expected exactly one accepted main target")
    main_physics = physics_config(Path(main_targets[0]["prepared_config"]), inputs)
    payloads = {}
    reference = None
    for model, batch_path in (("tree", args.tree_batch), ("bihalofit", args.bihalofit_batch)):
        plan = read_batch(batch_path, mapping, inputs, "fk")
        if plan.get("ladder_geometry") != "permaware34":
            raise ValueError("Revised cutoff products require the explicit restricted query contract")
        for cutoff in CUTS:
            targets = [target for target in plan["targets"] if target["id"] == f"{model}_{cutoff}"]
            if len(targets) != 1:
                raise ValueError(f"Missing or duplicate {model} cutoff {cutoff}")
            if physics_config(Path(targets[0]["prepared_config"]), inputs) != main_physics:
                raise ValueError("Cutoff fold physics or quadrature differs from the accepted main fold")
            table_path = Path(targets[0]["table"]).resolve()
            if str(table_path) not in inputs.records:
                raise ValueError("Cutoff table is not pinned in the prepared batch")
            with np.load(table_path, allow_pickle=False) as vertex:
                metadata = json.loads(np.asarray(vertex["cosmo_meta"], dtype=np.uint8).tobytes().decode())
                if int(vertex["ell_max"]) != cutoff or metadata["b_delta_model"] != model:
                    raise ValueError("Cutoff target and hashed vertex table model/cutoff disagree")
            x, y, gamma, lambdas, matrix = read_rows(targets[0], plan, inputs)
            if len(gamma) != 11 or not np.array_equal(lambdas, [lam]):
                raise ValueError("Cutoff product has wrong angle count or source endpoint")
            if not np.array_equal(x, contract["x"]) or not np.array_equal(y, contract["y"]):
                raise ValueError("Cutoff directions differ from the validated guard contract")
            current = (x, y, gamma)
            if reference is None:
                reference = current
            elif not all(np.array_equal(a, b) for a, b in zip(reference, current, strict=True)):
                raise ValueError("Cutoff products use inconsistent external directions")
            count = len(gamma)
            payloads[f"table_{model}_cut{cutoff}_r4_xi.npz"] = {
                "x": np.tile(x, (count, 1)), "y": y, "t_final": np.full(count, lam),
                "a": np.zeros(count, dtype=int), "b": np.zeros(count, dtype=int),
                "order": np.full(count, 2, dtype=int), "value": matrix[0],
            }
    acceptance = main_fold_gate(args, mapping, inputs, contract, reference,
                                payloads["table_tree_cut15360_r4_xi.npz"]["value"])
    details = {
        "source_redshift": 5.0, "lambda_mpc": lam, "cutoffs": list(CUTS),
        "angle_count": 11, "component_pair": [0, 0],
        "main_fold_acceptance": acceptance,
        "guard_validation": str(args.guard_validation.resolve()),
        "filename_note": "r4 filenames retain plotting compatibility; actual geometry is the validated permutation-aware 34-row contract",
    }
    if args.validate_only:
        return {"status": "validated", "details": details, "outputs": [],
                "inputs": sorted(inputs.records.values(), key=lambda item: item["path"])}
    return write_products(args.out_dir, payloads, inputs, details)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("multiz")
    command.add_argument("--o0-manifest", type=Path, default=ALIGNED / "source_folded_multiz_manifest.json")
    command.add_argument("--fk-batch", type=Path, default=HERE / "true_redshift_fk_gl24_corrected_multiz")
    command.add_argument("--ff-batch", type=Path, default=HERE / "true_redshift_ff_gl24_all")
    command = sub.add_parser("ladder")
    command.add_argument("--tree-batch", type=Path, default=HERE / "permaware_ladder" / "folds_tree")
    command.add_argument("--bihalofit-batch", type=Path, default=HERE / "permaware_ladder" / "folds_bihalofit")
    command.add_argument("--main-batch", type=Path, default=HERE / "true_redshift_fk_gl24_corrected_main")
    command.add_argument("--guard-validation", type=Path, default=HERE / "permaware_ladder" / "guard_validation.json")
    command.add_argument("--validate-only", action="store_true", help="Validate all inputs without creating products")
    for command in sub.choices.values():
        command.add_argument("--source-map", type=Path, default=HERE / "source_map.json")
        command.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = {"multiz": multiz, "ladder": ladder}[args.command](args)
    print(json.dumps({"status": manifest.get("status", "assembled"), "outputs": manifest["outputs"],
                      "details": manifest.get("details", {})}, indent=2))


if __name__ == "__main__":
    main()
