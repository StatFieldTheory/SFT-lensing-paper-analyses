"""Prepare separately checkpointed FF source planes and optional O0 diagnostics.

The reviewed FK replay module supplies provenance, path safety and the shared
serial runner. This module only defines companion inputs and prepares fresh
configs. Its plan and prepare commands never start a numerical fold.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil

import numpy as np
import yaml

import replay_suite as replay

HERE = Path(__file__).resolve().parent
FF_DIR = replay.PRODUCTS / "C_corr_op_K_limber_FF"
O0_DIR = replay.PRODUCTS / "C_corr_op_O0"
MULTIZ_DIR = replay.PRODUCTS / "multiz"
FF_SIDECARS = {
    1.0: "multiz_ff_real.npz",
    1.7: "multiz_ff_real_z1p7.npz",
    2.5: "multiz_ff_real_z2p5.npz",
    3.2: "multiz_ff_real_z3p2.npz",
    4.0: "multiz_ff_real_z4p0.npz",
}


def source_specs() -> list[dict]:
    """Read actual figure redshifts and the original per-source FF endpoints."""
    main_config = yaml.safe_load((FF_DIR / "config_L2.yaml").read_text())
    main_lam = main_config["sweep"]["t_final_grid"]
    if len(main_lam) != 1:
        raise ValueError("Expected exactly one historical main source")
    specs = [{
        "source_id": "main", "z": 5.0, "historical_lambda_mpc": float(main_lam[0]),
        "historical_product": str(FF_DIR / "xi_C_corr_op_K_limber_FF.npz"),
        "group": "main", "component_pairs": main_config["sweep"]["component_pairs"],
    }]
    with np.load(replay.locate(replay.MULTIZ), allow_pickle=False) as data:
        redshifts = np.asarray(data["z"], float)
        fk_lambdas = np.asarray(data["lam"], float)
    with np.load(replay.locate(MULTIZ_DIR / "multiz_ff_real_all5.npz"), allow_pickle=False) as real:
        if not np.array_equal(real["z"], redshifts):
            raise ValueError("Actual plotted FF and figure redshift arrays disagree")
        for index, redshift in enumerate(redshifts):
            if float(redshift) not in FF_SIDECARS:
                raise ValueError(f"No traced historical FF sidecar for z={redshift}")
            product = MULTIZ_DIR / FF_SIDECARS[float(redshift)]
            with np.load(replay.locate(product), allow_pickle=False) as source:
                old_lam = np.asarray(source["lam"], float).reshape(-1)
                ff = np.asarray(source["ff"], float).reshape(-1)
                if old_lam.size != 1:
                    raise ValueError(f"Expected one saved FF endpoint: {product}")
                if not np.array_equal(ff, real["ff"][index]):
                    raise ValueError(f"FF sidecar differs from the actual plotted slice: {product}")
                if not np.array_equal(source["gamma"], real["gamma"]):
                    raise ValueError(f"FF sidecar angle grid differs: {product}")
            source_id = "z" + f"{redshift:.1f}".replace(".", "p")
            specs.append({
                "source_id": source_id, "z": float(redshift),
                "historical_lambda_mpc": float(old_lam[0]),
                "figure_fk_lambda_mpc": float(fk_lambdas[index]),
                "historical_product": str(product), "group": "multiz",
                "component_pairs": [[0, 0]],
            })
    return specs


def build_plan(mode: str, source_map: Path | None, targets: list[str],
               channels: list[str], profile: str, n_gauss: int, n_jobs: int) -> dict:
    replay.verify_canoes()
    sources = source_specs()
    known = {source["source_id"] for source in sources} | {"all", "multiz"}
    if set(targets) - known:
        raise ValueError(f"Unknown source targets: {set(targets) - known}")
    sources = [source for source in sources if "all" in targets or source["source_id"] in targets
               or ("multiz" in targets and source["group"] == "multiz")]
    mapping = replay.read_source_map(source_map) if mode == "true-redshift" else None
    if mode == "historical" and source_map is not None:
        raise ValueError("Historical controls must retain their original response module")
    if profile == "full" and n_gauss < 24:
        raise ValueError("GL12 is a probe setting, not a full figure setting")
    if not set(channels) <= {"ff", "o0-diagnostic"} or len(channels) != len(set(channels)):
        raise ValueError("Channels must be unique ff and/or o0-diagnostic selections")
    dependencies = [Path(__file__).resolve(), HERE / "replay_suite.py", HERE / "framework_provenance.json",
                    replay.REBUILD / "run_fk_variant.py", replay.SACHS / "scripts" / "run_FF_single.py",
                    replay.SACHS / "scripts" / "D_callable.py", replay.SACHS / "scripts" / "kappa2_callable.py",
                    replay.SACHS / "scripts" / "inputs" / "F_tensor.npy",
                    replay.SACHS / "scripts" / "inputs" / "kappa2_grid_lcut0_low_empty.npz",
                    replay.MULTIZ, MULTIZ_DIR / "multiz_ff_real_all5.npz"]
    corr = replay.SACHS / "callables" / "C_propagator" / "corr_op"
    dependencies += [corr / "corr_op_C_callable.py", corr / "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"]
    if mapping:
        dependencies += [Path(mapping["source_file"]), Path(mapping["response_module"])]
    targets_out = []
    for channel in channels:
        base = FF_DIR if channel == "ff" else O0_DIR
        config_file = base / "config_L2.yaml"
        config = yaml.safe_load(config_file.read_text())
        dependencies.append(config_file)
        if channel == "ff":
            # Retain the historical declaration, although F-only selection
            # excludes every K diagram. Its import-time dependencies still exist.
            dependencies += [replay.PLAIN / "equal_time_limber_kappa3_callable.py",
                             replay.PLAIN / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"]
        for source in sources:
            endpoint = source["historical_lambda_mpc"]
            if channel == "o0-diagnostic" and source["group"] == "multiz":
                endpoint = source["figure_fk_lambda_mpc"]
            if mapping:
                matches = np.flatnonzero(np.isclose(mapping["z"], source["z"], rtol=0, atol=1e-12))
                if len(matches) != 1:
                    raise ValueError(f"No unique reviewed source-map entry for z={source['z']}")
                endpoint = float(mapping["lambda_mpc"][int(matches[0])])
            dependencies.append(Path(source["historical_product"]))
            angles = config["sweep"]["positions_grid"]["y"]
            angle_indices = [0, len(angles) - 1] if profile == "probe" else list(range(len(angles)))
            name = channel.replace("-", "_") + "_" + source["source_id"]
            targets_out.append({
                **source, "id": name, "channel": channel, "source_config": str(config_file),
                "target_lambdas_mpc": [float(endpoint)], "n_jobs": n_jobs, "n_gauss": n_gauss,
                "orders": ",".join(str(order) for order in config["expand"]["orders"]),
                "angle_indices": angle_indices, "output_name": f"xi_{name}.npz",
                "expected_order": 2 if channel == "ff" else 0,
                "publication_role": "FF component" if channel == "ff" else "Diagnostic only. This does not replace the separately validated matched Order-0 reference.",
            })
    return {
        "schema_version": 1, "suite": "companion", "framework": replay.framework(),
        "source_mode": mode, "source_mapping": mapping, "profile": profile,
        "targets": targets_out,
        "inputs": [replay.file_record(path) for path in sorted(set(dependencies))],
        "covariance_limitation": "The historical 21-source-node C table is unchanged. Source alignment and outer quadrature convergence do not establish convergence of this tabulated input.",
        "checkpoint_policy": "One source and channel per target, private caches, completed-result hash checked on resume. An interrupted target requires a fresh directory.",
    }


def prepare(plan: dict, work: Path) -> None:
    work = work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    helper = replay.variant_module()
    generated = []
    for target in plan["targets"]:
        run_dir = work / target["id"]
        _, config_file = helper.materialise_variant(
            None, run_dir, t_final=target["target_lambdas_mpc"], n_jobs=target["n_jobs"],
            base_config=Path(target["source_config"]), orders=target["orders"],
        )
        config = yaml.safe_load(config_file.read_text())
        mapping = plan["source_mapping"]
        if mapping:
            config["system"]["linear"]["R_time_module"] = mapping["response_module"]
            config["system"]["linear"]["R_time_attr"] = mapping["response_attr"]
            step = float(config["propagators"]["dt"])
            config["propagators"]["t_max"] = step * math.ceil(max(target["target_lambdas_mpc"]) / step)
        angles = config["sweep"]["positions_grid"]["y"]
        config["sweep"]["positions_grid"]["y"] = [angles[index] for index in target["angle_indices"]]
        config["sweep"]["component_pairs"] = target["component_pairs"]
        config["sweep"]["n_gauss"] = target["n_gauss"]
        config["expand"]["n_jobs"] = 1
        replay.validate_config_paths(config, work)
        config_file.write_text(yaml.safe_dump(config, sort_keys=False))
        shutil.copyfile(target["source_config"], run_dir / "historical_config.yaml")
        target["prepared_config"] = str(config_file)
        target["result"] = str(run_dir / target["output_name"])
        generated += [replay.file_record(config_file), replay.file_record(run_dir / "historical_config.yaml")]
    plan["work_directory"] = str(work)
    plan["generated_inputs"] = generated
    replay.dump_new(work / "prepared.json", plan)
    print(f"Prepared {len(plan['targets'])} independently checkpointed targets, no folds started: {work}")


def validate_results(work: Path) -> dict:
    """Check complete row coverage without claiming scientific convergence."""
    plan = json.loads((work / "prepared.json").read_text())
    summary = []
    for target in plan["targets"]:
        result = Path(target["result"])
        config = yaml.safe_load(Path(target["prepared_config"]).read_text())
        wanted_y = {tuple(direction) for direction in config["sweep"]["positions_grid"]["y"]}
        wanted_x = np.asarray(config["sweep"]["positions_grid"]["x"][0], float)
        with np.load(result, allow_pickle=True) as data:
            required = {"a", "b", "order", "t_final", "value", "x", "y"}
            if not required <= set(data.files):
                raise ValueError(f"Missing output fields: {result}")
            expected = len(target["angle_indices"]) * len(target["component_pairs"])
            values = np.asarray(data["value"], float)
            if values.shape != (expected,) or not np.all(np.isfinite(values)):
                raise ValueError(f"Wrong row count or non-finite values: {result}")
            if not np.all(data["order"] == target["expected_order"]):
                raise ValueError(f"Unexpected diagram order: {result}")
            if not np.allclose(data["t_final"], target["target_lambdas_mpc"][0], rtol=0, atol=1e-10):
                raise ValueError(f"Unexpected source endpoint: {result}")
            pairs = np.column_stack([data["a"], data["b"]])
            for pair in target["component_pairs"]:
                mask = np.all(pairs == pair, axis=1)
                if np.count_nonzero(mask) != len(target["angle_indices"]):
                    raise ValueError(f"Incomplete component pair {pair}: {result}")
                directions = np.array([np.asarray(y, float) for y in data["y"][mask]])
                if {tuple(direction) for direction in directions} != wanted_y:
                    raise ValueError(f"Missing or unexpected angular rows: {result}")
                x_directions = np.array([np.asarray(x, float) for x in data["x"][mask]])
                if not np.all(x_directions == wanted_x):
                    raise ValueError(f"Unexpected reference direction: {result}")
        summary.append({"id": target["id"], "rows": expected, "status": "finite and complete", "result": replay.file_record(result)})
    return {"targets": summary, "scientific_convergence": "Not established by this structural check"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "prepare"):
        command = sub.add_parser(name)
        command.add_argument("--source-mode", choices=("historical", "true-redshift"), required=True)
        command.add_argument("--source-map", type=Path)
        command.add_argument("--targets", nargs="+", default=["all"])
        command.add_argument("--channels", nargs="+", choices=("ff", "o0-diagnostic"), default=["ff"])
        command.add_argument("--profile", choices=("probe", "full"), default="full")
        command.add_argument("--n-gauss", type=int, choices=(12, 24, 48), default=24)
        command.add_argument("--n-jobs", type=int, choices=(1, 2, 4, 6), default=1)
        if name == "prepare":
            command.add_argument("--work-dir", type=Path, required=True)
    for name in ("run", "validate"):
        command = sub.add_parser(name)
        command.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "run":
        replay.run_batch(args.work_dir)
        print(json.dumps(validate_results(args.work_dir.resolve()), indent=2))
    elif args.command == "validate":
        print(json.dumps(validate_results(args.work_dir.resolve()), indent=2))
    else:
        plan = build_plan(args.source_mode, args.source_map, args.targets, args.channels,
                          args.profile, args.n_gauss, args.n_jobs)
        if args.command == "prepare":
            prepare(plan, args.work_dir)
        else:
            print(json.dumps(plan, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
