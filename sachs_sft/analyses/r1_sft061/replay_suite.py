"""Prepare and serially replay the historical FK inputs with sft-wick 0.6.1.

Planning and preparation never start a fold. The run subcommand is explicit,
requires a fresh prepared batch, and refuses to overlap another fold. All writes
stay inside a new batch directory. Historical production files remain inputs.
Those archived on 2026-10-02 are read through reproduce/archived_inputs.py,
which checks each archived copy against its recorded hash; a record keeps the
original path.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
SACHS = HERE.parents[1]
REPO = SACHS.parent
sys.path.insert(0, str(REPO / "reproduce"))
from archived_inputs import locate  # noqa: E402
PRODUCTS = SACHS / "sftwick_outputs" / "2PCF"
REBUILD = SACHS / "callables" / "kappa3_vertex" / "rebuild"
PERM = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"
PLAIN = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber"
MAIN = PRODUCTS / "C_corr_op_K_limber_FK_cut15360_permfix"
MULTIZ = SACHS / "analyses" / "analysis3" / "outputs" / "multiz_kappa_2pcf_5z.npz"
CUTS = (960, 1920, 3840, 7680, 15360)
PIN = "6945815ad9e6efa2b9e768891d6bf92b1214670d"
NP512 = REBUILD / "products" / "pieces_nphi512"
CORRECTED_LADDER = HERE / "corrected_kappa3"
RESTRICTED_LADDER = HERE / "corrected_kappa3_ladder52"
GUARD = HERE / "restricted_kappa3_guard.py"


def guard_module():
    spec = importlib.util.spec_from_file_location("r1_restricted_kappa3_guard", GUARD)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def file_record(path: Path) -> dict:
    """Record a file under its recorded path, hashing an archived input in the archive."""
    path = path.resolve()
    stored = locate(path)
    return {"path": str(path), "sha256": digest(stored), "bytes": stored.stat().st_size}


def dump_new(path: Path, value: dict) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def framework() -> dict:
    """Require the exact installed release already verified by the parent task."""
    package = importlib.util.find_spec("sft_wick")
    if package is None or package.origin is None:
        raise RuntimeError("sft_wick is not importable in this interpreter")
    source = Path(package.origin).resolve().parents[2]
    commit = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"],
        text=True,
    ).strip()
    version = importlib.metadata.version("sft-wick")
    if version != "0.6.1" or commit != PIN or dirty:
        raise RuntimeError(f"Expected clean sft-wick 0.6.1 at {PIN}, got {version}, {commit}, {dirty!r}")
    return {"version": version, "commit": commit, "source": str(source), "python": sys.executable}


def verify_canoes() -> Path:
    """Enforce the separately recorded loader state without modifying it."""
    expected = json.loads((HERE / "framework_provenance.json").read_text())["canoes"]
    source = Path(expected["source"]).resolve()
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    diff = subprocess.check_output(["git", "-C", str(source), "diff", "HEAD", "--binary"])
    if commit != expected["commit"] or hashlib.sha256(diff).hexdigest() != expected["working_state"]["tracked_diff_sha256"]:
        raise RuntimeError("canoes loader state changed from framework_provenance.json")
    override = os.environ.get("CANOES_ROOT")
    if override and Path(override).resolve() != source:
        raise RuntimeError(f"CANOES_ROOT disagrees with recorded loader source: {override}")
    return source


def nominal_specs() -> list[dict]:
    """Trace the figure inputs rather than the superseded compute defaults."""
    main_config = MAIN / "config_L2.yaml"
    baseline = yaml.safe_load(main_config.read_text())["sweep"]["t_final_grid"]
    with np.load(locate(MULTIZ), allow_pickle=False) as data:
        multi_lam = np.asarray(data["lam"], float).tolist()
        multi_z = np.asarray(data["z"], float).tolist()
    specs = [{
        "id": "main", "group": "main", "nominal_z": [5.0],
        "historical_lambdas_mpc": baseline,
        "source_config": str(main_config),
        "record_config": str(MAIN / "config_record_2026-08-26.yaml"),
        "table": str(PERM / "table_permclosed.npz"),
        "callable": str(PERM / "perm_aware_kappa3_callable.py"),
        "historical_product": str(MAIN / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"),
        "output_name": "xi_main.npz", "reproduction_figure_ids": [2, 3, 4, 6],
        "provenance": "Production run README and path-migrated config_L2.yaml.",
    }, {
        "id": "multiz", "group": "multiz", "nominal_z": multi_z,
        "historical_lambdas_mpc": multi_lam,
        "source_config": str(PRODUCTS / "multiz" / "configs" / "config_fk_cut15360_permfix_5planes_record_2026-08-26.yaml"),
        "record_config": str(PRODUCTS / "multiz" / "configs" / "config_fk_cut15360_permfix_5planes_record_2026-08-26.yaml"),
        "table": str(PERM / "table_permclosed.npz"),
        "callable": str(PERM / "perm_aware_kappa3_callable.py"),
        "historical_product": str(MULTIZ), "output_name": "xi_multiz_fk.npz",
        "reproduction_figure_ids": [5],
        "provenance": "multiz/README.md and revision_2026-08/regenerate_multiz.py identify the perm-aware production vertex. Exact runtime endpoints come from the deployed NPZ because the recorded YAML predates the runtime override.",
    }]
    for model in ("tree", "bihalofit"):
        for cutoff in CUTS:
            name = f"table_{model}_cut{cutoff}_r4"
            config = PRODUCTS / "cutoff_ladder" / "configs" / f"config_{model}_cut{cutoff}_r4_record.yaml"
            old_lambda = yaml.safe_load(config.read_text())["sweep"]["t_final_grid"]
            specs.append({
                "id": f"{model}_{cutoff}", "group": "ladder", "nominal_z": [5.0],
                "historical_lambdas_mpc": old_lambda,
                "source_config": str(config), "record_config": str(config),
                "table": str(REBUILD / "products" / f"{name}.npz"),
                "callable": str(PLAIN / "equal_time_limber_kappa3_callable.py"),
                "historical_product": str(PRODUCTS / "cutoff_ladder" / f"{name}_xi.npz"),
                "output_name": f"{name}_xi.npz", "reproduction_figure_ids": [17],
                "provenance": "cutoff_ladder/README.md and rebuild/sweep_cutoff.sh identify the sorting callable and each saved r4 table. These inputs differ from the main permutation-aware run.",
            })
    return specs


def read_source_map(path: Path | None) -> dict:
    if path is None:
        raise ValueError("true-redshift mode requires --source-map with reviewed geometry and response module")
    mapping = json.loads(path.read_text())
    z = np.asarray(mapping["z"], float)
    lam = np.asarray(mapping["lambda_mpc"], float)
    if (z.ndim != 1 or z.shape != lam.shape or len(z) < 1
            or not np.all(np.isfinite(z)) or not np.all(np.isfinite(lam))
            or not np.all(np.diff(z) > 0) or not np.all(np.diff(lam) > 0)
            or np.any(z <= 0) or np.any(lam <= 0)):
        raise ValueError("Source map requires finite, positive, strictly increasing z and lambda arrays")
    if not mapping.get("mapping_provenance"):
        raise ValueError("Source map must explain its coordinate and background provenance")
    response = Path(mapping["response_module"])
    if not response.is_absolute():
        response = path.resolve().parent / response
    mapping["response_module"] = str(response.resolve())
    file_record(response)
    mapping["response_attr"] = mapping.get("response_attr", "response_R")
    mapping["source_file"] = str(path.resolve())
    return mapping


def select_corrected_input(spec: dict, contract: dict | None = None) -> None:
    """Select only explicitly corrected tables with unchanged sampling geometry."""
    historical_table = Path(spec["table"])
    if spec["group"] in ("main", "multiz"):
        table = NP512 / "table_permclosed_np512.npz"
    else:
        model, cutoff = spec["id"].rsplit("_", 1)
        table = (RESTRICTED_LADDER / f"table_{model}_cut{cutoff}_np512.npz" if contract
                 else CORRECTED_LADDER / f"table_{model}_cut{cutoff}_r4_np512.npz")
    if not table.is_file():
        raise FileNotFoundError(
            f"Corrected input has not been prepared: {table}. "
            "See CORRECTED_INPUTS.md. Historical tables are never substituted."
        )
    with np.load(table, allow_pickle=False) as corrected, np.load(locate(historical_table), allow_pickle=False) as old:
        meta = json.loads(np.asarray(corrected["cosmo_meta"], dtype=np.uint8).tobytes().decode())
        if meta.get("band_quadrature") != [[96, 512]]:
            raise ValueError(f"Expected n_ell=96 and n_phi=512 in every HIGH band: {table}")
        for key in ("cosine_triples", "lambda_shells_Mpc"):
            expected = old[key]
            if contract and spec["group"] == "ladder" and key == "cosine_triples":
                expected = expected[contract["raw_row_indices"]]
            if not np.array_equal(corrected[key], expected):
                raise ValueError(f"Corrected input changes the recorded {key}: {table}")
        if int(corrected["ell_max"]) != int(old["ell_max"]):
            raise ValueError(f"Corrected input cutoff differs from the recorded cutoff: {table}")
        shape = (len(corrected["cosine_triples"]), len(corrected["lambda_shells_Mpc"]))
        for channel in ("TTT", "TTP", "TPP", "PPP", "Bmod", "Dmod"):
            values = corrected[f"zeta_{channel}"]
            if values.shape != shape or not np.all(np.isfinite(values)):
                raise ValueError(f"Invalid corrected {channel} cells: {table}")
        expected_model = "bihalofit" if spec["id"].startswith("bihalofit_") else "tree"
        if meta.get("b_delta_model") != expected_model:
            raise ValueError(f"Wrong corrected bispectrum model: {table}")
        low_piece = ((RESTRICTED_LADDER / "low_ladder52_np512.npz" if contract else
                      CORRECTED_LADDER / "low_r4_np512.npz") if expected_model == "bihalofit"
                     else NP512 / "low_permclosed_r4_np512.npz")
        if meta.get("assembled_from", {}).get("low") != low_piece.name:
            raise ValueError(f"Corrected table has an unrecognized LOW input: {table}")
        with np.load(low_piece, allow_pickle=False) as low:
            build = json.loads(np.asarray(low["build"], dtype=np.uint8).tobytes().decode())
            if build.get("Nmax") != 384:
                raise ValueError(f"Corrected LOW input must have Nmax=384: {low_piece}")
    if contract and spec["group"] == "ladder":
        guard_module().validate_table(table, contract)
        spec["query_contract"] = str(RESTRICTED_LADDER / "query_contract.json")
        spec["selected_angle_indices"] = contract["selected_indices"]
    spec["historical_table"] = str(historical_table)
    spec["table"] = str(table)
    spec["input_low_piece"] = str(low_piece)
    spec["input_correction"] = (
        "September pair-phase correction, HIGH n_phi=512 and LOW Nmax=384. " +
        ("Exact contracted query neighborhoods retain all raw duplicate members, radial shells, cutoff and sorting callable."
         if contract and spec["group"] == "ladder" else
         "Recorded spatial rows, radial shells, cutoff and callable are retained.")
    )


def build_plan(mode: str, source_map: Path | None, targets: list[str], n_jobs: int,
               input_mode: str = "historical", ladder_geometry: str = "historical") -> dict:
    verify_canoes()
    specs = nominal_specs()
    known = {spec["id"] for spec in specs} | {"all", "ladder"}
    if set(targets) - known:
        raise ValueError(f"Unknown targets: {set(targets) - known}")
    specs = [spec for spec in specs if "all" in targets or spec["id"] in targets
             or ("ladder" in targets and spec["group"] == "ladder")]
    contract = None
    if ladder_geometry == "contract52":
        if input_mode != "corrected-nphi512":
            raise ValueError("The restricted ladder contract requires corrected-nphi512 inputs")
        contract = guard_module().validate_contract(RESTRICTED_LADDER / "query_contract.json")
    elif ladder_geometry != "historical":
        raise ValueError(f"Unknown ladder geometry: {ladder_geometry}")
    if input_mode == "corrected-nphi512":
        with np.load(NP512 / "low_permclosed_r4_np512.npz", allow_pickle=False) as low:
            low_build = json.loads(np.asarray(low["build"], dtype=np.uint8).tobytes().decode())
            if low_build.get("Nmax") != 384:
                raise ValueError("Corrected LOW input must carry the recorded Nmax=384 build")
        for spec in specs:
            select_corrected_input(spec, contract if spec["group"] == "ladder" else None)
    elif input_mode != "historical":
        raise ValueError(f"Unknown input mode: {input_mode}")
    mapping = read_source_map(source_map) if mode == "true-redshift" else None
    if mode == "historical" and source_map is not None:
        raise ValueError("historical mode cannot accept a replacement source map")
    inputs = [HERE / "framework_provenance.json", REBUILD / "run_fk_variant.py",
              SACHS / "scripts" / "run_FF_single.py", SACHS / "scripts" / "D_callable.py",
              SACHS / "scripts" / "kappa2_callable.py", SACHS / "scripts" / "inputs" / "F_tensor.npy",
              SACHS / "scripts" / "inputs" / "kappa2_grid_lcut0_low_empty.npz",
              SACHS / "callables" / "C_propagator" / "corr_op" / "corr_op_C_callable.py",
              SACHS / "callables" / "C_propagator" / "corr_op" / "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"]
    if mapping:
        inputs += [Path(mapping["source_file"]), Path(mapping["response_module"])]
    if input_mode == "corrected-nphi512":
        inputs += [NP512 / "low_permclosed_r4_np512.npz",
                   REPO / "reproduce" / "provenance" / "vertex_rebuild_20260906.md",
                   REPO / "reorg_2026-09" / "nphi512_manifest_2026-09-06.txt"]
    if contract:
        inputs += [GUARD, RESTRICTED_LADDER / "query_contract.json"]
        inputs += [Path(contract[key]["path"]) for key in (
            "original_config", "original_geometry", "geometry_file", "low_piece", "framework_expansion_cache")]
    for spec in specs:
        if mapping:
            endpoints = []
            for redshift in spec["nominal_z"]:
                matches = np.flatnonzero(np.isclose(mapping["z"], redshift, rtol=0, atol=1e-12))
                if len(matches) != 1:
                    raise ValueError(f"Source map has no unique endpoint for z={redshift}")
                endpoints.append(float(mapping["lambda_mpc"][int(matches[0])]))
            spec["target_lambdas_mpc"] = endpoints
        else:
            spec["target_lambdas_mpc"] = list(spec["historical_lambdas_mpc"])
        config = yaml.safe_load(Path(spec["source_config"]).read_text())
        if max(spec["target_lambdas_mpc"]) > float(config["propagators"]["t_max"]):
            raise ValueError(f"{spec['id']} endpoints exceed the configured propagator interval")
        spec["n_gauss"] = int(config["sweep"]["n_gauss"])
        spec["n_jobs"] = n_jobs
        spec["orders"] = ",".join(str(x) for x in config["expand"]["orders"])
        spec["n_components"] = int(config["system"]["field"]["n_components"])
        spec["historical_component_pairs"] = config["sweep"]["component_pairs"]
        if "query_contract" in spec:
            positions = config["sweep"]["positions_grid"]
            if not np.array_equal(positions["x"], [contract["x"]]) or not np.array_equal(
                    np.asarray(positions["y"])[contract["selected_indices"]], contract["y"]):
                raise ValueError(f"Historical ladder angles differ from the query contract: {spec['id']}")
        spec["angle_count"] = (len(spec["selected_angle_indices"]) if "query_contract" in spec
                               else len(config["sweep"]["positions_grid"]["y"]))
        spec["component_pairs"] = (config["sweep"]["component_pairs"]
                                   if spec["group"] == "main" else [[0, 0]])
        spec["component_selection"] = (
            "All six historical observable pairs are needed by the main figures."
            if spec["group"] == "main" else
            "Only convergence (0,0) is plotted or quoted for this target. Internal fields retain all three components."
        )
        inputs += [Path(spec[key]) for key in ("source_config", "record_config", "table", "callable", "historical_product")]
        if "historical_table" in spec:
            inputs.append(Path(spec["historical_table"]))
            inputs.append(Path(spec["input_low_piece"]))
    return {
        "schema_version": 1, "framework": framework(), "source_mode": mode,
        "input_mode": input_mode,
        "ladder_geometry": ladder_geometry,
        "source_mapping": mapping, "targets": specs,
        "inputs": [file_record(path) for path in sorted(set(inputs))],
        "scope": "Historical source mode retains recorded endpoints and response. True-redshift mode replaces them with the explicit source map. K tables are selected separately by input_mode, with recorded callable variants and sampling geometry retained. Output-pair selection does not change the three internal field components.",
        "input_update": ("Explicit corrected-nphi512 selection. The historical framework-only control remains a separate archived run."
                         if input_mode == "corrected-nphi512" else
                         "The September n_phi=512 pair-phase K-table correction is not selected."),
        "prerequisites_for_publication": [
            "Compare the historical-endpoint replay against archived arrays before interpreting the version effect.",
            "True-redshift main and ladder plots require O0 and FF at the same new source endpoint and response convention.",
            "True-redshift multi-z plots require matching O0 and FF for every changed source endpoint.",
            "The multi-z cached ff field contains historical placeholders. Actual FF comes from multiz_ff_real_all5.npz and cannot be replaced with the placeholder.",
            "The figure-6 Monte Carlo comparison requires matching source coordinates and driving input before its overlay can validate a changed prediction.",
            "Multi-z and ladder outputs contain only the convergence pair. Their plot adapters must read that pair directly without using a loader that expects shear fields or inventing zero-valued missing channels.",
            "Re-render derived angular spectra and quoted ratios from the new arrays, and deploy only through reproduce/deploy.py.",
        ],
    }


def variant_module():
    spec = importlib.util.spec_from_file_location("r1_fk_variant", REBUILD / "run_fk_variant.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def validate_config_paths(config: dict, work: Path) -> None:
    for entry in config["output"]:
        if not Path(entry["path"]).resolve().is_relative_to(work):
            raise RuntimeError("Generated output escaped the new batch directory")
    for section in ("expand", "propagators"):
        cache = config.get(section, {}).get("cache_path")
        if cache and not Path(cache).resolve().is_relative_to(work):
            raise RuntimeError("Generated cache escaped the new batch directory")
    system = config["system"]
    paths = [system["linear"]["R_time_module"], system["noise"]["kappa2"]["module"],
             config["propagators"]["c_closed_form_module"]]
    paths += [item["coupling_path"] for item in system.get("vertices", [])]
    paths += [item["coupling_module"] for item in system.get("nonlocal_vertices", [])]
    for path in paths:
        if not Path(path).is_absolute() or not Path(path).is_file():
            raise FileNotFoundError(f"Generated config has an unresolved input: {path}")


def require_tables_in_place(plan: dict) -> None:
    """Refuse targets whose vertex table is archived.

    A generated callable opens its table by path at run time, so the archive
    resolver cannot serve it. Checked before anything is written.
    """
    absent = sorted({target["table"] for target in plan["targets"]
                     if not Path(target["table"]).is_file()})
    if absent:
        raise FileNotFoundError(
            "These vertex tables are archived and must be restored in place "
            "before this batch is prepared or run: " + ", ".join(absent))


def prepare(plan: dict, work: Path) -> None:
    work = work.resolve()
    require_tables_in_place(plan)
    work.mkdir(parents=True, exist_ok=False)
    helper = variant_module()
    generated = []
    for target in plan["targets"]:
        run_dir = work / target["id"]
        callable_path, config_path = helper.materialise_variant(
            Path(target["table"]), run_dir,
            t_final=target["target_lambdas_mpc"], n_jobs=target["n_jobs"],
            base_config=Path(target["source_config"]), orders=target["orders"],
            callable_src=Path(target["callable"]),
        )
        config = yaml.safe_load(config_path.read_text())
        mapping = plan["source_mapping"]
        if mapping:
            config["system"]["linear"]["R_time_module"] = mapping["response_module"]
            config["system"]["linear"]["R_time_attr"] = mapping["response_attr"]
            # The old grid extends past the affine horizon of the matched
            # response. Only the interval enclosing all requested sources is
            # needed, with the original one-Mpc spacing retained.
            config["propagators"]["t_max"] = float(math.ceil(max(target["target_lambdas_mpc"])))
        config["expand"]["n_jobs"] = 1
        config["sweep"]["n_gauss"] = target["n_gauss"]
        config["sweep"]["component_pairs"] = target["component_pairs"]
        if "query_contract" in target:
            contract_path = Path(target["query_contract"])
            contract = guard_module().validate_contract(contract_path)
            guard_module().validate_table(Path(target["table"]), contract)
            config["sweep"]["positions_grid"] = {"x": [contract["x"]], "y": contract["y"]}
            # Import the reviewed guard by explicit path and bind input hashes
            # into this fresh callable, preserving the historical implementation.
            with callable_path.open("a") as stream:
                stream.write("\n# Enforce the restricted ladder contract before any lookup.\n")
                stream.write("import importlib.util as _r1_importlib\n")
                stream.write(f"_r1_spec = _r1_importlib.spec_from_file_location('r1_ladder_guard', {str(GUARD)!r})\n")
                stream.write("_r1_guard = _r1_importlib.module_from_spec(_r1_spec)\n")
                stream.write("_r1_spec.loader.exec_module(_r1_guard)\n")
                stream.write(f"_r1_guard.install_guard(globals(), {str(contract_path)!r}, {digest(contract_path)!r}, {digest(Path(target['table']))!r})\n")
        validate_config_paths(config, work)
        config_path.write_text(yaml.safe_dump(config, sort_keys=False))
        shutil.copyfile(target["record_config"], run_dir / "historical_config.yaml")
        target["prepared_config"] = str(config_path)
        target["result"] = str(run_dir / target["output_name"])
        generated.extend([file_record(config_path), file_record(callable_path), file_record(run_dir / "historical_config.yaml")])
    plan["generated_inputs"] = generated
    plan["work_directory"] = str(work)
    dump_new(work / "prepared.json", plan)
    print(f"Prepared {len(plan['targets'])} targets without running any fold: {work}")


def check_other_folds() -> None:
    process_list = subprocess.check_output(["ps", "-axo", "pid=,command="], text=True)
    matches = [line.strip() for line in process_list.splitlines()
               if any(name in line for name in ("run_FF_single.py", "run_fk_variant.py", "run_ff_one_lambda.py"))]
    if matches:
        raise RuntimeError("Another fold or fold driver is active. Do not overlap folds:\n" + "\n".join(matches))


def run_batch(work: Path) -> None:
    work = work.resolve()
    plan = json.loads((work / "prepared.json").read_text())
    if framework() != plan["framework"]:
        raise RuntimeError("Framework or interpreter changed after preparation")
    os.environ["CANOES_ROOT"] = str(verify_canoes())
    for entry in plan["inputs"] + plan["generated_inputs"]:
        if file_record(Path(entry["path"])) != entry:
            raise RuntimeError(f"Input changed after preparation: {entry['path']}")
    require_tables_in_place(plan)
    # A shared advisory lock prevents two suite batches from running together.
    # The process check also catches historical drivers that do not use the lock.
    with (HERE / ".replay_suite.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another replay suite holds the fold lock") from exc
        for target in plan["targets"]:
            result = Path(target["result"])
            run_dir = result.parent
            if (run_dir / "completed.json").exists():
                completed = json.loads((run_dir / "completed.json").read_text())
                if file_record(result) != completed["result"]:
                    raise RuntimeError(f"Completed result changed: {result}")
                print(f"Already completed: {target['id']}", flush=True)
                continue
            if (run_dir / "started.json").exists() or result.exists():
                raise RuntimeError(f"Refusing to replace an existing or interrupted run: {run_dir}. Prepare a fresh batch.")
            check_other_folds()
            dump_new(run_dir / "started.json", {"start_unix": time.time(), "pid": os.getpid()})
            # The public variant entry point calls the same run() API as the
            # historical folds, but this batch already materialised its config.
            helper = variant_module()
            started = time.monotonic()
            for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
                os.environ[variable] = "1"
            helper.run(Path(target["prepared_config"]), result, target["n_gauss"], target["orders"])
            dump_new(run_dir / "completed.json", {
                "elapsed_seconds": time.monotonic() - started, "result": file_record(result),
                "framework": plan["framework"], "source_mode": plan["source_mode"],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("plan", "prepare"):
        child = sub.add_parser(command)
        child.add_argument("--source-mode", choices=("historical", "true-redshift"), required=True)
        child.add_argument("--input-mode", choices=("historical", "corrected-nphi512"), default="historical")
        child.add_argument("--ladder-geometry", choices=("historical", "contract52"), default="historical")
        child.add_argument("--source-map", type=Path)
        child.add_argument("--targets", nargs="+", default=["all"])
        child.add_argument("--n-jobs", type=int, choices=(4, 6), default=4)
        if command == "prepare":
            child.add_argument("--work-dir", type=Path, required=True)
        else:
            child.add_argument("--manifest", type=Path, help="New file only, otherwise print JSON")
    run = sub.add_parser("run")
    run.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "run":
        run_batch(args.work_dir)
    else:
        plan = build_plan(args.source_mode, args.source_map, args.targets, args.n_jobs,
                          args.input_mode, args.ladder_geometry)
        if args.command == "prepare":
            prepare(plan, args.work_dir)
        elif args.manifest:
            dump_new(args.manifest, plan)
        else:
            print(json.dumps(plan, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
