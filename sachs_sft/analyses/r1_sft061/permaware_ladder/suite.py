"""Assemble and prepare fresh permutation-aware cutoff batches, without folding."""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay_suite as replay

CUTS = (960, 1920, 3840, 7680, 15360)


def prepare(model):
    replay.verify_canoes()
    base = json.loads((HERE.parent / "true_redshift_fk_gl24_corrected_main/prepared.json").read_text())
    if replay.framework() != base["framework"]:
        raise ValueError("Framework differs from the completed main fold")
    for record in base["inputs"] + base["generated_inputs"]:
        if replay.file_record(Path(record["path"])) != record:
            raise ValueError(f"Main fold input changed: {record['path']}")
    contract_path = HERE / "contract.json"
    contract = json.loads(contract_path.read_text())
    work = HERE / f"folds_{model}"
    if work.exists():
        raise FileExistsError(work)
    tables = []
    for cutoff in CUTS:
        table = HERE / f"table_{model}_cut{cutoff}.npz"
        if table.exists():
            raise FileExistsError(table)
        subprocess.run([
            "/Users/zzhang/projects/angular_statistics/canoes/.venv/bin/python",
            str(replay.REBUILD / "assemble.py"), "--low", str(HERE / "low34.npz"),
            "--bands", str(HERE / f"pieces/band_{model}_34_*_np512.npz"),
            "--cutoff", str(cutoff), "--out", str(table),
        ], check=True)
        tables.append(table)
    work.mkdir()
    plan = copy.deepcopy(base)
    plan.update(ladder_geometry="permaware34", work_directory=str(work), targets=[], generated_inputs=[])
    extras = [Path(__file__), HERE / "pieces.py", HERE / "guard.py", contract_path,
              HERE / "low34.npz", HERE / "triples34.npz", *tables]
    extras += sorted((HERE / "pieces").glob(f"band_{model}_34_*_np512.npz"))
    plan["inputs"] += [replay.file_record(p) for p in extras]
    plan["scope"] = "Same main permutation-aware loader, GL24, true z=5, exact 34-row query restriction"
    for cutoff, table in zip(CUTS, tables, strict=True):
        target = copy.deepcopy(base["targets"][0])
        target.update(id=f"{model}_{cutoff}", group="ladder", table=str(table),
                      output_name=f"xi_{model}_{cutoff}.npz", component_pairs=[[0, 0]],
                      angle_count=11, selected_angle_indices=contract["selected_indices"])
        run = work / target["id"]
        callable_path, cfg_path = replay.variant_module().materialise_variant(
            table, run, t_final=target["target_lambdas_mpc"], n_jobs=4,
            base_config=Path(target["source_config"]), orders=target["orders"],
            callable_src=Path(contract["loader"]["path"]))
        with callable_path.open("a") as stream:
            stream.write("\nimport importlib.util as _r1_util\n")
            stream.write(f"_r1_spec = _r1_util.spec_from_file_location('r1_perm_guard', {str(HERE / 'guard.py')!r})\n")
            stream.write("_r1_guard = _r1_util.module_from_spec(_r1_spec)\n_r1_spec.loader.exec_module(_r1_guard)\n")
            stream.write(f"_r1_guard.install(globals(), {str(contract_path)!r}, {replay.digest(contract_path)!r}, {replay.digest(table)!r})\n")
        cfg = yaml.safe_load(cfg_path.read_text())
        mapping = base["source_mapping"]
        cfg["system"]["linear"].update(R_time_module=mapping["response_module"], R_time_attr=mapping["response_attr"])
        cfg["propagators"]["t_max"] = float(math.ceil(max(target["target_lambdas_mpc"])))
        cfg["expand"]["n_jobs"] = 1
        cfg["sweep"].update(n_gauss=24, component_pairs=[[0, 0]],
                              positions_grid={"x": [contract["x"]], "y": contract["y"]})
        replay.validate_config_paths(cfg, work)
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
        target.update(prepared_config=str(cfg_path), result=str(run / target["output_name"]))
        plan["targets"].append(target)
        plan["generated_inputs"] += [replay.file_record(cfg_path), replay.file_record(callable_path)]
    replay.dump_new(work / "prepared.json", plan)
    print(f"Prepared {work}, no fold started", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", choices=("tree", "bihalofit"))
    prepare(parser.parse_args().model)
