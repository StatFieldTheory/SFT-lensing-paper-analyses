"""Check restricted ladder preparation and actual callable parity, without folds."""

from __future__ import annotations

import importlib.util
import itertools
import json
from pathlib import Path
import tempfile

import numpy as np
import yaml

import replay_suite as suite
from restricted_kappa3_guard import validate_contract, validate_table


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def assert_rejected(fn, *args) -> None:
    try:
        fn(*args)
    except (ValueError, KeyError):
        return
    raise AssertionError("Unsupported query or incomplete table was accepted")


def main() -> None:
    contract = validate_contract(suite.RESTRICTED_LADDER / "query_contract.json")
    batch = suite.HERE / "true_redshift_fk_gl24_corrected_tree52"
    plan = json.loads((batch / "prepared.json").read_text())
    source = (suite.PLAIN / "equal_time_limber_kappa3_callable.py").read_text()
    marker = 'TABLE_PATH = _HERE / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"'
    with np.load(contract["original_geometry"]["path"], allow_pickle=False) as geometry:
        full_geometry = np.asarray(geometry["cosine_triples"])
    x, ys = np.asarray(contract["x"]), np.asarray(contract["y"])
    placements = np.asarray([[([x, y][i]) for i in ids] for y in ys
                             for ids in itertools.product((0, 1), repeat=3)])
    results = []
    with tempfile.TemporaryDirectory(prefix="r1_ladder_guard_") as temporary:
        temp = Path(temporary)
        for target in plan["targets"]:
            cutoff = target["id"].rsplit("_", 1)[1]
            table = Path(target["table"])
            validate_table(table, contract)
            config = yaml.safe_load(Path(target["prepared_config"]).read_text())
            assert config["system"]["field"]["n_components"] == 3
            assert config["sweep"]["component_pairs"] == [[0, 0]]
            assert config["sweep"]["positions_grid"] == {"x": [contract["x"]], "y": contract["y"]}
            for entry in plan["inputs"] + plan["generated_inputs"]:
                assert suite.file_record(Path(entry["path"])) == entry
            full_path = suite.NP512 / "ladder" / f"table_tree_np512_cut{cutoff}.npz"
            with np.load(full_path, allow_pickle=False) as data:
                arrays = dict(data)
            row_map = {tuple(row): i for i, row in enumerate(arrays["cosine_triples"])}
            rows = [row_map[tuple(row)] for row in full_geometry]
            for key, value in list(arrays.items()):
                if key.startswith("zeta_") or key == "cosine_triples":
                    arrays[key] = value[rows]
            reference_table = temp / f"full_r4_{cutoff}.npz"
            np.savez(reference_table, **arrays)
            reference_source = temp / f"full_callable_{cutoff}.py"
            reference_source.write_text(source.replace(marker, f"TABLE_PATH = Path({str(reference_table)!r})"))
            reference = load_module(reference_source, f"r1_full_{cutoff}")
            guarded = load_module(batch / target["id"] / f"callable_{target['id']}.py", f"r1_guard_{cutoff}")
            shells = arrays["lambda_shells_Mpc"]
            times = np.sort(np.concatenate([shells, (shells[:-1] + shells[1:]) / 2]))
            sample_count = 0
            for t in times:
                n_arr = np.moveaxis(placements, 0, 1)
                t_arr = np.full(n_arr.shape[:2], t)
                assert np.array_equal(guarded.coupling_fn_batch(n_arr, t_arr),
                                      reference.coupling_fn_batch(n_arr, t_arr))
                for n in placements:
                    assert np.array_equal(guarded.coupling_fn(n, [t] * 3),
                                          reference.coupling_fn(n, [t] * 3))
                    sample_count += 1
            t = float(shells[5])
            invalid = np.asarray([x, x, ys[0]]).copy()
            invalid[2, 0] = np.nextafter(invalid[2, 0], np.inf)
            for positions in (invalid, np.asarray([x, ys[0], ys[1]]), np.full((3, 3), np.nan)):
                assert_rejected(guarded.coupling_fn, positions, [t] * 3)
                assert_rejected(guarded.coupling_fn_batch, positions[:, None], np.full((3, 1), t))
            with np.load(table, allow_pickle=False) as data:
                damaged = dict(data)
            del damaged["zeta_Dmod"]
            damaged_path = temp / "missing_channel.npz"
            np.savez(damaged_path, **damaged)
            assert_rejected(validate_table, damaged_path, contract)
            results.append({"target": target["id"], "scalar_samples": sample_count,
                            "batch_samples": sample_count, "max_absolute_difference": 0.0,
                            "invalid_scalar_and_batch_queries_rejected": 6,
                            "missing_channel_rejected": True})
    report = {"status": "PASS", "targets": results,
              "contract": suite.file_record(suite.RESTRICTED_LADDER / "query_contract.json"),
              "guard": suite.file_record(suite.GUARD), "prepared": suite.file_record(batch / "prepared.json")}
    path = suite.RESTRICTED_LADDER / "guard_verification.json"
    suite.dump_new(path, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
