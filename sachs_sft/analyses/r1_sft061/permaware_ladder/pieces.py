"""Prepare exact-row inputs and merge independently built BiHalofit pieces.

This driver never starts a vertex integration. Its jobs JSON is dispatched
separately, serially, by the parent task. Existing artifacts are read only.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
R1 = HERE.parent
sys.path.insert(0, str(R1))
import replay_suite as replay  # noqa: E402

sys.path.insert(0, str(replay.REBUILD))
from _common import dump_meta, grid_fingerprint, load_meta  # noqa: E402

FIELDS = ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP", "bmod", "dmod")
FULL = replay.NP512 / "table_permclosed_np512.npz"
PERM_SOURCE = replay.PERM / "perm_aware_kappa3_callable.py"
OLD = R1 / "corrected_kappa3_ladder52"
PERMS = tuple(itertools.permutations(range(3)))


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save_npz(path: Path, arrays: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        np.savez(stream, **arrays)


def read_piece(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        arrays = dict(data)
    geometry, lam = arrays["cosine_triples"], arrays["lambda_shells"]
    if load_meta(arrays["fingerprint"]) != grid_fingerprint(geometry, lam):
        raise ValueError(f"Piece fingerprint does not match actual coordinates: {path}")
    for key in FIELDS:
        if arrays[key].shape != (len(geometry), 16) or not np.all(np.isfinite(arrays[key])):
            raise ValueError(f"Invalid piece channel {key}: {path}")
    if lam.shape != (16,) or not np.all(np.isfinite(lam)) or not np.all(np.diff(lam) > 0):
        raise ValueError(f"Invalid radial grid: {path}")
    return arrays


def subset(source: Path, rows: np.ndarray, output: Path) -> dict:
    arrays = read_piece(source)
    for key in (*FIELDS, "cosine_triples"):
        arrays[key] = arrays[key][rows]
    arrays["fingerprint"] = dump_meta(grid_fingerprint(arrays["cosine_triples"], arrays["lambda_shells"]))
    arrays["row_selection_provenance"] = dump_meta({"source": replay.file_record(source),
                                                  "rows": rows.tolist(), "operation": "exact raw-row selection"})
    save_npz(output, arrays)
    return replay.file_record(output)


def placements_and_queries(x: np.ndarray, ys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    placements = np.asarray([[([x, y][i]) for i in ids] for y in ys
                             for ids in itertools.product((0, 1), repeat=3)])
    nn = placements / np.linalg.norm(placements, axis=-1, keepdims=True)
    triples = np.stack([np.clip(np.einsum("ij,ij->i", nn[:, a], nn[:, b]), -1, 1)
                        for a, b in ((0, 1), (1, 2), (2, 0))], axis=-1)
    loader = load_module(PERM_SOURCE, "r1_geometry_loader")
    return placements, np.unique(np.concatenate([loader._permute_triple(triples, p) for p in PERMS]), axis=0)


def initialize() -> None:
    """Create geometry, exact LOW/tree subsets and eight missing-row jobs."""
    contract_path = HERE / "contract.json"
    outputs = [contract_path, HERE / "triples34.npz", HERE / "triples_missing22.npz",
               HERE / "low34.npz", HERE / "low_missing22.npz", HERE / "bihalofit_missing_jobs.json"]
    if any(path.exists() for path in outputs):
        raise FileExistsError("Permaware preparation already exists, refusing overwrite")
    old_contract = json.loads((OLD / "query_contract.json").read_text())
    x, ys = np.asarray(old_contract["x"]), np.asarray(old_contract["y"])
    placements, queries = placements_and_queries(x, ys)
    with np.load(FULL, allow_pickle=False) as full:
        full_triples = full["cosine_triples"]
        shells = full["lambda_shells_Mpc"]
    key, inverse = np.unique(np.round(full_triples, 12), axis=0, return_inverse=True)
    distances, nearest = cKDTree(key).query(queries, k=1)
    if not np.all(distances < 1e-10):
        raise ValueError("Not every main-loader query uses its exact-nearest shortcut")
    rows = np.flatnonzero(np.isin(inverse, np.unique(nearest)))
    with np.load(OLD / "triples_ladder52.npz", allow_pickle=False) as data:
        old_triples = data["cosine_triples"]
    old_keys = {tuple(row) for row in old_triples}
    missing = np.asarray([row for row in rows if tuple(full_triples[row]) not in old_keys])
    if len(rows) != 34 or len(missing) != 22:
        raise ValueError(f"Geometry changed: wanted 34 required and 22 missing rows, got {len(rows)}, {len(missing)}")
    save_npz(HERE / "triples34.npz", {"cosine_triples": full_triples[rows], "full_row_indices": rows})
    save_npz(HERE / "triples_missing22.npz", {"cosine_triples": full_triples[missing], "full_row_indices": missing})
    low = replay.NP512 / "low_permclosed_r4_np512.npz"
    subset(low, rows, HERE / "low34.npz")
    subset(low, missing, HERE / "low_missing22.npz")
    tree_records = []
    for source in sorted(replay.NP512.glob("band_tree_permclosed_*_np512.npz")):
        tree_records.append(subset(source, rows, HERE / "pieces" / source.name.replace("permclosed", "34")))
    if len(tree_records) != 8:
        raise ValueError("Expected eight corrected tree bands")
    old_jobs = json.loads((OLD / "bihalofit_jobs.json").read_text())
    jobs = []
    for old in old_jobs:
        job = {**old, "argv": list(old["argv"])}
        output = HERE / "missing_pieces" / Path(old["out"]).name.replace("r4", "missing22")
        job["out"] = str(output)
        job["log"] = str(HERE / "logs" / (output.stem + ".log"))
        job["name"] = output.stem
        for flag, value in (("--triples-npz", HERE / "triples_missing22.npz"), ("--out", output)):
            job["argv"][job["argv"].index(flag) + 1] = str(value)
        jobs.append(job)
    for path in (HERE / "missing_pieces", HERE / "logs"):
        path.mkdir(parents=True, exist_ok=True)
    replay.dump_new(HERE / "bihalofit_missing_jobs.json", jobs)
    contract = {
        "kind": "permutation-aware exact-query restriction", "source": replay.file_record(FULL),
        "loader": replay.file_record(PERM_SOURCE), "source_map": replay.file_record(R1 / "source_map.json"),
        "original_config": old_contract["original_config"],
        "framework_expansion_cache": old_contract["framework_expansion_cache"],
        "selected_indices": old_contract["selected_indices"], "x": x.tolist(), "y": ys.tolist(),
        "queries": queries.tolist(), "nearest_keys": key[nearest].tolist(),
        "nearest_distances": distances.tolist(), "raw_full_indices": rows.tolist(),
        "missing_full_indices": missing.tolist(), "key_decimals": 12, "exact_threshold": 1e-10,
        "row_count": 34, "lambda_shells": shells.tolist(),
        "geometry": replay.file_record(HERE / "triples34.npz"),
        "missing_geometry": replay.file_record(HERE / "triples_missing22.npz"),
        "low": replay.file_record(HERE / "low34.npz"), "tree_pieces": tree_records,
        "required_channels": ["zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP", "zeta_Bmod", "zeta_Dmod"],
        "n_ell": 96, "n_phi": 512, "low_Nmax": 384,
        "scope": "Exact external positions only, no nonexact nearest-neighbor blending is allowed",
    }
    replay.dump_new(contract_path, contract)
    print(f"Prepared {len(placements)} external placements, 34 full rows, 22 missing rows, eight serial build jobs")


def merge_bihalofit() -> None:
    """Combine the 12 previously built rows and 22 missing rows per band."""
    contract = json.loads((HERE / "contract.json").read_text())
    with np.load(HERE / "triples34.npz", allow_pickle=False) as data:
        wanted = data["cosine_triples"]
    jobs = json.loads((HERE / "bihalofit_missing_jobs.json").read_text())
    sources = []
    for job in jobs:
        missing_path = Path(job["out"])
        old_path = OLD / "pieces" / missing_path.name.replace("missing22", "r4")
        output = HERE / "pieces" / missing_path.name.replace("missing22", "34")
        if output.exists():
            raise FileExistsError(output)
        old, missing = read_piece(old_path), read_piece(missing_path)
        for key in ("lambda_shells", "chi_shells", "z_shells", "ell_max", "cosmo_meta"):
            if not np.array_equal(old[key], missing[key]):
                raise ValueError(f"Band physical metadata differs: {key}, {missing_path}")
        old_build, new_build = load_meta(old["build"]), load_meta(missing["build"])
        if {k: v for k, v in old_build.items() if k != "triples_npz"} != {
                k: v for k, v in new_build.items() if k != "triples_npz"}:
            raise ValueError("Cannot merge builds with different quadrature/model/window provenance")
        if (new_build["n_ell"], new_build["n_phi"], new_build["b_model"]) != (96, 512, "bihalofit"):
            raise ValueError("Unexpected BiHalofit build precision")
        pools = [{tuple(row): i for i, row in enumerate(data["cosine_triples"])} for data in (old, missing)]
        assembled = dict(old)
        for key in FIELDS:
            selected = []
            for row in wanted:
                matches = [(data, pool[tuple(row)]) for data, pool in zip((old, missing), pools, strict=True)
                           if tuple(row) in pool]
                if len(matches) != 1:
                    raise ValueError("Each required ordered row must be supplied exactly once")
                data, index = matches[0]
                selected.append(data[key][index])
            assembled[key] = np.asarray(selected)
        assembled["cosine_triples"] = wanted
        assembled["fingerprint"] = dump_meta(grid_fingerprint(wanted, assembled["lambda_shells"]))
        assembled["build"] = dump_meta({**new_build, "triples_npz": "triples34.npz"})
        provenance = {"old52": replay.file_record(old_path), "missing22": replay.file_record(missing_path),
                      "contract": replay.file_record(HERE / "contract.json"),
                      "operation": "exact ordered-row union, 12 existing plus 22 new, no averaging"}
        assembled["merge_provenance"] = dump_meta(provenance)
        sources += list(provenance[key] for key in ("old52", "missing22"))
        save_npz(output, assembled)
    replay.dump_new(HERE / "bihalofit_merge_manifest.json", {
        "inputs": sources, "contract": replay.file_record(HERE / "contract.json"),
        "outputs": [replay.file_record(p) for p in sorted((HERE / "pieces").glob("band_bihalofit_34_*.npz"))],
        "lambda_shells": contract["lambda_shells"],
    })
    print("Merged eight BiHalofit bands, all artifacts preserved")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("initialize", "merge-bihalofit"))
    args = parser.parse_args()
    {"initialize": initialize, "merge-bihalofit": merge_bihalofit}[args.command]()


if __name__ == "__main__":
    main()
