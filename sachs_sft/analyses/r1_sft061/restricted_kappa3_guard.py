"""Validate and enforce the explicit 11-angle ladder query contract.

The underlying historical sorting callable remains responsible for numerical
interpolation and tensor reconstruction. This module validates its complete
query neighborhoods and prevents use outside the retained external positions.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import yaml

CHANNELS = ("TTT", "TTP", "TPP", "PPP", "Bmod", "Dmod")


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _check_record(record: dict) -> Path:
    path = Path(record["path"])
    if sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise ValueError(f"Contract input changed: {path}")
    return path


def _queries(placements: np.ndarray, batch: bool) -> np.ndarray:
    unit = placements / np.linalg.norm(placements, axis=-1, keepdims=True)
    if batch:
        cosines = [np.einsum("ij,ij->i", unit[:, a], unit[:, b])
                   for a, b in ((0, 1), (1, 2), (2, 0))]
        return np.sort(np.clip(np.stack(cosines, axis=-1), -1, 1), axis=1)[:, ::-1]
    return np.asarray([sorted([float(np.clip(np.dot(row[a], row[b]), -1, 1))
                               for a, b in ((0, 1), (1, 2), (2, 0))], reverse=True)
                       for row in unit])


def validate_contract(path: Path, expected_hash: str | None = None) -> dict:
    """Recheck provenance, complete duplicate groups and exact KD neighbors."""
    if expected_hash is not None and sha256(path) != expected_hash:
        raise ValueError("Restricted query contract changed after preparation")
    contract = json.loads(path.read_text())
    for key in ("original_config", "original_geometry", "geometry_file",
                "low_piece", "framework_expansion_cache"):
        _check_record(contract[key])
    config = yaml.safe_load(Path(contract["original_config"]["path"]).read_text())
    positions = config["sweep"]["positions_grid"]
    x, ys = np.asarray(contract["x"]), np.asarray(contract["y"])
    if not np.array_equal(positions["x"], [x]) or not np.array_equal(
            np.asarray(positions["y"])[contract["selected_indices"]], ys):
        raise ValueError("Contract positions differ from the original angle selection")
    with np.load(contract["original_geometry"]["path"], allow_pickle=False) as data:
        full = np.asarray(data["cosine_triples"])
    with np.load(contract["geometry_file"]["path"], allow_pickle=False) as data:
        selected = np.asarray(data["cosine_triples"])
        rows = np.asarray(data["original_r4_row_indices"])
    if not np.array_equal(rows, contract["raw_row_indices"]) or not np.array_equal(full[rows], selected):
        raise ValueError("Restricted raw rows differ from the recorded full geometry")
    canonical = np.round(np.sort(full, axis=1)[:, ::-1], 8)
    full_unique, inverse = np.unique(canonical, axis=0, return_inverse=True)
    groups = np.unique(inverse[rows])
    complete_rows = np.flatnonzero(np.isin(inverse, groups))
    if not np.array_equal(rows, complete_rows):
        raise ValueError("Restricted geometry omits raw members of a duplicate group")
    reduced_unique = np.unique(np.round(np.sort(selected, axis=1)[:, ::-1], 8), axis=0)
    placements = np.asarray([[([x, y][i]) for i in ids] for y in ys
                             for ids in itertools.product((0, 1), repeat=3)])
    scalar_queries = _queries(placements, batch=False)
    if not np.array_equal(np.unique(scalar_queries, axis=0), contract["sorted_queries"]):
        raise ValueError("Scalar cosine queries differ from the recorded contract")
    full_tree, reduced_tree = cKDTree(full_unique), cKDTree(reduced_unique)
    for queries in (scalar_queries, _queries(placements, batch=True)):
        distances, indices = full_tree.query(queries, k=4)
        reduced_distances, reduced_indices = reduced_tree.query(queries, k=4)
        if not np.array_equal(distances, reduced_distances) or not np.array_equal(
                full_unique[indices], reduced_unique[reduced_indices]):
            raise ValueError("Restricted geometry changes a scalar or batch query neighborhood")
    recorded_queries = np.asarray(contract["sorted_queries"])
    distances, indices = full_tree.query(recorded_queries, k=4)
    if not np.array_equal(distances, contract["neighbor_distances"]) or not np.array_equal(
            full_unique[indices], contract["neighbor_coordinates"]):
        raise ValueError("Full geometry no longer matches the recorded neighbors")
    return contract


def validate_table(table: Path, contract: dict, expected_hash: str | None = None) -> None:
    """Require every retained row, radial shell and tensor input channel."""
    if expected_hash is not None and sha256(table) != expected_hash:
        raise ValueError("Restricted vertex table changed after preparation")
    with np.load(contract["geometry_file"]["path"], allow_pickle=False) as geometry, \
            np.load(contract["low_piece"]["path"], allow_pickle=False) as low, \
            np.load(table, allow_pickle=False) as data:
        if not np.array_equal(data["cosine_triples"], geometry["cosine_triples"]):
            raise ValueError("Restricted table row order differs from the contract")
        if not np.array_equal(data["lambda_shells_Mpc"], low["lambda_shells"]):
            raise ValueError("Restricted table radial shells differ from the corrected LOW input")
        shape = (contract["selected_raw_rows"], contract["radial_shell_count"])
        for channel in CHANNELS:
            values = data[f"zeta_{channel}"]
            if values.shape != shape or not np.all(np.isfinite(values)):
                raise ValueError(f"Restricted table has missing or invalid {channel} values")


def install_guard(namespace: dict, contract_path: str, contract_hash: str,
                  table_hash: str) -> None:
    """Wrap both public entry points before the historical lookup can run."""
    contract = validate_contract(Path(contract_path), contract_hash)
    validate_table(Path(namespace["TABLE_PATH"]), contract, table_hash)
    directions = np.asarray([contract["x"], *contract["y"]], dtype=float)
    original_scalar = namespace["coupling_fn"]
    original_batch = namespace["coupling_fn_batch"]

    def check_positions(arr: np.ndarray) -> None:
        # The framework passes the external config positions verbatim.
        # Exact identity prevents even nearby unsupported angles being used.
        matches = np.all(arr[..., None, :] == directions, axis=-1)
        if not np.all(np.sum(matches, axis=-1) == 1):
            raise ValueError("K query contains a direction outside the explicit ladder contract")
        ids = np.argmax(matches, axis=-1)
        largest = np.max(ids, axis=0)
        if not np.all((ids == 0) | (ids == largest[None, :])):
            raise ValueError("K query combines different ladder angle pairs")

    def scalar(n_list, t_list):
        arr = np.asarray(n_list, dtype=float)
        if arr.shape != (3, 3):
            raise ValueError(f"Expected three three-dimensional vertex directions, got {arr.shape}")
        check_positions(arr[:, None, :])
        return original_scalar(n_list, t_list)

    def batch(n_arr, t_arr):
        arr = np.asarray(n_arr, dtype=float)
        if arr.ndim != 3 or arr.shape[0] != 3 or arr.shape[2] != 3:
            raise ValueError(f"Expected directions shaped (3, samples, 3), got {arr.shape}")
        check_positions(arr)
        return original_batch(n_arr, t_arr)

    scalar.__dict__.update(original_scalar.__dict__)
    batch.__dict__.update(original_batch.__dict__)
    namespace["coupling_fn"] = scalar
    namespace["coupling_fn_batch"] = batch
