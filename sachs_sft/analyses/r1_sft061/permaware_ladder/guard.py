"""Restrict the unchanged permutation-aware loader to verified exact queries."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def install(namespace, contract_path, contract_hash, table_hash):
    if digest(contract_path) != contract_hash or digest(namespace["TABLE_PATH"]) != table_hash:
        raise ValueError("Permutation-aware ladder input changed")
    contract = json.loads(Path(contract_path).read_text())
    for key in ("source", "loader", "geometry", "source_map"):
        record = contract[key]
        if digest(record["path"]) != record["sha256"]:
            raise ValueError(f"Contract provenance changed: {key}")
    with np.load(contract["geometry"]["path"]) as geometry, np.load(namespace["TABLE_PATH"]) as table:
        if not np.array_equal(table["cosine_triples"], geometry["cosine_triples"]):
            raise ValueError("Wrong ordered geometry")
        if not np.array_equal(table["lambda_shells_Mpc"], contract["lambda_shells"]):
            raise ValueError("Wrong affine shells")
        for key in contract["required_channels"]:
            if table[key].shape != (34, 16) or not np.all(np.isfinite(table[key])):
                raise ValueError(f"Missing or nonfinite channel: {key}")
    directions = np.asarray([contract["x"], *contract["y"]])
    original_batch = namespace["coupling_fn_batch"]
    original_locate = namespace["_locate"]

    def locate(queries, lam):
        result = original_locate(queries, lam)
        if not np.all(result["exact"]):
            raise ValueError("Restricted table cannot blend a nonexact query")
        return result

    def batch(n_arr, t_arr):
        arr = np.asarray(n_arr, dtype=float)
        times = np.asarray(t_arr, dtype=float)
        if arr.ndim != 3 or arr.shape[0] != 3 or arr.shape[2] != 3:
            raise ValueError("Expected directions with shape (3, samples, 3)")
        if times.shape != arr.shape[:2] or not np.all(np.isfinite(times)):
            raise ValueError("Invalid affine coordinates")
        matches = np.all(arr[..., None, :] == directions, axis=-1)
        if not np.all(np.sum(matches, axis=-1) == 1):
            raise ValueError("Direction outside the explicit external grid")
        ids = np.argmax(matches, axis=-1)
        largest = np.max(ids, axis=0)
        if not np.all((ids == 0) | (ids == largest[None, :])):
            raise ValueError("Mixed external angular separations")
        return original_batch(n_arr, t_arr)

    batch.__dict__.update(original_batch.__dict__)
    namespace["_locate"] = locate
    namespace["coupling_fn_batch"] = batch
    # The existing scalar entry point calls coupling_fn_batch through its
    # module globals and therefore receives both checks without a new formula.
