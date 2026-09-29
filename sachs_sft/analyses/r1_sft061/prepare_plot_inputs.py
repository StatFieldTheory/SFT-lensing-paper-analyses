"""Adapt validated numerical products to the existing manuscript plot formats.

These operations evaluate no new physical model. Historical products are read
only, and output paths must be fresh. Every output records its input hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent


def record(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def fresh(path: Path) -> None:
    if path.suffix != ".npz":
        raise ValueError("Output must have an explicit .npz suffix")
    if path.exists() or path.with_suffix(".json").exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)


def write(path: Path, payload: dict, provenance: dict) -> None:
    fresh(path)
    np.savez(path, **payload, metadata_json=json.dumps(provenance, sort_keys=True))
    provenance["output"] = record(path)
    path.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")


def directions(config: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions = yaml.safe_load(config.read_text())["sweep"]["positions_grid"]
    x = np.asarray(positions["x"], float)
    y = np.asarray(positions["y"], float)
    if x.shape != (1, 3) or y.ndim != 2 or y.shape[1] != 3:
        raise ValueError("Expected one reference direction and a set of three-vector directions")
    unit_x = x[0] / np.linalg.norm(x[0])
    unit_y = y / np.linalg.norm(y, axis=1)[:, None]
    gamma = np.arccos(np.clip(unit_y @ unit_x, -1, 1)) * 10800 / np.pi
    order = np.argsort(gamma)
    if np.any(np.diff(gamma[order]) <= 0):
        raise ValueError("Directions do not give distinct angular separations")
    return x[0], y[order], gamma[order]


def order0(args: argparse.Namespace) -> None:
    x, y, gamma = directions(args.config)
    with np.load(args.source, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        if metadata["observable_operator"] != "full":
            raise ValueError("The manuscript baseline requires the full Ricci source")
        if not np.allclose(gamma, data["theta_arcmin"], rtol=0, atol=1e-7):
            raise ValueError("Source and plot angular grids differ")
        lam = float(metadata["source_lambda_mpc"])
        expected = json.loads(args.source_map.read_text())
        match = np.flatnonzero(np.asarray(expected["z"]) == 5)
        if match.size != 1 or not np.isclose(lam, expected["lambda_mpc"][int(match[0])], rtol=0, atol=1e-10):
            raise ValueError("Order-0 is not on the reviewed z=5 source plane")
        kk, plus, minus, cross = [np.asarray(data[k], float) for k in ("kk", "xip", "xim", "kgt")]
    for array in (kk, plus, minus, cross):
        if array.shape != gamma.shape or not np.all(np.isfinite(array)):
            raise ValueError("Invalid source correlation function")
    # Invert the same observable combinations used by the existing plotters.
    # The two parity-odd correlations vanish for the scalar reference model.
    groups = {(0, 0): kk, (0, 1): -cross, (0, 2): np.zeros_like(kk),
              (1, 1): (plus + minus) / 2, (1, 2): np.zeros_like(kk),
              (2, 2): (plus - minus) / 2}
    if not np.allclose(groups[(1, 1)] + groups[(2, 2)], plus, rtol=2e-15, atol=1e-20):
        raise ValueError("Spin-trace round trip failed")
    if not np.allclose(groups[(1, 1)] - groups[(2, 2)], minus, rtol=2e-10, atol=1e-19):
        raise ValueError("Spin-difference round trip failed")
    n = len(gamma)
    payload = {
        "x": np.tile(x, (n * len(groups), 1)), "y": np.tile(y, (len(groups), 1)),
        "a": np.repeat([pair[0] for pair in groups], n),
        "b": np.repeat([pair[1] for pair in groups], n),
        "order": np.zeros(n * len(groups), dtype=int),
        "t_final": np.full(n * len(groups), lam),
        "value": np.concatenate(list(groups.values())),
    }
    provenance = {"kind": "Order-0 plot-format adapter", "source": record(args.source),
                  "config": record(args.config), "source_map": record(args.source_map),
                  "lambda_source_mpc": lam, "z_source": 5,
                  "parity_policy": "The scalar Gaussian reference has zero kappa-cross and plus-cross correlations.",
                  "calculation": "Continuous outer-response projection in the source product. This file is not a framework rerun."}
    write(args.out, payload, provenance)


def zeta(args: argparse.Namespace) -> None:
    x, y, gamma = directions(args.config)
    keep = gamma <= 85
    x = x / np.linalg.norm(x)
    cosine = (y[keep] / np.linalg.norm(y[keep], axis=1)[:, None]) @ x
    queries = np.column_stack([np.ones_like(cosine), cosine, cosine])
    indices, mismatches = [], []
    with np.load(args.source, allow_pickle=False) as data:
        if int(data["ell_max"]) != 15360:
            raise ValueError("Expected the reviewed ell_max=15360 corrected table")
        triples = data["cosine_triples"]
        for query in queries:
            difference = np.max(np.abs(triples - query), axis=1)
            matches = np.flatnonzero(difference <= 1e-12)
            if matches.size != 1:
                raise ValueError(f"Expected exactly one stored collapsed row for {query}")
            index = int(matches[0])
            indices.append(index)
            mismatches.append(float(difference[index]))
        payload = {"gamma_arcmin": gamma[keep], "lambda_phys": data["lambda_shells_Mpc"],
                   "z_shells": data["z_shells"]}
        for channel in ("TTT", "TTP", "TPP", "PPP", "Bmod", "Dmod"):
            values = np.asarray(data[f"zeta_{channel}"][indices], float)
            if values.shape != (keep.sum(), len(payload["lambda_phys"])) or not np.all(np.isfinite(values)):
                raise ValueError(f"Invalid stored cumulant channel: {channel}")
            payload[f"full_{channel}"] = values
    provenance = {"kind": "Stored corrected-table collapsed slices", "source": record(args.source),
                  "config": record(args.config), "ell_max": 15360,
                  "table_rows": indices, "max_cosine_difference": max(mismatches),
                  "geometry": "(1, cos(gamma), cos(gamma)), unchanged stored values, no interpolation",
                  "angle_policy": "Original 40-point observation grid restricted to gamma <= 85 arcmin",
                  "scope": "Driving-field shell statistics, not a source-plane lensing observable"}
    write(args.out, payload, provenance)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("order0")
    command.add_argument("--source", type=Path, required=True)
    command.add_argument("--config", type=Path, required=True)
    command.add_argument("--source-map", type=Path, default=HERE / "source_map.json")
    command.add_argument("--out", type=Path, required=True)
    command = sub.add_parser("zeta")
    command.add_argument("--source", type=Path, required=True)
    command.add_argument("--config", type=Path, required=True)
    command.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    {"order0": order0, "zeta": zeta}[args.command](args)


if __name__ == "__main__":
    main()
