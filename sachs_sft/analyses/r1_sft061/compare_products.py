"""Compare trusted local SFT-Wick products on identical component/angle grids."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def load_product(path: Path) -> dict:
    """Load locally generated object arrays and verify their scientific row keys."""
    with np.load(path, allow_pickle=True) as saved:
        x = np.stack(saved["x"]).astype(float)
        y = np.stack(saved["y"]).astype(float)
        x /= np.linalg.norm(x, axis=1)[:, None]
        y /= np.linalg.norm(y, axis=1)[:, None]
        angle = np.arccos(np.clip(np.sum(x * y, axis=1), -1, 1)) * 10800 / np.pi
        values = np.real_if_close(np.asarray(saved["value"]))
        if np.iscomplexobj(values) or not np.all(np.isfinite(values)):
            raise ValueError(f"Nonreal or nonfinite values in {path}")
        result = {}
        for a, b, order, endpoint, theta, value in zip(
            saved["a"], saved["b"], saved["order"], saved["t_final"], angle, values,
            strict=True,
        ):
            key = (int(a), int(b), int(order), float(endpoint), float(theta))
            if key in result:
                raise ValueError(f"Duplicate scientific row key in {path}: {key}")
            if not np.all(np.isfinite(key)):
                raise ValueError(f"Nonfinite row coordinates in {path}")
            result[key] = float(value)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, required=True)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Comparison output exists. Choose a fresh path.")
    old, new = load_product(args.old), load_product(args.new)
    if set(old) != set(new):
        raise ValueError(f"Scientific row grids differ: {len(old)} old, {len(new)} new rows")
    groups = sorted({key[:4] for key in old})
    records = []
    for group in groups:
        keys = sorted(key for key in old if key[:4] == group)
        before = np.array([old[key] for key in keys])
        after = np.array([new[key] for key in keys])
        theta = np.array([key[4] for key in keys])
        delta = after - before
        peak = float(np.max(np.abs(before)))
        keep = np.abs(before) >= 0.001 * peak if peak else np.zeros(len(keys), bool)
        records.append({
            "a": group[0], "b": group[1], "order": group[2],
            "lambda_source_mpc": group[3], "samples": len(keys),
            "array_equal": bool(np.array_equal(before, after)),
            "max_abs_difference": float(np.max(np.abs(delta))),
            "max_abs_difference_over_old_peak_percent":
                float(100 * np.max(np.abs(delta)) / peak) if peak else None,
            "max_abs_relative_percent_above_0p1_percent_old_peak":
                float(100 * np.max(np.abs(delta[keep] / before[keep]))) if keep.any() else None,
            "relative_mask_excluded_fraction": float(1 - keep.mean()),
            "theta_arcmin": theta.tolist(), "old": before.tolist(), "new": after.tolist(),
        })
    report = {
        "old_path": str(args.old.resolve()), "new_path": str(args.new.resolve()),
        "old_sha256": hashlib.sha256(args.old.read_bytes()).hexdigest(),
        "new_sha256": hashlib.sha256(args.new.read_bytes()).hexdigest(),
        "row_grid_exactly_equal": True, "row_count": len(old), "components": records,
        "interpretation": "Matched-input historical-versus-current framework replay. Input-table changes are excluded.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for row in records:
        print({key: value for key, value in row.items() if key not in {"theta_arcmin", "old", "new"}})


if __name__ == "__main__":
    main()
