"""Compare FK two-point results across kappa3 table variants.

Answers the two questions the vertex-level diagnostics left open:

1. Does the observable inherit the vertex's ell_high_max sensitivity?
   Compare the same bispectrum at cutoff 1000 and 2000. If xi_FK moves like
   the vertex did, the FK amplitude is UV-conditioned and the paper must say
   so; if it is stable, the fold suppresses the UV-sensitive part.
2. How much of the deployed FK curve is the interpolation grid?
   Compare the deployed covgrid16 table against a tree table densified on
   the collapsed (1, c, c) family that the fold actually queries.

Usage (from ``driver_field_emulators/code``)::

    python compare_fk_variants.py --label deployed <deployed_xi.npz> \\
                                  --label dense-tree-1000 <variant_xi.npz> ...
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

ARCMIN = 180.0 * 60.0 / np.pi
#: Component-pair indices of the (Phi_00, Re Psi_0, Im Psi_0) basis.
KK = (0, 0)


def load_curve(path: Path, pair: tuple[int, int] = KK, order: int = 2):
    """Return ``(gamma_arcmin, xi)`` for one component pair and order.

    The sweep result stores the direction vectors as object arrays, so it
    can only be read with ``allow_pickle``. These files are written by the
    local sft-wick run moments earlier in the same pipeline; they are never
    downloaded or user-supplied.
    """
    data = np.load(path, allow_pickle=True)
    select = (data["a"] == pair[0]) & (data["b"] == pair[1]) & (data["order"] == order)
    x = np.stack([np.asarray(v, dtype=float) for v in data["x"][select]])
    y = np.stack([np.asarray(v, dtype=float) for v in data["y"][select]])
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    y /= np.linalg.norm(y, axis=1, keepdims=True)
    cos_gamma = np.clip(np.einsum("ij,ij->i", x, y), -1.0, 1.0)
    gamma = np.arccos(cos_gamma) * ARCMIN
    value = np.asarray(data["value"][select], dtype=float)
    order_index = np.argsort(gamma)
    return gamma[order_index], value[order_index]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", action="append", default=[], required=True)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--order", type=int, default=2)
    args = parser.parse_args()
    if len(args.label) != len(args.paths):
        raise SystemExit("give one --label per path")

    curves = {}
    for label, path in zip(args.label, args.paths):
        gamma, value = load_curve(path, order=args.order)
        curves[label] = (gamma, value)
        print(f"[{label}] {path.name}: {value.size} points, "
              f"xi(min gamma) = {value[0]:+.6e}")

    labels = list(curves)
    reference = labels[0]
    gamma_ref, value_ref = curves[reference]
    print(f"\n=== FK xi ratios against '{reference}' (order {args.order}) ===")
    header = "  gamma[']  " + "".join(f"{lab:>16s}" for lab in labels)
    print(header)
    for i, gamma in enumerate(gamma_ref):
        row = ""
        for label in labels:
            g, v = curves[label]
            j = int(np.argmin(np.abs(g - gamma)))
            if abs(g[j] - gamma) > 0.05 * max(gamma, 1e-6):
                row += f"{'--':>16s}"
            elif label == reference:
                row += f"{v[j]:16.6e}"
            else:
                row += f"{v[j] / value_ref[i]:16.4f}"
        print(f"  {gamma:8.2f}  " + row)

    print("\n(first column is the absolute reference value, the rest are ratios)")


if __name__ == "__main__":
    main()
