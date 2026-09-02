"""Figure: the published FK curve against the interpolation weight.

Two panels. Left, the normalised FK curve and the four-neighbour
inverse-distance weight on the (1, 1, 1) node, drawn together. Right, their
difference, to show that the agreement is not merely visual.

Everything shown comes from the deployed table and the deployed folded
result, so the figure does not depend on any rebuild.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    d = np.load(args.npz)
    gamma, shape, weight = d["gamma_arcmin"], d["fk_shape"], d["corner_weight"]

    plotstyle.apply()
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    left.plot(gamma, shape, "-", color=plotstyle.PALETTE[1], linewidth=2.4,
              label=r"published FK curve, $\xi^{FK}_\kappa(\gamma)/\xi^{FK}_\kappa(0.5')$")
    left.plot(gamma, weight, "o", color=plotstyle.PALETTE[0], markersize=5.0,
              markerfacecolor="none", markeredgewidth=1.3,
              label=r"interpolation weight on the $(1,1,1)$ node")
    left.set_xscale("log")
    left.set_xlabel(r"separation $\gamma$ [arcmin]")
    left.set_ylabel("normalised amplitude")
    left.legend(loc="lower left")
    left.set_title("the plotted angular dependence")

    inside = weight > 1e-6
    right.plot(gamma[inside], np.abs(shape[inside] - weight[inside]), "o-",
               color=plotstyle.PALETTE[0], markersize=4.2)
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlabel(r"separation $\gamma$ [arcmin]")
    right.set_ylabel("absolute difference")
    right.set_title("difference between the two")
    worst = np.max(np.abs(shape[inside] - weight[inside]))
    right.axhline(worst, color=plotstyle.PALETTE[1], linestyle="--", linewidth=1.0)
    right.text(gamma[1], worst * 1.35, f"largest difference {worst:.1e}",
               color=plotstyle.PALETTE[1], fontsize=9.5)

    leaves = gamma[~inside]
    if leaves.size:
        for ax in (left, right):
            ax.axvline(leaves.min(), color="0.5", linestyle=":", linewidth=1.2)
        left.text(leaves.min() * 0.92, 0.22,
                  f"  corner node leaves the\n  four-neighbour set\n  ({leaves.min() / 60:.1f} deg)",
                  fontsize=8.5, color="0.35", ha="right")

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[fig] -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
