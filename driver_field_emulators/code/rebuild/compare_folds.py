"""Compare folded FK results across table variants, in all four observables.

Produces the figure the revision needs: the same four observables the paper
plots, with the deployed FK curve and one or more corrected curves on the
same axes, and Order-0 for scale. Everything except the kappa3 table is held
fixed, so a difference between curves is the table and nothing else.

Usage::

    python compare_folds.py --order 2 \
        --curve deployed  <deployed FK xi.npz> \
        --curve corrected <variant xi.npz> \
        --o0 <O0 xi.npz> --out ../../figures/fk_grid_comparison.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from observables import OBSERVABLES, load_observables

TITLES = {
    "xi_kappa": r"$\xi_\kappa$",
    "xi_plus": r"$\xi_+$",
    "xi_minus": r"$\xi_-$",
    "xi_kappa_gamma_t": r"$\xi_{\kappa\gamma_t}$",
}
#: Below this the curves are at the cancellation floor and ratios are noise.
FLOOR = 1e-18


def signed_loglog(ax, gamma, value, color, label, marker):
    """Plot |value| on log axes, filled markers positive, open negative.

    The FK curves change sign inside the plotted range, so a linear or
    symlog axis compresses three decades of decay into nothing. This is the
    convention the paper's own figures use.
    """
    magnitude = np.abs(value)
    ok = magnitude > FLOOR
    ax.plot(gamma[ok], magnitude[ok], "-", color=color, lw=1.4, label=label)
    positive = ok & (value > 0)
    negative = ok & (value < 0)
    ax.plot(gamma[positive], magnitude[positive], marker, color=color, ms=4.5,
            mfc=color, ls="none")
    ax.plot(gamma[negative], magnitude[negative], marker, color=color, ms=4.5,
            mfc="none", ls="none")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curve", nargs=2, action="append", metavar=("LABEL", "PATH"),
                    required=True)
    ap.add_argument("--o0", type=Path, default=None)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--t-final", type=float, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    curves = {}
    for label, path in args.curve:
        gamma, obs = load_observables(Path(path), args.order, args.t_final)
        curves[label] = (gamma, obs)
        print(f"[cmp] {label:24s} {Path(path).name}")

    o0 = None
    if args.o0 is not None:
        g0, o0 = load_observables(args.o0, 0, args.t_final)

    labels = list(curves)
    reference = labels[0]
    gamma_ref = curves[reference][0]

    for name in OBSERVABLES:
        print(f"\n=== {name} (order {args.order}) ===")
        header = "  gamma[']" + "".join(f"{l:>15s}" for l in labels)
        if o0 is not None:
            header += f"{'O0':>15s}"
        print(header + "".join(f"{l + '/' + reference:>16s}" for l in labels[1:]))
        for i, g in enumerate(gamma_ref):
            if not (g <= 1.0 or 4.0 < g < 6.0 or 20 < g < 24 or 40 < g < 46
                    or 110 < g < 120 or 180 < g < 190 or 900 < g < 1000):
                continue
            row = f"  {g:8.2f}"
            for label in labels:
                value = curves[label][1][name]
                row += f"{value[i]:+15.4e}" if value is not None else f"{'-':>15s}"
            if o0 is not None and o0[name] is not None:
                row += f"{o0[name][i]:+15.4e}"
            base = curves[reference][1][name]
            for label in labels[1:]:
                other = curves[label][1][name]
                if base is None or other is None or abs(base[i]) < FLOOR:
                    row += f"{'-':>16s}"
                else:
                    row += f"{other[i] / base[i]:+16.4f}"
            print(row)

    if args.out is None:
        return 0

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = ["#444444", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    markers = ["o", "s", "^", "D", "v"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6), sharex=True)
    for ax, name in zip(axes.ravel(), OBSERVABLES):
        for k, label in enumerate(labels):
            gamma, obs = curves[label]
            if obs[name] is None:
                continue
            signed_loglog(ax, gamma, obs[name], palette[k % len(palette)],
                          label, markers[k % len(markers)])
        if o0 is not None and o0[name] is not None:
            ax.plot(g0, np.abs(o0[name]), ":", color="0.55", lw=1.6,
                    label="Order 0" if name == "xi_kappa" else None)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(TITLES[name])
        ax.grid(alpha=0.25, which="both", lw=0.4)
        ax.set_xlabel(r"$\gamma$ [arcmin]")
    axes[0, 0].legend(fontsize=9, frameon=False)
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    print(f"\n[cmp] -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
