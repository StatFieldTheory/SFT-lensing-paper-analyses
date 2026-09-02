"""Convergence figure: several folded FK runs against the last one.

Serves two of the three analyses, because both have the same shape: vary
one numerical parameter, fold, and ask whether the answer has stopped
moving. Used for the grid-density study (r1, r2, r4) and for the multipole
cutoff study.

Top row: the curves themselves. Bottom row: each curve divided by the
reference, which is the last one given. A converged parameter shows a
bottom row flat at one over the range that carries signal.

Usage::

    python fig_convergence.py --out fig.pdf --observable xi_kappa \
        --curve "r1=..._xi.npz" --curve "r2=..._xi.npz" --curve "r4=..._xi.npz"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle
from observables import load_observables

FLOOR = 1e-15


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curve", action="append", required=True,
                    help="label=path/to/xi.npz; the LAST one is the reference")
    ap.add_argument("--observable", nargs="+", default=["xi_kappa"])
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--title", default=None)
    ap.add_argument("--ratio-label", default="ratio to reference")
    ap.add_argument("--gamma-max", type=float, default=700.0,
                    help="ratios are only shown below this separation, where "
                         "the curves are above the cancellation floor")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    curves = []
    for spec in args.curve:
        label, _, path = spec.partition("=")
        gamma, obs = load_observables(Path(path), order=args.order)
        curves.append((label, gamma, obs))
    reference = curves[-1]

    plotstyle.apply()
    import matplotlib.pyplot as plt

    n = len(args.observable)
    fig, axes = plt.subplots(2, n, figsize=(5.4 * n, 7.2), sharex=True,
                             squeeze=False,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
    for c, name in enumerate(args.observable):
        top, bottom = axes[0][c], axes[1][c]
        for i, (label, gamma, obs) in enumerate(curves):
            if obs[name] is None:
                continue
            plotstyle.signed_loglog(top, gamma, obs[name], label=label,
                                    color=plotstyle.PALETTE[i % len(plotstyle.PALETTE)],
                                    marker="oDs^v"[i % 5], floor=FLOOR)
        top.set_xscale("log")
        top.set_yscale("log")
        top.set_title(args.title or plotstyle.OBSERVABLE_TITLES.get(name, name))
        top.set_ylabel("magnitude")
        top.legend(loc="best")

        ref_values = reference[2][name]
        for i, (label, gamma, obs) in enumerate(curves[:-1]):
            if obs[name] is None or ref_values is None:
                continue
            keep = (np.abs(ref_values) > FLOOR) & (gamma <= args.gamma_max)
            bottom.plot(gamma[keep], obs[name][keep] / ref_values[keep],
                        color=plotstyle.PALETTE[i % len(plotstyle.PALETTE)],
                        marker="oDs^v"[i % 5], markersize=4.0,
                        label=f"{label} / {reference[0]}")
        bottom.axhline(1.0, color="0.4", linewidth=0.9, linestyle="--")
        bottom.set_xscale("log")
        bottom.set_xlabel(r"separation $\gamma$ [arcmin]")
        bottom.set_ylabel(args.ratio_label)
        bottom.legend(loc="best", fontsize=8.5)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[fig] -> {args.out}")

    name = args.observable[0]
    ref_values = reference[2][name]
    print(f"\n{name}: values and ratio to {reference[0]}")
    header = "gamma[']  " + "  ".join(f"{lab:>13}" for lab, _, _ in curves)
    print(header + "   " + "  ".join(f"{lab}/{reference[0]:>6}"
                                     for lab, _, _ in curves[:-1]))
    gamma = curves[0][1]
    for j, g in enumerate(gamma):
        if g > args.gamma_max:
            continue
        row = f"{g:8.2f}  " + "  ".join(f"{c[2][name][j]:+13.5e}" for c in curves)
        row += "   " + "  ".join(
            f"{c[2][name][j] / ref_values[j]:8.4f}" if abs(ref_values[j]) > FLOOR
            else "     n/a" for c in curves[:-1])
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
