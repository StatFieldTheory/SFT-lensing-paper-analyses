#!/usr/bin/env python
"""Clean, self-explanatory appendix validation figure: the cosmic-shear 3PCF
predicted by our SFT driving-field 3-cumulant (zeta_D) vs the independent
external code fastnc (Sugiyama+2024), as a function of triangle opening angle.

Single panel, publication quality. Reads outputs/appendix_fastnc_curve.npz.

Run:  <sft-wick python> make_appendix_fastnc_figure.py [--quantity fi|apex]
Writes (for the chosen quantity):
  figures/appendix_shear3pcf_fastnc.pdf   (paper figure, also copied to <paper>/figures/)
  outputs/appendix_shear3pcf_fastnc.png   (for on-screen review)
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
FIG = os.path.join(HERE, "figures")
PAPER_FIG = "/Users/zzhang/Documents/MyDrafts/STF_lensing/figures"
NPZ = os.path.join(OUT, "appendix_fastnc_curve.npz")

# ---- shared house style (same module the paper's analysis figures use) ----
sys.path.insert(0, HERE)
from _plot_style import PALETTE, apply_rcparams  # noqa: E402

apply_rcparams()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator  # noqa: E402

# This work = blue (PALETTE[0], the same colour used for the SFT/O0 baseline in
# the analysis figures); fastnc reference = vermilion (PALETTE[1]).
C_OURS = PALETTE[0]
C_FNC = PALETTE[1]


def frame_invariant(d, side):
    return d[f"{side}_fi"]            # (n_gamma, n_phi)


def apex_modulus(d, side):
    return np.abs(d[f"{side}_G1"])   # (n_gamma, n_phi)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quantity", choices=["fi", "apex"], default="fi")
    ap.add_argument("--gamma", type=float, default=50.0,
                    help="which scale to plot as the primary curve (arcmin)")
    ap.add_argument("--second-scale", action="store_true",
                    help="overplot the other gamma as a faint second pair")
    ap.add_argument("--no-band", action="store_true",
                    help="omit the x2 agreement band around fastnc")
    args = ap.parse_args()

    d = np.load(NPZ, allow_pickle=True)
    phi = d["phi_deg"]
    gammas = d["gamma_arcmin"]
    ig = int(np.argmin(np.abs(gammas - args.gamma)))

    if args.quantity == "fi":
        ours = frame_invariant(d, "ours")
        fnc = frame_invariant(d, "fnc")
        # Descriptive, symbol-free label: the precise definition (root-sum-square
        # of the four natural components) is given in the caption, so we do NOT
        # introduce the natural-component symbol on the axis.
        ylabel = "cosmic-shear 3PCF amplitude"
        qtag = "frame-invariant 3PCF"
    else:
        ours = apex_modulus(d, "ours")
        fnc = apex_modulus(d, "fnc")
        ylabel = "cosmic-shear 3PCF (apex channel)"
        qtag = "apex natural component"

    fig, ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)

    g = gammas[ig]

    # subtle "within a factor of 2" band around the external reference: makes the
    # order-of-magnitude agreement readable at a glance, no caption required.
    if not args.no_band:
        ax.fill_between(phi, fnc[ig] / 2.0, fnc[ig] * 2.0, color=C_FNC,
                        alpha=0.085, lw=0, zorder=1,
                        label=r"fastnc $\times\,2^{\pm1}$")

    # Two methods distinguished by colour + marker shape (both filled: this is a
    # positive amplitude, so we avoid hollow markers, which mean "negative" in the
    # companion figure's sign convention).
    ax.plot(phi, fnc[ig], color=C_FNC, lw=1.8, ls="--", marker="s", ms=5.2,
            mfc=C_FNC, mec=C_FNC, zorder=3, label="fastnc (tree)")
    ax.plot(phi, ours[ig], color=C_OURS, lw=2.0, ls="-", marker="o", ms=5.2,
            mfc=C_OURS, mec=C_OURS, zorder=4,
            label=r"This work (SFT $\zeta_D$)")

    if args.second_scale:
        ig2 = 1 - ig
        g2 = gammas[ig2]
        ax.plot(phi, fnc[ig2], color=C_FNC, lw=1.3, ls="--", marker="s", ms=3.5,
                mfc="white", mec=C_FNC, mew=1.0, alpha=0.45, zorder=2)
        ax.plot(phi, ours[ig2], color=C_OURS, lw=1.3, ls="-", marker="o", ms=3.5,
                mfc=C_OURS, mec=C_OURS, alpha=0.45, zorder=2,
                label=fr"(faint: $\gamma={g2:.0f}'$)")

    ax.set_yscale("log")
    ax.set_xlabel(r"triangle opening angle  $\varphi$  [deg]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(phi.min() - 3, phi.max() + 3)
    # compact config annotation in lieu of a title (the caption carries the
    # description; gamma is fixed here so it is shown in-axes for a standalone read)
    ax.text(0.035, 0.05, fr"$\gamma={g:.0f}'$,  $z_s=5$",
            transform=ax.transAxes, ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="0.78", alpha=0.92))

    # legend order: This work, fastnc, band (focus first, reference shading last)
    handles, labels = ax.get_legend_handles_labels()
    order_keys = ["This work", "fastnc (tree)", r"fastnc $\times"]
    idx = []
    for key in order_keys:
        for i, lab in enumerate(labels):
            if lab.startswith(key) and i not in idx:
                idx.append(i)
                break
    idx += [i for i in range(len(labels)) if i not in idx]
    ax.legend([handles[i] for i in idx], [labels[i] for i in idx],
              loc="upper right")
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=12))

    ratio = ours[ig] / fnc[ig]
    stem = f"appendix_shear3pcf_fastnc_{args.quantity}"
    png = os.path.join(OUT, stem + ".png")
    pdf = os.path.join(FIG, stem + ".pdf")
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(png, dpi=200)
    fig.savefig(pdf)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    print(f"[{qtag}] gamma={g:.0f}'  ratio ours/fnc range "
          f"[{np.nanmin(ratio):.3f}, {np.nanmax(ratio):.3f}], "
          f"median {np.nanmedian(ratio):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
