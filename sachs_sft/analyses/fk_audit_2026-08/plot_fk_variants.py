"""Figures for the FK variant study, in the paper's own visual language.

Reuses ``analysis3/_plot_style.py`` so these panels can be read against the
draft's Order-0/FF/FK figure without a style mismatch: the same signed
log-log convention (filled markers positive, hollow negative), the same
palette, the same gamma clipping.

Two figures:

* ``fk_variants_kappa.pdf`` - xi_FK(gamma) for the kk pair across the four
  table variants, with the Order-0 curve as the grey reference, plus the
  ratio panel that makes the crossover question explicit.
* ``fk_variant_effects.pdf`` - the three effects isolated as ratios to the
  grid-fixed tree baseline: grid, cutoff, bispectrum model.

Outputs land in ``figures``. The draft's own
``figures/`` directory is not touched.

Usage (from ``sachs_sft/analyses/fk_audit_2026-08``)::

    SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    $SFT plot_fk_variants.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSES = (
    HERE.parents[1] / "analyses"
    / "analysis3"
)
sys.path.insert(0, str(ANALYSES))

from _plot_style import (  # noqa: E402
    GAMMA_X_MAX_ARCMIN,
    MARKERS,
    PALETTE,
    annotate_sign_legend,
    apply_rcparams,
    clip_gamma_axis,
    plot_signed_line,
    plot_signed_markers,
)

from compare_fk_variants import load_curve  # noqa: E402

RUNS = (
    HERE.parents[1]
    / "sftwick_outputs" / "2PCF"
)
PRODUCTS = HERE / "products"
FIGURES = HERE / "figures"

O0_NPZ = RUNS / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
VARIANTS = [
    ("Deployed grid (tree, $\\ell_{\\max}=1000$)",
     RUNS / "C_corr_op_K_limber_FK" / "xi_C_corr_op_K_limber_FK.npz"),
    ("Dense grid (tree, $\\ell_{\\max}=1000$)",
     PRODUCTS / "collapsed_dense_tree_cut1000_xi.npz"),
    ("Dense grid (tree, $\\ell_{\\max}=2000$)",
     PRODUCTS / "collapsed_dense_tree_cut2000_xi.npz"),
    ("Dense grid (BiHalofit, $\\ell_{\\max}=1000$)",
     PRODUCTS / "collapsed_dense_bihalofit_cut1000_xi.npz"),
    ("Dense grid (tree, $\\ell_{\\max}=8192$, banded)",
     PRODUCTS / "collapsed_dense_tree_conv8192_xi.npz"),
]


def _load_all():
    gamma, o0 = load_curve(O0_NPZ, order=0)
    curves = []
    for label, path in VARIANTS:
        g, v = load_curve(path, order=2)
        if not np.allclose(g, gamma):
            raise ValueError(f"gamma grid of {path.name} differs from Order-0")
        curves.append((label, v))
    return gamma, o0, curves


def figure_curves(gamma, o0, curves) -> None:
    import matplotlib.pyplot as plt

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.2, 8.6), sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0]},
    )

    plot_signed_line(ax_top, gamma, o0, label="Order-0 (reference)")
    for i, (label, value) in enumerate(curves):
        plot_signed_markers(ax_top, gamma, value, color=PALETTE[i],
                            marker=MARKERS[i], label=label)
    ax_top.set_ylabel(r"$|\xi_{\kappa\kappa}(\gamma)|$")
    ax_top.set_ylim(1e-11, 3e-3)
    annotate_sign_legend(ax_top)
    ax_top.legend(loc="lower left", frameon=False, fontsize=10.5)

    for i, (label, value) in enumerate(curves):
        ax_bot.plot(gamma, np.abs(value / o0), color=PALETTE[i],
                    marker=MARKERS[i], markersize=4.5, linewidth=1.0,
                    markerfacecolor=PALETTE[i], label=label)
    ax_bot.axhline(1.0, color="0.35", linewidth=1.0, linestyle=":")
    ax_bot.text(gamma[1], 1.25, "FK = Order-0", color="0.35", fontsize=10.5)
    ax_bot.set_xscale("log")
    ax_bot.set_yscale("log")
    ax_bot.set_ylim(1e-3, 1e2)
    ax_bot.set_ylabel(r"$|\xi_{\rm FK}\,/\,\xi_{\rm Order\text{-}0}|$")
    ax_bot.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
    clip_gamma_axis(ax_bot, gamma_min=gamma.min(), gamma_max=GAMMA_X_MAX_ARCMIN)

    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(FIGURES / f"fk_variants_kappa{ext}", dpi=180)
    plt.close(fig)
    print(f"-> {FIGURES / 'fk_variants_kappa.pdf'}")


def figure_effects(gamma, curves) -> None:
    """The three effects as ratios to the grid-fixed tree baseline."""
    import matplotlib.pyplot as plt

    labels = [label for label, _ in curves]
    values = {label: value for label, value in curves}
    baseline = values[labels[1]]  # dense tree, cutoff 1000

    effects = [
        ("Interpolation grid (deployed / dense)", values[labels[0]] / baseline),
        (r"Cutoff $\ell_{\max}$: 2000 / 1000", values[labels[2]] / baseline),
        ("Bispectrum: BiHalofit / tree", values[labels[3]] / baseline),
        (r"Cutoff $\ell_{\max}$: 8192 / 1000", values[labels[4]] / baseline),
    ]

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for i, (label, ratio) in enumerate(effects):
        ax.plot(gamma, np.abs(ratio), color=PALETTE[i], marker=MARKERS[i],
                markersize=4.8, linewidth=1.1, label=label)
    ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1e3)
    ax.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
    ax.set_ylabel("ratio to dense tree, $\\ell_{\\max}=1000$")
    ax.legend(loc="upper left", frameon=False, fontsize=10.5)
    clip_gamma_axis(ax, gamma_min=gamma.min(), gamma_max=GAMMA_X_MAX_ARCMIN)
    fig.tight_layout()
    for ext in (".png", ".pdf"):
        fig.savefig(FIGURES / f"fk_variant_effects{ext}", dpi=180)
    plt.close(fig)
    print(f"-> {FIGURES / 'fk_variant_effects.pdf'}")


def main() -> None:
    apply_rcparams()
    gamma, o0, curves = _load_all()
    figure_curves(gamma, o0, curves)
    figure_effects(gamma, curves)


if __name__ == "__main__":
    main()
