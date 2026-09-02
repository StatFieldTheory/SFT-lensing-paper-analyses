"""Selection-rule "mixing hierarchy" ladder (paper Fig.~\\ref{fig: selection rule}).

A compact, single-column line-art schematic in the paper's muted schematic
palette (matching Figs.~1, 3, 4): which driving-field cumulant $\\mathcal{K}^{(n)}$
(left) reaches which observable (right), at which leading order and through how
many nonlinear-propagation vertices $F$.

Visual grammar (two independent axes):
  * colour  -- building block: blue = free Gaussian 2-cumulant (Order 0),
    teal = a pure source vertex $K^{(n)}$ (no $F$), coral = $F$-mixing;
  * line style -- number of $F$-gates: solid 0, dashed 1, dotted 2.

No boxes/blocks: symbols are typeset directly as nodes.

Output:
    figures/selection_rule_schematic.{pdf,png}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
FIGURES = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py

# ----- muted schematic palette (cf. make_screen_basis_spin2_figure.py) -----
INK = "#16202a"
BLUE = "#2f6f9f"      # free / Gaussian (Order 0)
TEAL = "#2a8580"      # pure source vertex K (direct, no F)
CORAL = "#b95f58"     # F-mixing
GRAY = "#68717d"      # captions / vdots / legend

# line style = number of F-gates
LS_DIRECT = "solid"
LS_1F = (0, (5, 2.6))
LS_2F = (0, (1.2, 2.4))
LW = 1.4


def setup_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 360,
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 9.0,
        "axes.linewidth": 0.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _node(ax, x, y, s, *, size=13.0, color=INK):
    ax.text(x, y, s, color=color, fontsize=size, ha="center", va="center",
            zorder=6)


def _arrow(ax, p0, p1, *, color, ls):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=9,
        connectionstyle="arc3,rad=0.0", color=color, lw=LW, ls=ls,
        shrinkA=2.0, shrinkB=2.0, zorder=4,
    ))


def _label(ax, x, y, s, *, color):
    t = ax.text(x, y, s, color=color, fontsize=8.6, ha="center", va="center",
                zorder=7)
    t.set_path_effects([pe.withStroke(linewidth=2.4, foreground="white",
                                      alpha=0.95)])


def render(seed: int, out_stem: Path) -> None:
    del seed
    setup_style()
    # Full-width, wide-and-short band so it places like the Feynman figures.
    fig, ax = plt.subplots(figsize=(8.2, 2.9))
    ax.set_xlim(0, 18)
    ax.set_ylim(0.0, 6.4)
    ax.axis("off")

    # ----- nodes (no boxes) -----
    cum = {"k2": 4.2, "k3": 2.75, "k4": 1.3}
    cum_lab = {"k2": r"$\mathcal{K}^{(2)}$", "k3": r"$\mathcal{K}^{(3)}$",
               "k4": r"$\mathcal{K}^{(4)}$"}
    CX = 2.2
    for key, y in cum.items():
        _node(ax, CX, y, cum_lab[key])
    ax.text(CX, 5.25, "driving cumulants", color=GRAY, fontsize=8.5, ha="center")
    ax.text(CX, 0.35, r"$\vdots$", color=GRAY, fontsize=14, ha="center", va="center")

    obs = {"xi": 3.5, "zeta": 1.55}
    OX = 15.8
    _node(ax, OX, obs["xi"], r"$\xi$", size=14)
    ax.text(OX, obs["xi"] + 0.6, "2PCF", color=GRAY, fontsize=8.0, ha="center")
    _node(ax, OX, obs["zeta"], r"$\zeta$", size=14)
    ax.text(OX, obs["zeta"] + 0.6, "3PCF", color=GRAY, fontsize=8.0, ha="center")
    ax.text(OX, 5.25, "observables", color=GRAY, fontsize=8.5, ha="center")
    ax.text(OX, 0.55, r"$\vdots$", color=GRAY, fontsize=14, ha="center", va="center")

    xs, xe = 3.3, 14.7   # arrow start/end x (clear of the symbols)

    # ----- channels.  number = leading perturbative order. -----
    _arrow(ax, (xs, cum["k2"]), (xe, obs["xi"]), color=BLUE, ls=LS_DIRECT)
    _label(ax, 7.6, 4.18, r"$\mathrm{O}0\,\cdot\,\mathcal{K}^{(2)}$", color=BLUE)

    _arrow(ax, (xs, cum["k3"]), (xe, obs["zeta"]), color=TEAL, ls=LS_DIRECT)
    _label(ax, 7.1, 2.12, r"$\mathrm{O}1\,\cdot\,K^{(3)}$", color=TEAL)

    _arrow(ax, (xs, cum["k3"]), (xe, obs["xi"]), color=CORAL, ls=LS_1F)
    _label(ax, 8.9, 3.42, r"$\mathrm{O}2\,\cdot\,FK^{(3)}$", color=CORAL)

    _arrow(ax, (xs, cum["k4"]), (xe, obs["zeta"]), color=CORAL, ls=LS_1F)
    _label(ax, 10.1, 1.25, r"$\mathrm{O}2\,\cdot\,FK^{(4)}$", color=CORAL)

    _arrow(ax, (xs, cum["k4"]), (xe, obs["xi"]), color=CORAL, ls=LS_2F)
    _label(ax, 11.7, 2.62, r"$\mathrm{O}3\,\cdot\,FFK^{(4)}$", color=CORAL)

    # ----- legend: line style = number of F-gates (single horizontal row, top) -----
    y_leg = 5.9
    for x0, x1, xt, ls, txt in (
            (2.5, 3.5, 3.7, LS_DIRECT, "direct"),
            (7.2, 8.2, 8.4, LS_1F, r"one $F$"),
            (11.8, 12.8, 13.0, LS_2F, r"two $F$")):
        ax.plot([x0, x1], [y_leg, y_leg], color=GRAY, ls=ls, lw=LW,
                solid_capstyle="round", zorder=4)
        ax.text(xt, y_leg, txt, color=GRAY, fontsize=8.0, ha="left", va="center")

    pdf = out_stem.with_suffix(".pdf")
    png = out_stem.with_suffix(".png")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-stem", type=Path,
                        default=FIGURES / "selection_rule_schematic")
    args = parser.parse_args()
    render(args.seed, args.out_stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
