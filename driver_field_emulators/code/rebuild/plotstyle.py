"""Shared plotting style for the note's figures.

Deliberately independent of the paper's `_plot_style.py`, which lives in an
analysis directory and pulls in that package's imports. The note is a
standalone document, so its figures carry their own style.

Sign handling follows the paper's convention: the quantities plotted change
sign, so curves are drawn as |value| on log axes with filled markers where
the value is positive and open markers where it is negative. A plain log
plot would silently drop half the curve.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PALETTE = ("#1f4e79", "#c1121f", "#2a9d8f", "#e09f3e", "#6a4c93", "#495057")

RC = {
    "figure.dpi": 130,
    "savefig.dpi": 200,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "legend.framealpha": 0.9,
    "lines.linewidth": 1.6,
}


def apply() -> None:
    plt.rcParams.update(RC)


def signed_loglog(ax, x, y, label=None, color=None, marker="o",
                  markersize=4.2, linestyle="-", floor=None):
    """Plot |y| against x on log axes, marking the sign.

    Filled markers are y > 0, open markers y < 0. Values whose magnitude is
    below `floor` are dropped rather than drawn, because below the
    cancellation floor the sign is numerical noise and drawing it invites a
    reader to interpret it.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    magnitude = np.abs(y)
    keep = np.isfinite(magnitude) & (magnitude > 0)
    if floor is not None:
        keep &= magnitude >= floor
    ax.plot(x[keep], magnitude[keep], linestyle=linestyle, color=color,
            label=label, zorder=3)
    positive = keep & (y > 0)
    negative = keep & (y < 0)
    ax.plot(x[positive], magnitude[positive], marker, color=color,
            markersize=markersize, linestyle="none", zorder=4)
    ax.plot(x[negative], magnitude[negative], marker, color=color,
            markersize=markersize, linestyle="none", markerfacecolor="none",
            markeredgewidth=1.1, zorder=4)


def sign_note(ax, loc="lower left") -> None:
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", color="0.35", linestyle="none",
               markersize=4.5, label="positive"),
        Line2D([], [], marker="o", color="0.35", linestyle="none",
               markersize=4.5, markerfacecolor="none", label="negative"),
    ]
    legend = ax.legend(handles=handles, loc=loc, fontsize=8.5,
                       framealpha=0.85, handletextpad=0.4)
    ax.add_artist(legend)
    return legend


OBSERVABLE_TITLES = {
    "xi_kappa": r"$\xi_\kappa(\gamma)=\langle\kappa\kappa\rangle$",
    "xi_plus": r"$\xi_+(\gamma)=\langle\gamma_+\gamma_+\rangle+\langle\gamma_\times\gamma_\times\rangle$",
    "xi_minus": r"$\xi_-(\gamma)=\langle\gamma_+\gamma_+\rangle-\langle\gamma_\times\gamma_\times\rangle$",
    "xi_kappa_gamma_t": r"$\xi_{\kappa\gamma_t}(\gamma)$",
}
