"""Shared matplotlib styling and primitives for the STF_lensing paper figures.

The four paper plot scripts (``plot_paper_analysis{1,3,4}.py`` and
``plot_paper_lcut.py``) share the same visual language:

* A single thick semi-transparent gray reference line plotted FIRST as
  background -- the "Total" or "PyCCL" curve.  The line is sign-aware:
  solid for positive segments, dashed for negative segments, on log-log
  ``|value|``.
* Component curves plotted as markers only (no connecting line),
  filled for positive values and hollow for negative.

This module centralises the mpl rcParams, the color palette, and the
two drawing primitives so the four figures stay visually consistent.
"""
from __future__ import annotations

import matplotlib as mpl
import numpy as np

# ----- visual constants -----
GAMMA_X_MAX_ARCMIN = 2000.0  # paper-figure x-axis clip

LINE_BG_COLOR = "0.30"
LINE_BG_LW = 3.4
LINE_BG_ALPHA = 0.32
SIGN_MARKER_SIZE = 6.5  # gray markers on top of the |value| line
SIGN_MARKER_ALPHA = 0.55

MARKER_SIZE = 5.2  # coloured component markers
MARKER_LW_FILLED = 0.6
MARKER_LW_HOLLOW = 1.3

# Components are plotted in the order they appear in the legend; this
# palette is colour-blind safe (Tableau 10 / Wong) and intentionally
# avoids gray which is reserved for the background reference line.
PALETTE = (
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
    "#56B4E9",  # sky blue
)
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")


# ----- rcParams -----
def apply_rcparams() -> None:
    """Apply the paper-figure mpl rcParams.  Idempotent."""
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 13.5,
        "axes.titlesize": 13,
        "axes.titleweight": "normal",
        "axes.linewidth": 0.8,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.fontsize": 10.5,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "0.82",
        "legend.fancybox": False,
        "legend.borderpad": 0.5,
        "legend.borderaxespad": 0.7,
        "legend.handlelength": 2.0,
        "legend.handletextpad": 0.6,
        "legend.columnspacing": 1.2,
        "lines.linewidth": 1.5,
        "lines.markersize": MARKER_SIZE,
        "mathtext.fontset": "stix",
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.4,
        "grid.color": "0.80",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
    })


# ----- drawing primitives -----
def plot_signed_line(
    ax,
    gamma,
    arr,
    *,
    color=LINE_BG_COLOR,
    lw=LINE_BG_LW,
    alpha=LINE_BG_ALPHA,
    marker="o",
    sign_marker_size=SIGN_MARKER_SIZE,
    sign_marker_alpha=SIGN_MARKER_ALPHA,
    label=None,
    zorder_line=1,
    zorder_markers=2,
    faint_outside=None,
    faint_alpha=0.28,
):
    """Continuous ``|arr|`` line plus same-color sign markers.

    The line is a single continuous curve on log-log ``|arr|`` (no gaps
    at sign transitions).  Sign is encoded by overlaid markers at every
    data point: filled = positive, hollow = negative.  The line carries
    the legend handle; the sign markers do not produce extra legend
    entries.
    """
    gamma = np.asarray(gamma, dtype=float)
    arr = np.asarray(arr, dtype=float)
    abs_arr = np.abs(arr)

    # ``faint_outside=(lo, hi)`` draws the curve outside that window dashed and
    # translucent: the values are shown but are not the ones being quoted.
    # Used for the FK channel, whose multipole sum stops converging beyond
    # about a degree.  Default None reproduces the previous single-style call
    # exactly.
    keep = np.ones(gamma.size, dtype=bool)
    if faint_outside is not None:
        lo, hi = faint_outside
        if lo is not None:
            keep &= gamma >= lo
        if hi is not None:
            keep &= gamma <= hi
        ax.loglog(
            gamma, abs_arr, linestyle=(0, (3, 2)), color=color,
            lw=lw, alpha=alpha * faint_alpha, solid_capstyle="round",
            zorder=zorder_line,
        )
        g_main = np.where(keep, gamma, np.nan)
        ax.loglog(
            g_main, abs_arr, linestyle="-", color=color,
            lw=lw, alpha=alpha, solid_capstyle="round",
            label=label, zorder=zorder_line,
        )
    else:
        ax.loglog(
            gamma, abs_arr, linestyle="-", color=color,
            lw=lw, alpha=alpha, solid_capstyle="round",
            label=label, zorder=zorder_line,
        )
    pos = (arr > 0) & keep
    neg = (arr < 0) & keep
    if pos.any():
        ax.loglog(
            gamma[pos], abs_arr[pos],
            marker=marker, linestyle="none", color=color,
            mfc=color, mec=color, mew=0.6,
            alpha=sign_marker_alpha, markersize=sign_marker_size,
            zorder=zorder_markers,
        )
    if neg.any():
        ax.loglog(
            gamma[neg], abs_arr[neg],
            marker=marker, linestyle="none", color=color,
            mfc="white", mec=color, mew=1.2,
            alpha=sign_marker_alpha, markersize=sign_marker_size,
            zorder=zorder_markers,
        )


def plot_signed_markers(
    ax,
    gamma,
    arr,
    *,
    color,
    marker,
    label=None,
    size=MARKER_SIZE,
    alpha=0.92,
    zorder=3,
    faint_outside=None,
    faint_alpha=0.25,
):
    """Sign-aware ``|arr|`` markers (no connecting line).

    Positive values: filled marker.  Negative values: hollow marker
    (white face, coloured edge).  Legend label is attached to whichever
    of the two flavours appears first in the data.
    """
    gamma = np.asarray(gamma, dtype=float)
    arr = np.asarray(arr, dtype=float)
    # ``faint_outside=(lo, hi)`` renders points outside that window at reduced
    # opacity: shown, but not the values being quoted. Default None is the
    # previous behaviour exactly.
    keep = np.ones(gamma.size, dtype=bool)
    if faint_outside is not None:
        lo, hi = faint_outside
        if lo is not None:
            keep &= gamma >= lo
        if hi is not None:
            keep &= gamma <= hi
    for sel, a, lbl in ((keep, alpha, label),
                        (~keep, alpha * faint_alpha, None)):
        if not sel.any():
            continue
        pos = (arr > 0) & sel
        neg = (arr < 0) & sel
        used = False
        if pos.any():
            ax.loglog(
                gamma[pos], arr[pos],
                marker=marker, linestyle="none", color=color,
                mfc=color, mec=color, mew=MARKER_LW_FILLED,
                markersize=size, alpha=a,
                label=lbl, zorder=zorder,
            )
            used = True
        if neg.any():
            ax.loglog(
                gamma[neg], -arr[neg],
                marker=marker, linestyle="none", color=color,
                mfc="white", mec=color, mew=MARKER_LW_HOLLOW,
                markersize=size, alpha=a,
                label=None if used else lbl, zorder=zorder,
            )


def annotate_sign_legend(ax, *, loc=(0.98, 0.98), fontsize=None):
    """Standard 'filled marker: positive,  hollow: negative' note.

    Font size defaults to the active ``legend.fontsize`` rcParam so the
    annotation matches the legend across all paper figures.
    """
    if fontsize is None:
        fontsize = mpl.rcParams["legend.fontsize"]
    ax.text(
        loc[0], loc[1],
        "filled markers: $+$    hollow markers: $-$",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=fontsize, color="0.25",
        bbox=dict(
            facecolor="white", edgecolor="0.78",
            boxstyle="round,pad=0.35", alpha=0.92,
        ),
    )


def clip_gamma_axis(ax, gamma_min, gamma_max=GAMMA_X_MAX_ARCMIN, pad=0.05):
    """Set log-x limits with a small margin both sides; clip max at gamma_max."""
    lo = max(gamma_min * (1.0 - pad), 1e-3)
    hi = gamma_max * (1.0 + pad)
    ax.set_xlim(lo, hi)
