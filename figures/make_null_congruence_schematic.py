from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py
OUT_PDF = FIG_DIR / "null_geodesic_congruence.pdf"
OUT_PNG = FIG_DIR / "null_geodesic_congruence.png"


NAVY = "#1f2d3a"
BLUE = "#6f8fa3"
TEAL = "#2a8580"
MINT = "#d8ece8"
GOLD = "#b9872f"
CORAL = "#b85e4d"
GRAY = "#6f7782"
LIGHT_GRAY = "#d8dde3"
INK = "#17212b"


def center_y(x: np.ndarray | float) -> np.ndarray | float:
    return 0.18 * np.sin(0.72 * (np.asarray(x) - 1.2)) - 0.018 * (np.asarray(x) - 5.0)


def dy_dx(x: float) -> float:
    return 0.18 * 0.72 * np.cos(0.72 * (x - 1.2)) - 0.018


def beam_width(x: np.ndarray) -> np.ndarray:
    t = np.clip((x - 0.55) / 8.85, 0.0, 1.0)
    return 0.04 + 1.22 * t**1.08


def arrow(
    ax: plt.Axes,
    xy0: tuple[float, float],
    xy1: tuple[float, float],
    *,
    color: str,
    lw: float = 1.6,
    ms: float = 12.0,
    linestyle: str = "-",
    alpha: float = 1.0,
    zorder: int = 5,
) -> FancyArrowPatch:
    patch = FancyArrowPatch(
        xy0,
        xy1,
        arrowstyle="-|>",
        mutation_scale=ms,
        linewidth=lw,
        linestyle=linestyle,
        color=color,
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def label(
    ax: plt.Axes,
    x: float,
    y: float,
    s: str,
    *,
    size: float = 10,
    color: str = INK,
    ha: str = "center",
    va: str = "center",
    weight: str | None = None,
    zorder: int = 10,
) -> None:
    text = ax.text(
        x,
        y,
        s,
        fontsize=size,
        color=color,
        ha=ha,
        va=va,
        fontweight=weight,
        zorder=zorder,
    )
    text.set_path_effects([pe.withStroke(linewidth=3.2, foreground="white", alpha=0.92)])


def tangent_frame(x0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    p = np.array([x0, float(center_y(x0))])
    t = np.array([1.0, dy_dx(x0)])
    t = t / np.linalg.norm(t)
    n = np.array([-t[1], t[0]])
    angle = np.degrees(np.arctan2(t[1], t[0]))
    return p, t, n, angle


def draw_congruence(ax: plt.Axes) -> None:
    x = np.linspace(0.55, 9.42, 500)
    yc = center_y(x)
    w = beam_width(x)
    offsets = [-1.0, -0.62, -0.25, 0.25, 0.62, 1.0]
    for s in offsets:
        y = yc + s * w + 0.035 * s * np.sin(1.4 * x + 0.6)
        ax.plot(x, y, color=BLUE, lw=0.88, alpha=0.34, solid_capstyle="round", zorder=1)
    ax.plot(x, yc, color=NAVY, lw=2.05, solid_capstyle="round", zorder=3)

    # Observer vertex. The screen basis is local and is drawn at the tetrad,
    # not as a labelled source boundary.
    ax.add_patch(Circle((0.55, float(center_y(0.55))), 0.06, facecolor=INK, edgecolor="white", lw=0.8, zorder=6))
    label(ax, 0.78, -0.30, r"observer", size=7.8, ha="left", color=GRAY)

    # Past-directed affine-parameter direction.
    p0 = np.array([1.30, float(center_y(1.30)) - 0.34])
    p1 = np.array([2.18, float(center_y(2.18)) - 0.34])
    arrow(ax, tuple(p0), tuple(p1), color=GRAY, lw=0.9, ms=8.4, alpha=0.66, zorder=2)
    label(ax, p1[0] + 0.16, p1[1] - 0.02, r"$\lambda$", size=8.6, color=GRAY, ha="left")


def draw_tetrad(ax: plt.Axes) -> None:
    x0 = 5.1
    p, t, n, angle = tangent_frame(x0)
    screen_angle_deg = angle + 94.0
    screen_angle = np.deg2rad(screen_angle_deg)
    basis_origin = p
    basis_major = np.array([0.22, 0.98])
    basis_major = basis_major / np.linalg.norm(basis_major)
    basis_minor = np.array([0.58, -0.82])
    basis_minor = basis_minor / np.linalg.norm(basis_minor)

    # Local screen disk behind the tetrad arrows.
    ax.add_patch(
        Ellipse(
            p,
            1.46,
            0.62,
            angle=screen_angle_deg,
            facecolor=MINT,
            edgecolor=TEAL,
            lw=0.98,
            alpha=0.62,
            zorder=3,
        )
    )
    ax.add_patch(Circle(tuple(p), 0.035, facecolor=INK, edgecolor="white", lw=0.6, zorder=7))

    arrow(ax, tuple(p), tuple(p + 0.90 * t), color=GOLD, lw=1.85, ms=12.5, zorder=8)
    label(ax, *(p + 1.03 * t + np.array([0.0, 0.05])), r"$k^\mu$", size=9.8, color=GOLD)

    arrow(ax, tuple(p), tuple(p - 0.78 * t), color=GRAY, lw=1.18, ms=10.5, linestyle="--", alpha=0.86, zorder=8)
    label(ax, *(p - 0.92 * t + np.array([-0.03, 0.10])), r"$e^\mu$", size=9.5, color=GRAY)

    arrow(ax, tuple(basis_origin), tuple(basis_origin + 0.40 * basis_major), color=TEAL, lw=1.28, ms=8.2, alpha=0.94, zorder=8)
    arrow(ax, tuple(basis_origin), tuple(basis_origin + 0.23 * basis_minor), color=TEAL, lw=1.08, ms=7.0, alpha=0.88, zorder=8)
    label(ax, p[0] + 0.52, p[1] - 0.44, r"$Z^\mu,\bar Z^\mu$", size=8.4, color=TEAL)


def draw_theta_icon(ax: plt.Axes, cx: float, cy: float) -> None:
    ax.add_patch(Circle((cx, cy), 0.19, facecolor="none", edgecolor=GRAY, lw=0.88, alpha=0.58, zorder=5))
    ax.add_patch(Circle((cx, cy), 0.37, facecolor=MINT, edgecolor=TEAL, lw=1.25, alpha=0.82, zorder=4))
    for a in np.linspace(0, 2 * np.pi, 6, endpoint=False):
        start = (cx + 0.47 * np.cos(a), cy + 0.47 * np.sin(a))
        end = (cx + 0.66 * np.cos(a), cy + 0.66 * np.sin(a))
        arrow(ax, start, end, color=TEAL, lw=0.86, ms=7.0, zorder=6)
    label(ax, cx, cy - 0.72, r"$\theta$", size=10.2, color=TEAL, weight="bold")


def draw_sigma_plus_icon(ax: plt.Axes, cx: float, cy: float) -> None:
    ax.add_patch(Ellipse((cx, cy), 0.92, 0.38, facecolor="#f8e9e4", edgecolor=CORAL, lw=1.30, alpha=0.88, zorder=4))
    ax.add_patch(Ellipse((cx, cy), 0.40, 0.66, facecolor="none", edgecolor=LIGHT_GRAY, lw=0.84, alpha=0.82, zorder=3))
    arrow(ax, (cx - 0.30, cy), (cx - 0.58, cy), color=CORAL, lw=0.90, ms=7.0, zorder=6)
    arrow(ax, (cx + 0.30, cy), (cx + 0.58, cy), color=CORAL, lw=0.90, ms=7.0, zorder=6)
    arrow(ax, (cx, cy + 0.38), (cx, cy + 0.20), color=GRAY, lw=0.80, ms=6.5, zorder=6)
    arrow(ax, (cx, cy - 0.38), (cx, cy - 0.20), color=GRAY, lw=0.80, ms=6.5, zorder=6)
    label(ax, cx, cy - 0.72, r"$\sigma_{+}$", size=10.2, color=CORAL, weight="bold")


def draw_sigma_cross_icon(ax: plt.Axes, cx: float, cy: float) -> None:
    ax.add_patch(Ellipse((cx, cy), 0.92, 0.38, angle=45, facecolor="#f8e9e4", edgecolor=CORAL, lw=1.30, alpha=0.88, zorder=4))
    ax.add_patch(Ellipse((cx, cy), 0.40, 0.66, angle=45, facecolor="none", edgecolor=LIGHT_GRAY, lw=0.84, alpha=0.82, zorder=3))
    u = np.array([1.0, 1.0]) / np.sqrt(2)
    for sign in [-1, 1]:
        start = np.array([cx, cy]) + sign * 0.28 * u
        end = np.array([cx, cy]) + sign * 0.56 * u
        arrow(ax, tuple(start), tuple(end), color=CORAL, lw=0.90, ms=7.0, zorder=6)
    label(ax, cx, cy - 0.72, r"$\sigma_{\times}$", size=10.2, color=CORAL, weight="bold")


def draw_sachs_scalars(ax: plt.Axes) -> None:
    y0 = -2.78
    draw_theta_icon(ax, 2.65, y0)
    draw_sigma_plus_icon(ax, 5.05, y0)
    draw_sigma_cross_icon(ax, 7.45, y0)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 0.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.0, 3.55), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    draw_congruence(ax)
    draw_tetrad(ax)
    draw_sachs_scalars(ax)

    ax.set_xlim(0.25, 9.85)
    ax.set_ylim(-3.78, 1.88)
    ax.axis("off")

    FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.025)
    fig.savefig(OUT_PNG, bbox_inches="tight", pad_inches=0.025, dpi=300)
    plt.close(fig)
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
