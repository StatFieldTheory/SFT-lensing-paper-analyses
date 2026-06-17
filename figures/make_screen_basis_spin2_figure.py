from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "figures"

INK = "#16202a"
BLUE = "#2f6f9f"
GOLD = "#c58a2a"
CORAL = "#b95f58"
TEAL = "#2a8580"
GRAY = "#68717d"
LIGHT_GRAY = "#dfe5ec"
SIGN_GRAY = "#4f5964"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 360,
            "font.family": "serif",
            "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 11.3,
            "axes.linewidth": 0.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = INK,
    lw: float = 1.0,
    mutation_scale: float = 8.0,
    alpha: float = 1.0,
    zorder: int = 4,
    style: str = "-|>",
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            color=color,
            linewidth=lw,
            mutation_scale=mutation_scale,
            alpha=alpha,
            shrinkA=0.0,
            shrinkB=0.0,
            zorder=zorder,
        )
    )


def double_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = GRAY,
    lw: float = 0.9,
    mutation_scale: float = 7.5,
    alpha: float = 1.0,
    zorder: int = 4,
) -> None:
    arrow(
        ax,
        start,
        end,
        color=color,
        lw=lw,
        mutation_scale=mutation_scale,
        alpha=alpha,
        zorder=zorder,
        style="<->",
    )


def text_with_halo(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    color: str = INK,
    fontsize: float = 11.0,
    ha: str = "center",
    va: str = "center",
    weight: str | None = None,
    zorder: int = 10,
) -> None:
    t = ax.text(
        x,
        y,
        text,
        color=color,
        fontsize=fontsize,
        ha=ha,
        va=va,
        fontweight=weight,
        zorder=zorder,
    )
    t.set_path_effects([pe.withStroke(linewidth=2.6, foreground="white", alpha=0.92)])


def panel_label(ax: plt.Axes, label: str, x: float, y: float, *, transform=None) -> None:
    ax.text(
        x,
        y,
        label,
        transform=transform,
        ha="left",
        va="top",
        fontsize=12.6,
        fontweight="bold",
        color=INK,
    )


def parity_label(ax: plt.Axes, x: float, y: float, parity: str) -> None:
    ax.text(
        x,
        y,
        rf"$\mathrm{{parity}}\ \mathrm{{{parity}}}$",
        ha="center",
        va="center",
        fontsize=9.4,
        color=GRAY,
    )


def draw_pair_basis(ax: plt.Axes, cx: float, cy: float) -> None:
    left, right = cx - 0.93, cx + 0.93
    ax.plot([left - 0.36, right + 0.36], [cy, cy], color="#c8d0d9", lw=0.72, ls=(0, (2.0, 2.0)))
    ax.scatter([left, right], [cy, cy], s=18, color=BLUE, edgecolor="white", lw=0.55, zorder=5)
    double_arrow(ax, (left, cy - 0.34), (right, cy - 0.34), color=BLUE, lw=0.74, mutation_scale=7.0)
    ax.text(cx, cy - 0.52, r"$\gamma$", ha="center", va="top", fontsize=9.5, color=BLUE)

    origin = (cx, cy)
    arrow(ax, origin, (cx + 0.88, cy), lw=0.86, mutation_scale=7.0)
    arrow(ax, origin, (cx, cy + 0.72), lw=0.86, mutation_scale=7.0)
    ax.text(cx + 0.95, cy - 0.02, r"$\hat x\parallel\gamma$", ha="left", va="center", fontsize=9.5, color=INK)
    ax.text(cx + 0.06, cy + 0.78, r"$\hat y$", ha="left", va="bottom", fontsize=9.5, color=INK)
    ax.text(cx + 2.40, cy + 0.23, r"$\mathsf{P}_{\gamma}:\ \hat y\mapsto-\hat y$", ha="center", va="center", fontsize=9.4, color=GRAY)
    ax.plot([cx + 1.55, cx + 3.25], [cy, cy], color="#d7dde4", lw=0.72, ls=(0, (2.0, 2.0)))


def draw_trace_response(
    ax: plt.Axes,
    cx: float,
    cy: float,
    *,
    inward: bool = False,
    show_circle: bool = True,
) -> None:
    if show_circle:
        ax.add_patch(Circle((cx, cy), 0.36, facecolor="none", edgecolor=LIGHT_GRAY, lw=0.85, zorder=1))
        ax.add_patch(Circle((cx, cy), 0.48, facecolor="#fff6dd", edgecolor=GOLD, lw=1.15, alpha=0.92, zorder=2))
    inner_radius = 0.50 if show_circle else 0.16
    outer_radius = 0.72 if show_circle else 0.56
    for angle in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        inner = (cx + inner_radius * np.cos(angle), cy + inner_radius * np.sin(angle))
        outer = (cx + outer_radius * np.cos(angle), cy + outer_radius * np.sin(angle))
        start, end = (outer, inner) if inward else (inner, outer)
        arrow(ax, start, end, color=GOLD, lw=0.86, mutation_scale=6.4, alpha=0.88, zorder=4)


def draw_shear_response(
    ax: plt.Axes,
    cx: float,
    cy: float,
    angle: float,
    color: str,
    value_label: str | None = None,
) -> None:
    u = np.array([np.cos(angle), np.sin(angle)])
    v = np.array([-np.sin(angle), np.cos(angle)])
    center = np.array([cx, cy])
    for sign in (-1.0, 1.0):
        start = center + sign * 0.13 * u
        end = center + sign * 0.46 * u
        arrow(ax, tuple(start), tuple(end), color=color, lw=1.15, mutation_scale=7.6, alpha=0.95, zorder=4)
    # Contraction (minor) axis: inward gray arrows. Kept bold enough to stay
    # legible at print size for the cross-type glyphs, where the gray arrows lie
    # on the diagonal crossing the colored stretch arrows.
    for sign in (-1.0, 1.0):
        start = center + sign * 0.44 * v
        end = center + sign * 0.13 * v
        arrow(ax, tuple(start), tuple(end), color=GRAY, lw=1.02, mutation_scale=7.2, alpha=0.9, zorder=4)
    if value_label is not None:
        ax.text(cx, cy - 0.56, value_label, ha="center", va="top", fontsize=9.3, color=color)


def draw_component_column(
    ax: plt.Axes,
    x: float,
    *,
    color: str,
    name: str,
    formula: str,
    target: str,
    mode: str,
) -> None:
    ax.text(x, 2.28, name, ha="center", va="center", fontsize=11.9, color=INK)
    ax.text(x, 2.00, formula, ha="center", va="center", fontsize=9.1, color=INK)

    if mode == "trace":
        draw_trace_response(ax, x, 1.23)
    elif mode == "plus":
        draw_shear_response(ax, x - 0.40, 1.23, 0.0, color, r"$>0$")
        draw_shear_response(ax, x + 0.40, 1.23, np.pi / 2, color, r"$<0$")
    elif mode == "cross":
        draw_shear_response(ax, x - 0.40, 1.23, np.pi / 4, color, r"$>0$")
        draw_shear_response(ax, x + 0.40, 1.23, -np.pi / 4, color, r"$<0$")

    ax.text(x, 0.28, target, ha="center", va="center", fontsize=10.5, color=color)


def draw_panel_a(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 3.35)
    ax.set_ylim(-0.25, 5.98)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.text(0.50, 5.17, r"$\Phi_{00}<0$", ha="center", va="center", fontsize=11.6, color=GOLD)
    parity_label(ax, 0.50, 4.84, "even")
    draw_trace_response(ax, 1.66, 5.10, inward=True, show_circle=False)
    ax.text(2.72, 5.17, r"$\dot\theta\downarrow$", ha="center", va="center", fontsize=11.6, color=GOLD)

    plus_x, minus_x = 1.48, 2.42
    text_with_halo(ax, 1.52, 4.05, r"$+$", fontsize=11.9, color=SIGN_GRAY)
    text_with_halo(ax, 2.36, 4.05, r"$-$", fontsize=11.9, color=SIGN_GRAY)
    draw_source_row(ax, 3.68, r"$\Psi_+$", BLUE, 0.0, np.pi / 2)
    parity_label(ax, 0.62, 3.38, "even")
    ax.text(0.55, 2.88, r"$\dot\sigma_+$", ha="center", va="center", fontsize=11.5, color=BLUE)
    draw_shear_response(ax, plus_x, 2.84, 0.0, BLUE)
    draw_shear_response(ax, minus_x, 2.84, np.pi / 2, BLUE)

    text_with_halo(ax, 1.52, 1.91, r"$+$", fontsize=11.9, color=SIGN_GRAY)
    text_with_halo(ax, 2.36, 1.91, r"$-$", fontsize=11.9, color=SIGN_GRAY)
    draw_source_row(ax, 1.54, r"$\Psi_\times$", CORAL, np.pi / 4, -np.pi / 4)
    parity_label(ax, 0.62, 1.24, "odd")
    ax.text(0.55, 0.54, r"$\dot\sigma_\times$", ha="center", va="center", fontsize=11.5, color=CORAL)
    draw_shear_response(ax, plus_x, 0.48, np.pi / 4, CORAL)
    draw_shear_response(ax, minus_x, 0.48, -np.pi / 4, CORAL)


def draw_shared_basis(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    origin = (0.16, 0.28)
    arrow(ax, origin, (0.54, 0.28), lw=0.78, mutation_scale=6.4)
    arrow(ax, origin, (0.16, 0.72), lw=0.78, mutation_scale=6.4)
    ax.text(0.58, 0.26, r"$\hat x$", ha="left", va="center", fontsize=9.9, color=INK)
    ax.text(0.18, 0.78, r"$\hat y$", ha="left", va="bottom", fontsize=9.9, color=INK)


def draw_source_row(
    ax: plt.Axes,
    y: float,
    label: str,
    color: str,
    positive_angle: float,
    negative_angle: float,
) -> None:
    ax.text(0.62, y, label, ha="center", va="center", fontsize=11.6, color=color)
    draw_line_key_pair(ax, 1.94, y, positive_angle, negative_angle, color)


def draw_line_key_pair(
    ax: plt.Axes,
    cx: float,
    cy: float,
    positive_angle: float,
    negative_angle: float,
    color: str,
) -> None:
    for x, angle in ((cx - 0.42, positive_angle), (cx + 0.42, negative_angle)):
        length = 0.48
        dx = 0.5 * length * np.cos(angle)
        dy = 0.5 * length * np.sin(angle)
        ax.plot([x - dx, x + dx], [cy - dy, cy + dy], color=color, lw=2.0, alpha=0.96, solid_capstyle="round")


def phi_colormap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        "phi00_soft",
        ["#ffffff", "#fff7cc", "#f1c64b", "#d8902f", "#b95d2a"],
    )


def draw_component_segment(
    ax: plt.Axes,
    x: float,
    y: float,
    value: float,
    *,
    positive_angle: float,
    negative_angle: float,
    color: str,
    scale: float,
    max_len: float,
    cutoff: float = 0.06,
) -> None:
    # Length is strictly proportional to the local component amplitude.
    # ``scale`` is the strongest drawn component value, so ``strength`` already
    # lies in [0, 1]: no length floor and no saturation clamp compress the
    # dynamic range, and the strongest stick draws exactly ``max_len``.
    strength = abs(value) / scale
    if strength < cutoff:
        return
    angle = positive_angle if value >= 0.0 else negative_angle
    length = max_len * strength
    dx = 0.5 * length * np.cos(angle)
    dy = 0.5 * length * np.sin(angle)
    ax.plot(
        [x - dx, x + dx],
        [y - dy, y + dy],
        color=color,
        lw=0.5 + 1.3 * strength,
        alpha=0.32 + 0.56 * np.sqrt(strength),
        solid_capstyle="round",
        zorder=4,
    )


def draw_field_pattern(ax: plt.Axes) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-3.40, 3.40)
    ax.set_ylim(-2.55, 2.55)
    ax.axis("off")
    x = np.linspace(-3.55, 3.55, 420)
    y = np.linspace(-2.55, 2.55, 300)
    xx, yy = np.meshgrid(x, y)
    rr2 = xx * xx + yy * yy
    sigma = 0.92
    phi00 = np.exp(-0.5 * rr2 / sigma**2)

    levels = np.linspace(0.08, 1.001, 16)
    ax.contourf(xx, yy, phi00, levels=levels, cmap=phi_colormap(), alpha=0.88, antialiased=True, zorder=1)
    ax.contour(xx, yy, phi00, levels=[0.18, 0.38, 0.62, 0.82], colors="#d4d8dc", linewidths=0.72, zorder=2)
    ax.add_patch(Circle((0.0, 0.0), 0.095, facecolor=GOLD, edgecolor="white", lw=0.7, zorder=6))
    text_with_halo(ax, -0.28, 0.33, r"$\Phi_{00}<0$", color=GOLD, fontsize=11.0, ha="right")

    # Weyl shear Psi_0 sticks.  The radial amplitude is strongest near the mass
    # and decays outward (regularized point-mass tidal shear, |gamma| ~ 1/r^2),
    # so the sticks intensify toward the overdensity rather than vanishing at
    # its centre.  The cos2theta / sin2theta angular factors carry the spin-2
    # orientation pattern (the pedagogical content of the panel).
    r_core = 0.9

    def shear_amp(r2: float) -> float:
        return 1.0 / (1.0 + r2 / r_core**2)

    xs = np.linspace(-3.10, 3.10, 17)
    ys = np.linspace(-2.05, 2.05, 11)
    # Inner cutoff small enough to keep the four axial cells (r~0.39-0.41), where
    # Re Psi_0 (cos2theta) is maximal, so the innermost ring shows BOTH Re (blue,
    # on the axes) and Im (red, on the diagonals at r~0.56), not Im alone.
    r2_min, r2_max = 0.14, 18.6

    # First pass: fix one shared length scale from the strongest drawn component
    # value (over both Re and Im) so that stick length is strictly proportional
    # to the local amplitude and the longest stick draws exactly max_len.
    cells: list[tuple[float, float, float, float]] = []
    scale = 1e-12
    for yy0 in ys:
        for xx0 in xs:
            r2 = xx0 * xx0 + yy0 * yy0
            if r2 < r2_min or r2 > r2_max:
                continue
            local_amp = shear_amp(r2)
            local_theta = np.arctan2(yy0, xx0) + 0.5 * np.pi
            re_val = local_amp * np.cos(2.0 * local_theta)
            im_val = local_amp * np.sin(2.0 * local_theta)
            cells.append((xx0, yy0, re_val, im_val))
            scale = max(scale, abs(re_val), abs(im_val))

    # Cap the longest stick below the grid pitch so nothing overlaps, even in
    # the dense, high-amplitude region now sitting near the centre.
    max_len = 0.85 * min(xs[1] - xs[0], ys[1] - ys[0])

    for xx0, yy0, re_val, im_val in cells:
        draw_component_segment(
            ax,
            xx0,
            yy0,
            re_val,
            positive_angle=0.0,
            negative_angle=np.pi / 2,
            color=BLUE,
            scale=scale,
            max_len=max_len,
        )
        draw_component_segment(
            ax,
            xx0,
            yy0,
            im_val,
            positive_angle=np.pi / 4,
            negative_angle=-np.pi / 4,
            color=CORAL,
            scale=scale,
            max_len=max_len,
        )

def render() -> None:
    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.35, 4.02), facecolor="white")
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=(0.82, 1.78),
        left=0.030,
        right=0.992,
        bottom=0.040,
        top=0.975,
        wspace=0.040,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_a(ax_a)
    draw_field_pattern(ax_b)
    basis_ax = fig.add_axes([0.304, 0.790, 0.092, 0.132])
    draw_shared_basis(basis_ax)
    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"screen_basis_spin2_pattern.{suffix}", bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


if __name__ == "__main__":
    render()
