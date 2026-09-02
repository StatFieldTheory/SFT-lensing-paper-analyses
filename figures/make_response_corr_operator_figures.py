from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm
from scipy.ndimage import gaussian_filter, gaussian_filter1d


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py

SACHS_SCRIPT_DIR = ROOT / "SFT-lensing-paper-analyses" / "sachs_sft" / "scripts"
CORR_OP_DIR = ROOT / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables" / "C_propagator" / "corr_op"
for path in (SACHS_SCRIPT_DIR, CORR_OP_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from D_callable import D_at, _CACHE as D_CACHE, z_of_lambda  # noqa: E402
from corr_op_C_callable import C_fn_batch, TABLE_PATH  # noqa: E402


LAM_SOURCE_Z5 = 2318.0
GAMMA_RANGE_ARCMIN = (0.5, 300.0)

INK = "#16202a"
MUTED_BLUE = "#2f6f9f"
LIGHT_BLUE = "#dbeaf4"
TEAL = "#248f8d"
GOLD = "#c58a2a"
CORAL = "#b95f58"
VIOLET = "#7e6aa8"
GRAY = "#68717d"
PALE_GRAY = "#eef1f4"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 360,
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "dejavuserif",
            "font.size": 9.5,
            "axes.labelsize": 10.5,
            "axes.titlesize": 10.5,
            "axes.edgecolor": "#2b3138",
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": "#c8ced6",
            "grid.alpha": 0.28,
            "grid.linewidth": 0.55,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def lambda_of_z(z: np.ndarray | float) -> np.ndarray | float:
    z_arr = np.asarray(z, dtype=float)
    lam = np.interp(z_arr, D_CACHE["z"], D_CACHE["lambda"])
    return float(lam) if np.ndim(z) == 0 else lam


def z_label_axis(ax: plt.Axes, *, ticks_z: tuple[float, ...]) -> plt.Axes:
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    ticks_lam = [lambda_of_z(z) for z in ticks_z]
    top.set_xticks(ticks_lam)
    top.set_xticklabels([f"{z:g}" for z in ticks_z])
    top.set_xlabel("")
    top.text(
        0.0,
        1.035,
        r"$z(\lambda)$",
        transform=top.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
    )
    top.tick_params(direction="in", top=True, bottom=False, labelsize=8.5)
    return top


def panel_label(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.022,
        0.96,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=INK,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#d9dde3", lw=0.6, alpha=0.94),
        zorder=10,
    )


def signed_power_ticks(scale: float) -> tuple[list[float], list[str]]:
    exponent = int(np.floor(np.log10(max(scale, 1e-300))))
    tick = 10.0**exponent
    return [-tick, 0.0, tick], [rf"$-10^{{{exponent}}}$", r"$0$", rf"$10^{{{exponent}}}$"]


def smooth_curve(values: np.ndarray) -> np.ndarray:
    if float(np.nanmax(np.abs(values))) < 1e-20:
        return values
    return gaussian_filter1d(values, sigma=2.2, mode="nearest")


def smooth_image(values: np.ndarray) -> np.ndarray:
    if float(np.nanmax(np.abs(values))) < 1e-20:
        return values
    return gaussian_filter(values, sigma=1.05, mode="nearest")


def directions_for_gamma(gamma_arcmin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    theta = np.deg2rad(np.asarray(gamma_arcmin, dtype=float) / 60.0)
    n = theta.size
    n1 = np.repeat(np.array([[0.0, 0.0, 1.0]]), n, axis=0)
    n2 = np.column_stack((np.sin(theta), np.zeros(n), np.cos(theta)))
    return n1, n2


def corr_batch(
    gamma_arcmin: np.ndarray | float,
    lam1: np.ndarray | float,
    lam2: np.ndarray | float,
) -> np.ndarray:
    gamma_b, lam1_b, lam2_b = np.broadcast_arrays(
        np.asarray(gamma_arcmin, dtype=float),
        np.asarray(lam1, dtype=float),
        np.asarray(lam2, dtype=float),
    )
    n1, n2 = directions_for_gamma(gamma_b.ravel())
    out = C_fn_batch(n1, lam1_b.ravel(), n2, lam2_b.ravel())
    return out.reshape(gamma_b.shape + (3, 3))


def make_response_figure() -> None:
    lam_s = LAM_SOURCE_Z5
    z_s = float(z_of_lambda(lam_s))
    lam = np.linspace(0.1, lam_s, 1400)
    response = (D_at(lam) / D_at(lam_s)) ** 2
    i_peak = int(np.argmax(response))
    lam_peak = float(lam[i_peak])
    z_peak = float(z_of_lambda(lam_peak))

    fig, ax = plt.subplots(figsize=(5.55, 3.45), constrained_layout=True)
    ax.set_facecolor("#fbfcfd")
    ax.fill_between(lam, 0.0, response, color=LIGHT_BLUE, alpha=0.86, zorder=1)
    ax.plot(lam, response, color=MUTED_BLUE, lw=2.25, zorder=3)
    ax.plot(lam_peak, response[i_peak], marker="o", ms=4.8, color=GOLD, zorder=5)
    ax.axvline(lam_peak, color=GOLD, lw=0.9, ls=(0, (3, 2)), alpha=0.78, zorder=2)
    ax.axvline(lam_s, color=CORAL, lw=1.0, ls=(0, (4, 2)), alpha=0.82, zorder=2)
    ax.axhline(1.0, color=GRAY, lw=0.8, ls=":", alpha=0.8, zorder=2)

    ax.annotate(
        rf"peak $z\simeq {z_peak:.2f}$" "\n" rf"$\mathcal{{R}}={response[i_peak]:.2f}$",
        xy=(lam_peak, response[i_peak]),
        xytext=(lam_peak - 475, response[i_peak] - 0.33),
        arrowprops=dict(arrowstyle="-|>", color=GOLD, lw=0.9, shrinkA=2, shrinkB=4),
        ha="right",
        va="top",
        fontsize=8.7,
        color="#5b4217",
        bbox=dict(boxstyle="round,pad=0.24", fc="#fff8e7", ec="#e1c377", lw=0.55),
    )
    ax.text(
        lam_s - 15,
        0.11,
        rf"$\lambda_s={lam_s:.0f}\,\mathrm{{Mpc}}$" "\n" rf"$z_s\simeq {z_s:.2f}$",
        ha="right",
        va="bottom",
        fontsize=8.7,
        color="#723b36",
        bbox=dict(boxstyle="round,pad=0.24", fc="#fff3f1", ec="#e1b8b3", lw=0.55),
    )
    ax.set_xlim(0.0, lam_s)
    ax.set_ylim(0.0, max(2.32, float(response.max()) * 1.09))
    ax.set_xlabel(r"$\lambda$ [Mpc]")
    ax.set_ylabel(r"$\mathcal{R}(\lambda,\lambda_s)$")
    z_label_axis(ax, ticks_z=(0.5, 1, 2, 5))
    ax.tick_params(labelsize=8.9)

    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"response_operator_draft.{suffix}", bbox_inches="tight")
    plt.close(fig)


def make_corr_operator_figure() -> None:
    table = np.load(TABLE_PATH)
    lam_grid = np.asarray(table["lambda_grid"], dtype=float)
    lam_grid = lam_grid[lam_grid <= LAM_SOURCE_Z5]
    lam_s = float(lam_grid[-1])

    gamma = np.geomspace(*GAMMA_RANGE_ARCMIN, 170)
    shell_lams = np.array([1300.0, 1900.0, lam_s])
    shell_colors = [TEAL, GOLD, CORAL]
    shell_markers = ["s", "^", "D"]
    shell_linestyles = ["solid", (0, (1.2, 1.4)), (0, (5.0, 1.7, 1.1, 1.7))]

    equal_shell = {}
    for shell in shell_lams:
        equal_shell[shell] = gaussian_filter1d(
            corr_batch(gamma, shell, shell),
            sigma=2.2,
            axis=0,
            mode="nearest",
        )

    matrix_gamma = 5.0
    radial_axis = np.linspace(lam_grid[0], lam_s, 82)
    l1_mat, l2_mat = np.meshgrid(radial_axis, radial_axis, indexing="ij")
    radial_tensor = corr_batch(matrix_gamma, l1_mat, l2_mat)

    gamma_slice = np.geomspace(*GAMMA_RANGE_ARCMIN, 118)
    lambda_slice = np.linspace(lam_grid[0], lam_s, 86)
    gam_mesh, lam_mesh = np.meshgrid(gamma_slice, lambda_slice, indexing="xy")
    angular_tensor = corr_batch(gam_mesh, lam_mesh, lam_s)

    channels = [
        ((0, 0), r"$\Phi_{00}\Phi_{00}$"),
        ((0, 1), r"$\Phi_{00}\,\Psi_+$"),
        ((1, 1), r"$\Psi_+\,\Psi_+$"),
        ((2, 2), r"$\Psi_\times\,\Psi_\times$"),
    ]

    fig = plt.figure(figsize=(9.35, 8.55))
    bottom, top = 0.075, 0.955
    row_gap = 0.034
    row_h = (top - bottom - row_gap * (len(channels) - 1)) / len(channels)
    curve_x, curve_w = 0.065, 0.350
    heat1_x = curve_x + curve_w + 0.067
    heat_size = row_h
    heat_gap = 0.010
    heat2_x = heat1_x + heat_size + heat_gap
    cbar_x = heat2_x + heat_size + 0.014
    axs = np.empty((len(channels), 3), dtype=object)
    for row in range(len(channels)):
        y = top - (row + 1) * row_h - row * row_gap
        axs[row, 0] = fig.add_axes(
            [curve_x, y, curve_w, row_h],
            sharex=axs[0, 0] if row else None,
        )
        axs[row, 1] = fig.add_axes(
            [heat1_x, y, heat_size, row_h],
            sharex=axs[0, 1] if row else None,
        )
        axs[row, 2] = fig.add_axes(
            [heat2_x, y, heat_size, row_h],
            sharex=axs[0, 2] if row else None,
            sharey=axs[row, 1],
        )
    for ax in axs.ravel():
        ax.set_facecolor("#fbfcfd")
    for ax in axs[:, 1:].ravel():
        ax.grid(False, which="both")

    column_titles = [
        r"$\mathcal{C}_{XY}(\gamma;\lambda,\lambda)$",
        rf"$\mathcal{{C}}_{{XY}}(\lambda_1,\lambda_2;\gamma={matrix_gamma:g}')$",
        rf"$\mathcal{{C}}_{{XY}}(\gamma,\lambda_1;\lambda_2={lam_s:.0f}\,\mathrm{{Mpc}})$",
    ]
    for col, title in enumerate(column_titles):
        axs[0, col].set_title(title, loc="left", pad=8)

    def signed_norm(values: np.ndarray, *, linthresh_frac: float = 1e-3) -> SymLogNorm:
        absmax = float(np.nanpercentile(np.abs(values), 99.7))
        absmax = max(absmax, 1e-20)
        linthresh = max(absmax * linthresh_frac, 1e-18)
        return SymLogNorm(linthresh=linthresh, linscale=0.55, vmin=-absmax, vmax=absmax, base=10)

    smoothed_radial: dict[tuple[int, int], np.ndarray] = {}
    smoothed_angular: dict[tuple[int, int], np.ndarray] = {}
    global_chunks: list[np.ndarray] = []
    for (i, j), _label in channels:
        key = (i, j)
        smoothed_radial[key] = smooth_image(radial_tensor[..., i, j])
        smoothed_angular[key] = smooth_image(angular_tensor[..., i, j])
        global_chunks.extend([smoothed_radial[key].ravel(), smoothed_angular[key].ravel()])
    global_norm = signed_norm(np.concatenate(global_chunks), linthresh_frac=1e-3)
    last_im = None

    for row, ((i, j), label) in enumerate(channels):
        ax_cut, ax_lam, ax_gam_lam = axs[row]
        equal_values = [equal_shell[shell][..., i, j] for shell in shell_lams]
        radial_values = smoothed_radial[(i, j)]
        angular_values = smoothed_angular[(i, j)]

        curve_absmax = max(float(np.nanmax(np.abs(np.concatenate(equal_values)))), 1e-20)
        curve_linthresh = max(curve_absmax * 1e-3, 1e-22)

        for style_i, (shell, color, values) in enumerate(zip(shell_lams, shell_colors, equal_values)):
            z_shell = float(z_of_lambda(shell))
            ax_cut.plot(
                gamma,
                values,
                color=color,
                lw=1.25,
                ls=shell_linestyles[style_i],
                marker=shell_markers[style_i],
                markevery=(style_i * 6, 34),
                ms=2.5,
                mec="white",
                mew=0.35,
                label=rf"$\lambda={shell:.0f}$ Mpc, $z\simeq{z_shell:.2g}$",
            )
        ax_cut.axhline(0.0, color="#39424e", lw=0.65, alpha=0.8)
        ax_cut.grid(False)
        ax_cut.set_xscale("log")
        ax_cut.set_yscale("symlog", linthresh=curve_linthresh, linscale=0.55)
        ax_cut.set_xlim(*GAMMA_RANGE_ARCMIN)
        ax_cut.set_ylim(-1.25 * curve_absmax, 1.25 * curve_absmax)
        ticks, tick_labels = signed_power_ticks(curve_absmax)
        ax_cut.set_yticks(ticks)
        ax_cut.set_yticklabels(tick_labels)
        ax_cut.tick_params(axis="y", labelsize=8.7)
        ax_cut.text(
            0.035,
            0.90,
            label,
            transform=ax_cut.transAxes,
            ha="left",
            va="top",
            fontsize=8.2,
            color=INK,
            bbox=dict(boxstyle="round,pad=0.17", fc="white", ec="#d9dde3", lw=0.55, alpha=0.92),
        )
        if row == 0:
            ax_cut.legend(loc="lower left", fontsize=8.3, handlelength=1.65, frameon=False)

        im_lam = ax_lam.imshow(
            radial_values,
            origin="lower",
            extent=(radial_axis[0], radial_axis[-1], radial_axis[0], radial_axis[-1]),
            cmap="RdBu_r",
            norm=global_norm,
            interpolation="bicubic",
            aspect="equal",
            rasterized=True,
        )
        ax_lam.plot(radial_axis, radial_axis, color="0.18", lw=0.65, alpha=0.52)
        ax_lam.grid(False)
        ax_lam.set_xlim(radial_axis[0], radial_axis[-1])
        ax_lam.set_ylim(radial_axis[0], radial_axis[-1])
        ax_lam.set_box_aspect(1)

        log_gamma_slice = np.log10(gamma_slice)
        ax_gam_lam.imshow(
            angular_values,
            origin="lower",
            extent=(
                log_gamma_slice[0],
                log_gamma_slice[-1],
                lambda_slice[0],
                lambda_slice[-1],
            ),
            cmap="RdBu_r",
            norm=global_norm,
            interpolation="bicubic",
            aspect="auto",
            rasterized=True,
        )
        last_im = im_lam
        ax_gam_lam.grid(False)
        ax_gam_lam.set_xlim(log_gamma_slice[0], log_gamma_slice[-1])
        ax_gam_lam.set_ylim(lambda_slice[0], lambda_slice[-1])
        ax_gam_lam.set_box_aspect(1)

        if max(float(np.nanmax(np.abs(radial_values))), float(np.nanmax(np.abs(angular_values)))) < 1e-20:
            for ax in (ax_lam, ax_gam_lam):
                ax.text(
                    0.5,
                    0.5,
                    r"zero in this basis",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=7.2,
                    color=GRAY,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#d9dde3", lw=0.5, alpha=0.85),
                )

        if row < len(channels) - 1:
            ax_cut.tick_params(labelbottom=False)
            ax_lam.tick_params(labelbottom=False)
            ax_gam_lam.tick_params(labelbottom=False)
        ax_gam_lam.tick_params(axis="y", left=False, right=True, labelleft=False, labelright=False)

    for row in range(len(channels)):
        axs[row, 0].set_ylabel(r"$\mathcal{C}_{XY}$")
    axs[-1, 0].set_xlabel(r"$\gamma$ [arcmin]")
    axs[-1, 1].set_xlabel(r"$\lambda_2$ [Mpc]")
    axs[-1, 2].set_xlabel(r"$\gamma$ [arcmin]")
    gamma_ticks = np.array([1.0, 10.0, 100.0])
    axs[-1, 2].set_xticks(np.log10(gamma_ticks))
    axs[-1, 2].set_xticklabels([r"$10^0$", r"$10^1$", r"$10^2$"])
    fig.text(heat1_x - 0.040, 0.515, r"$\lambda_1$ [Mpc]", rotation=90, ha="center", va="center")
    cax = fig.add_axes([cbar_x, bottom + 0.055, 0.018, top - bottom - 0.110])
    cbar = fig.colorbar(last_im, cax=cax)
    hi_exp = int(np.floor(np.log10(max(float(global_norm.vmax), 1e-300))))
    mid_exp = hi_exp - 2
    c_ticks = [-10.0**hi_exp, -10.0**mid_exp, 0.0, 10.0**mid_exp, 10.0**hi_exp]
    c_tick_labels = [
        rf"$-10^{{{hi_exp}}}$",
        rf"$-10^{{{mid_exp}}}$",
        r"$0$",
        rf"$10^{{{mid_exp}}}$",
        rf"$10^{{{hi_exp}}}$",
    ]
    cbar.set_ticks(c_ticks)
    cbar.set_ticklabels(c_tick_labels)
    cbar.ax.tick_params(labelsize=7.2, length=2.5)

    for suffix in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"corr_operator_slices_draft.{suffix}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    setup_style()
    FIG_DIR.mkdir(exist_ok=True)
    make_response_figure()
    make_corr_operator_figure()
    print("wrote figures/response_operator_draft.pdf")
    print("wrote figures/response_operator_draft.png")
    print("wrote figures/corr_operator_slices_draft.pdf")
    print("wrote figures/corr_operator_slices_draft.png")


if __name__ == "__main__":
    main()
