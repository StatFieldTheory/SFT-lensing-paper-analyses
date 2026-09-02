"""Draft paper figure for the equal-shell driving-field cumulant ζ.

Reads ell_band_decomp_results.npz (direct-eval ζ on the squeezed-isoceles family
(1,cosγ,cosγ), 4 channels × γ × λ-shells). Produces ONE figure:

  ζ(γ, shell) for the 4 channels (TTT,TTP,TPP,PPP) as square imshow panels, shared
  y-axis with dual ticks (λ [Mpc] left, source redshift z right), one shared
  per-panel-normalised diverging-symlog colorbar. The source shells are the 16
  discrete λ-rows of the table (highly non-uniform Δλ → shown in shell-index
  space; the dual ticks give the physical λ and z).

The ℓ-band decomposition is presented as a LaTeX table (see emit_ell_band_table.py)
rather than a figure.

Run (PyCCL env has matplotlib):
  /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python plot_zeta_figures.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm

_HERE = Path(__file__).resolve().parent
# Converged cutoff (eight HIGH windows to ell = 15360, clipped at 85 arcmin), 2026-08-27;
# the June cut1000 ell_band_decomp_results.npz next to this script is superseded.
_NPZ = _HERE / "outputs" / "zeta_bands_cut15360_gmax85.npz"
_FIGDIR = _HERE / "outputs"   # the paper copy is deployed by reproduce/deploy.py
_FIGDIR.mkdir(parents=True, exist_ok=True)
# Panels show the spin-0 channels (TTT, TTP) and, for the spin-2 sector, the
# conjugate-helicity MODULUS channels Bmod=<Phi00 |Psi0|^2> and Dmod=<Psi0 |Psi0|^2>
# that actually source the FK coupling (the all-+2 zeta_TPP/zeta_PPP are
# gamma^4-suppressed; see Appendix D / derive_kappa3_spin2_helicity.wl).
_CHAN = ("TTT", "TTP", "Bmod", "Dmod")
_CHAN_TEX = {"TTT": r"$\zeta_{TTT}=\langle\Phi_{00}^3\rangle$",
             "TTP": r"$\zeta_{TTP}=\langle\Phi_{00}^2\,\Psi_+\rangle$",
             "Bmod": r"$\zeta_{B}=\langle\Phi_{00}\,|\Psi_0|^2\rangle$",
             "Dmod": r"$\zeta_{D}=\langle\Psi_0\,|\Psi_0|^2\rangle$"}
_GMAX = 350.0  # arcmin — the κκ-2PCF-relevant range (beyond, ζ is near-zero)

# The figure is rendered ~13.6 in wide but placed at \textwidth (~7 in) in a
# two-column figure*, i.e. downscaled ~0.5x; sizes below are chosen so the
# on-page text lands near the 9-10 pt body size after that downscale without
# the (dense) x- and colorbar-ticks overlapping.
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 17, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130, "savefig.dpi": 200, "legend.fontsize": 13,
})


def _load():
    d = np.load(_NPZ, allow_pickle=False)
    g = np.asarray(d["gamma_arcmin"]); z = np.asarray(d["z_shells"])
    lam = np.asarray(d["lambda_phys"])
    full = {c: np.asarray(d[f"full_{c}"]) for c in _CHAN}
    return g, z, lam, full


def fig_slices(g, z, lam, full):
    mg = g <= _GMAX
    gg = g[mg]               # uniform in log γ (geomspace) → imshow pixels OK
    ng, nshell = gg.size, lam.size

    fig, axes = plt.subplots(1, 4, figsize=(13.6, 4.9), sharey=True,
                             constrained_layout=True)
    # symmetric-log, normalised to [-1,1]; linthresh=1e-5 so the colour scale
    # resolves structure down to 10^-5 of each panel's peak (≈5 decades).
    norm = SymLogNorm(linthresh=1e-5, vmin=-1.0, vmax=1.0, base=10)

    # log-γ tick positions in index space.
    logg = np.log10(gg)
    def gidx(val):  # fractional index of γ=val on the uniform-log grid
        return float(np.interp(np.log10(val), logg, np.arange(ng)))
    gticks = [1, 3, 10, 30, 100, 300]
    gtick_pos = [gidx(v) for v in gticks if gg[0] <= v <= gg[-1]]
    gtick_lab = [f"{v:g}" for v in gticks if gg[0] <= v <= gg[-1]]

    # shell-index tick selections for the dual y-axis (shells cluster near the
    # source, so a sparse set keeps the λ labels legible).
    yi = [0, 4, 8, 12, 15]

    im = None
    for i, c in enumerate(_CHAN):
        ax = axes[i]
        arr = full[c][mg].T                      # (nshell, ng), row 0 = lowest z
        peak = float(np.nanmax(np.abs(arr)))
        im = ax.imshow(arr / peak, origin="lower", aspect="auto",
                       interpolation="bilinear", cmap="RdBu_r", norm=norm,
                       extent=[-0.5, ng - 0.5, -0.5, nshell - 0.5])
        ax.set_box_aspect(1.0)
        ax.set_title(_CHAN_TEX[c], pad=6)
        ax.set_xticks(gtick_pos); ax.set_xticklabels(gtick_lab)
        ax.set_xlabel(r"$\gamma$ [arcmin]")
        ax.text(0.04, 0.96, fr"peak $|\zeta|$={peak:.1e}", transform=ax.transAxes,
                fontsize=11, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", alpha=0.8))
        if i == 0:                               # left axis: λ [Mpc]
            ax.set_yticks(yi)
            ax.set_yticklabels([f"{lam[k]:.0f}" for k in yi])
            ax.set_ylabel(r"affine shell $\lambda$ [Mpc]")
        if i == 3:                               # right axis: source redshift z
            # secondary_yaxis shares the SAME axes box (unlike twinx), so all
            # four panels keep an identical square height under constrained_layout.
            axr = ax.secondary_yaxis("right")
            axr.set_yticks(yi)
            axr.set_yticklabels([f"{z[k]:.2f}" for k in yi])
            axr.set_ylabel(r"source redshift $z$", rotation=270, labelpad=16)

    cbticks = [-1, -1e-1, -1e-3, -1e-5, 0, 1e-5, 1e-3, 1e-1, 1]
    cb = fig.colorbar(im, ax=axes, location="bottom", orientation="horizontal",
                      fraction=0.10, pad=0.14, aspect=55, ticks=cbticks)
    cb.set_label(r"$\zeta / \max|\zeta|$  per panel "
                 r"(symmetric-log; blue $<0$, red $>0$; resolves to $10^{-5}$)",
                 fontsize=14)
    cb.ax.tick_params(labelsize=11)

    for ext in ("pdf", "png"):
        fig.savefig(_FIGDIR / f"zeta_driving_field_slices_draft.{ext}",
                    bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] wrote {_FIGDIR/'zeta_driving_field_slices_draft.pdf'}")


def main():
    g, z, lam, full = _load()
    fig_slices(g, z, lam, full)


if __name__ == "__main__":
    main()
