"""Multi-redshift convergence figure: 2 rows x 5 columns.

  row 1 : convergence 2PCF  xi_kappa(gamma)
  row 2 : angular power spectrum  C_ell^{kappa kappa}  (curved-sky Wigner-d
          transform of row 1, with the ell=0 monopole removed; see
          plot_analysis3_cl_decomposition)
  cols  : source redshift z_s = 1, 1.7, 2.5, 3.2, 4.0

Each panel shows ONLY Order-0 (blue markers) and Full = O0+FF+FK (grey line),
in the marker/colour convention of Fig. 11 / Fig. 12.  No FF/FK breakdown, no
fractional-residual sub-panels.

Input : outputs/multiz_kappa_2pcf_5z.npz  (from compute_multiz_kappa_2pcf.py)
Output: outputs/multiz_kappa_xi_cl_2x5.{png,pdf}

Run with the PyCCL interpreter (needs scipy + matplotlib).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _plot_style import (  # noqa: E402
    GAMMA_X_MAX_ARCMIN, PALETTE, apply_rcparams, clip_gamma_axis,
    plot_signed_line, plot_signed_markers,
)
import plot_analysis3_cl_decomposition as C  # noqa: E402

NPZ = HERE / "outputs" / "multiz_kappa_2pcf_5z.npz"
OUT_STEM = HERE / "outputs" / "multiz_kappa_xi_cl_2x5"
LOAD_KW: dict[str, Any] = {"allow_pickle": True}

COLOR_O0 = PALETTE[0]
MARKER_O0 = "o"
ELL = C.ELL


def main() -> int:
    d = np.load(NPZ, **LOAD_KW)
    z = np.asarray(d["z"], float)
    gamma = np.asarray(d["gamma"], float)              # arcmin
    o0, ff, fk = (np.asarray(d[k], float) for k in ("o0", "ff", "fk"))
    full = o0 + ff + fk
    nz = len(z)

    # one curved-sky operator (kappa-kappa -> Legendre d^l_{00}); gamma shared.
    setup = C.build_curved_matrix(gamma, ELL, 0, 0)
    pref = ELL * (ELL + 1.0) / (2.0 * np.pi)           # band power D_ell

    apply_rcparams()
    import matplotlib.pyplot as plt
    # large fonts: the figure spans \textwidth (a ~2x downscale in the PDF), so
    # source-point sizes need to be roughly doubled to read at print size.
    panel_rc = {"axes.labelsize": 24, "axes.titlesize": 23,
                "xtick.labelsize": 19, "ytick.labelsize": 19,
                "legend.fontsize": 22}
    saved = {k: plt.rcParams[k] for k in panel_rc}
    plt.rcParams.update(panel_rc)

    fig, axes = plt.subplots(2, nz, figsize=(3.1 * nz, 7.2),
                             sharey="row")
    YLIM_XI = (1e-9, 2e-3)      # row 1: |xi_kappa|
    YLIM_CL = (3e-8, 5e-4)      # row 2: band power D_ell
    legend_handles = None
    cl_o0_all, cl_full_all = [], []
    for j in range(nz):
        # --- row 1: 2PCF xi_kappa(gamma) ---
        ax = axes[0, j]
        plot_signed_line(ax, gamma, full[j], label="Full (O0+FF+FK)")
        plot_signed_markers(ax, gamma, o0[j], color=COLOR_O0, marker=MARKER_O0,
                            label="Order-0")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(*YLIM_XI)
        clip_gamma_axis(ax, gamma_min=gamma.min(), gamma_max=GAMMA_X_MAX_ARCMIN)
        ax.set_title(rf"$z_s = {z[j]:.1f}$")
        ax.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
        if j == 0:
            ax.set_ylabel(r"$\xi_\kappa(\gamma)$")

        # --- row 2: C_ell^{kk} (curved-sky transform, monopole removed) ---
        axc = axes[1, j]
        cl_o0 = C.forward_curved(o0[j], setup)
        cl_full = C.forward_curved(full[j], setup)
        cl_o0_all.append(cl_o0); cl_full_all.append(cl_full)
        plot_signed_line(axc, ELL, pref * cl_full, label="Full (O0+FF+FK)")
        plot_signed_markers(axc, ELL, pref * cl_o0, color=COLOR_O0,
                            marker=MARKER_O0, label="Order-0")
        axc.set_xscale("log"); axc.set_yscale("log")
        axc.set_xlim(ELL.min(), ELL.max())
        axc.set_ylim(*YLIM_CL)
        axc.set_xlabel(r"$\ell$")
        if j == 0:
            axc.set_ylabel(r"$\ell(\ell+1)\,C_\ell^{\kappa\kappa}/2\pi$")
        if legend_handles is None:
            h, l = ax.get_legend_handles_labels()
            seen, uh, ul = set(), [], []
            for hh, ll in zip(h, l):
                if ll not in seen:
                    seen.add(ll); uh.append(hh); ul.append(ll)
            legend_handles = (uh, ul)
            # filled/hollow = +/- convention is stated in the caption (no in-panel
            # annotation, which would overflow the narrow panels).

    if legend_handles is not None:
        hh, ll = legend_handles
        fig.legend(hh, ll, loc="upper center", ncol=len(hh),
                   bbox_to_anchor=(0.5, 0.995), frameon=False, fontsize=22)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(OUT_STEM.with_suffix(ext))
    plt.close(fig)
    plt.rcParams.update(saved)
    print(f"-> {OUT_STEM.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
