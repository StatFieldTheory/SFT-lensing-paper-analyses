"""Polarization decomposition of the two Order-2 corrections at z_s = 5: the
FF (nonlinear-propagation) and FK (three-point-cumulant) contributions to the
shear angular power spectra C_ell^{EE}, C_ell^{BB}, C_ell^{EB} (1x3 panels).

  EE+BB  <- xi_+ = <g+ g+> + <gx gx>   (pairs [1,1]+[2,2], kernel d^l_{2, 2})
  EE-BB  <- xi_- = <g+ g+> - <gx gx>   (pairs [1,1]-[2,2], kernel d^l_{2,-2})
  EE = (EE+BB + EE-BB)/2,  BB = (EE+BB - EE-BB)/2
  EB    <- <g+ gx>                      (pair [1,2],        kernel d^l_{2,-2})

The cross <g+ gx> carries C^EB via Im<gg> = 2 C^EB (same spin-4 kernel as xi_-).
Both EB channels vanish by parity (an odd number of Psi_x legs reflects to minus
itself). The FK cross [1,2] is identically zero in the input construction.
A nonzero transformed FF cross is a numerical parity residual whose amplitude
is measured from the selected product. In the scalar example FK has one linear pure-E external
leg, so its physical BB is zero. Its transformed nonzero BB is displayed
as a numerical residual, not a signal. FF can generate BB at this order.

Reuses the curved-sky Wigner-d transform of plot_analysis3_cl_decomposition.
Run with the PyCCL interpreter.  Output: outputs/cl_EB_polarization.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _plot_style import PALETTE, apply_rcparams, plot_signed_markers  # noqa: E402
import plot_analysis3_cl_decomposition as C  # noqa: E402

OUT_STEM = HERE / "outputs" / "cl_EB_polarization"
ELL = C.ELL
PREF = ELL * (ELL + 1.0) / (2.0 * np.pi)
COL_FF, COL_FK = PALETTE[1], PALETTE[2]
MK_FF, MK_FK = "^", "s"


def polarization(grouped, s22, s2m2):
    """Return (C_EE, C_BB, C_EB, raw_cross_max) for one Order-2 term."""
    ee_pbb = C.forward_curved(C._combine(grouped, [((1, 1), +1.0), ((2, 2), +1.0)]), s22, dc_subtract=False)
    ee_mbb = C.forward_curved(C._combine(grouped, [((1, 1), +1.0), ((2, 2), -1.0)]), s2m2, dc_subtract=False)
    cross = C._combine(grouped, [((1, 2), +1.0)])
    eb = C.forward_curved(cross, s2m2, dc_subtract=False)
    return 0.5 * (ee_pbb + ee_mbb), 0.5 * (ee_pbb - ee_mbb), eb, float(np.max(np.abs(cross)))


def main() -> int:
    g_fk, fkg = C.load_sweep_order(C.FK_NPZ, 2)
    g_ff, ffg = C.load_sweep_order(C.FF_NPZ, 2)
    assert np.allclose(g_fk, g_ff), "gamma grids differ"
    s22 = C.build_curved_matrix(g_fk, ELL, 2, 2)
    s2m2 = C.build_curved_matrix(g_fk, ELL, 2, -2)

    ee_ff, bb_ff, eb_ff, cross_ff = polarization(ffg, s22, s2m2)
    ee_fk, bb_fk, eb_fk, cross_fk = polarization(fkg, s22, s2m2)

    print("=== Order-2 polarization C_ell (z_s=5) ===")
    for nm, ee, bb, eb, cr in (("FF", ee_ff, bb_ff, eb_ff, cross_ff),
                               ("FK", ee_fk, bb_fk, eb_fk, cross_fk)):
        print(f"  {nm}: max|D_EE|={np.max(np.abs(PREF*ee)):.2e}  "
              f"max|D_BB|={np.max(np.abs(PREF*bb)):.2e}  "
              f"max|D_EB|={np.max(np.abs(PREF*eb)):.2e}  "
              f"raw|<g+gx>|max={cr:.2e}")

    apply_rcparams()
    import matplotlib.pyplot as plt
    rc = {"axes.labelsize": 21, "axes.titlesize": 21,
          "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 18}
    saved = {k: plt.rcParams[k] for k in rc}
    plt.rcParams.update(rc)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.4), sharey=True)
    YLIM = (1e-16, 1e-5)
    handles = None

    # --- EE and BB panels: FF and FK signed markers ---
    for ax, ff, fk, title in ((axes[0], ee_ff, ee_fk, r"$\Delta C_\ell^{EE}$"),
                              (axes[1], bb_ff, bb_fk, r"$\Delta C_\ell^{BB}$")):
        plot_signed_markers(ax, ELL, PREF * ff, color=COL_FF, marker=MK_FF, label="FF")
        # Below ell ~ 50 the FK values come from the large-separation part of
        # its 2PCF, where the multipole sum no longer converges; drawn faint
        # and not quoted.
        is_bb = ax is axes[1]
        plot_signed_markers(ax, ELL, PREF * fk,
                            color="0.55" if is_bb else COL_FK, marker=MK_FK,
                            label="FK numerical residual" if is_bb else "FK",
                            faint_outside=(50.0, None))
        ax.set_title(title)
        if is_bb:
            ax.legend(loc="best", fontsize=12)
        if handles is None:
            handles = ax.get_legend_handles_labels()

    # --- EB panel: parity null. ------------------------------------------
    # FK is identically zero here as a THEOREM, not as a measurement: under a
    # reflection about the separation axis Re Psi_0 is even and Im Psi_0 is
    # odd, so <Phi_00 Psi_+ Psi_x> equals minus itself. The vertex callable
    # encodes that structurally, leaving every coupling-tensor entry with an
    # odd number of Im slots at zero. The panel therefore demonstrates that
    # the transform and the fold do not leak E or B power into EB, which is a
    # useful check but a different statement; the label says so.
    axeb = axes[2]
    axeb.loglog(ELL, PREF * np.abs(eb_ff), color=COL_FF, lw=1.6, ls="--",
                marker=MK_FF, ms=5, label="FF")
    axeb.set_title(r"$\Delta C_\ell^{EB}$")
    axeb.text(0.5, 0.5,
              "parity null\n"
              r"FK:  $\Delta C_\ell^{EB}=0$ by parity" + "\n"
              "FF: numerical residual",
              transform=axeb.transAxes, ha="center", va="center", fontsize=14,
              bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    for ax in axes:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(ELL.min(), ELL.max()); ax.set_ylim(*YLIM)
        ax.set_xlabel(r"$\ell$")
    axes[0].set_ylabel(r"$\ell(\ell+1)\,\Delta C_\ell/2\pi$")

    if handles is not None:
        fig.legend(handles[0], handles[1], loc="upper center", ncol=2,
                   bbox_to_anchor=(0.5, 1.02), frameon=False, fontsize=18)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(OUT_STEM.with_suffix(ext), bbox_inches="tight")
    plt.close(fig)
    plt.rcParams.update(saved)
    print(f"-> {OUT_STEM.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
