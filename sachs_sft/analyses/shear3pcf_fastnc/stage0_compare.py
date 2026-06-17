"""STAGE 0 comparison + plot: ours (zeta-fold) vs SPT-Limber reference.

Interpreter (ABSOLUTE):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python

Loads outputs/stage0_ours.npz + outputs/stage0_spt_reference.npz, computes the
|ratio| over the gamma sweep, writes outputs/stage0_results.npz +
outputs/stage0_compare.png + figures/stage0_convergence_3pcf.pdf.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
_OUT = _HERE / "outputs"
_FIG = _HERE / "figures"


def main() -> int:
    ours = np.load(_OUT / "stage0_ours.npz")
    ref = np.load(_OUT / "stage0_spt_reference.npz")

    g = ours["gamma_arcmin"]
    assert np.allclose(g, ref["gamma_arcmin"]), "gamma grids differ"
    Z_ours = ours["Z_lambda"]          # lambda-route fold (primary)
    Z_ours_chi = ours["Z_chi"]         # chi-route cross-check
    Z_ref = ref["Z_ref"]

    ratio = np.abs(Z_ours) / np.abs(Z_ref)
    sign_ours = np.sign(Z_ours)
    sign_ref = np.sign(Z_ref)

    print("gamma['] |   ours (Z_lambda) |   ref (Z_ref)   | |ratio| | sgn(ours) sgn(ref)")
    for i in range(g.size):
        print(f"{g[i]:8.3f} | {Z_ours[i]:+.5e} | {Z_ref[i]:+.5e} | "
              f"{ratio[i]:6.3f} |   {int(sign_ours[i]):+d}      {int(sign_ref[i]):+d}")

    i1 = int(np.argmin(np.abs(g - 1.0)))
    print()
    print(f"SINGLE-GAMMA gamma=1': ours={Z_ours[i1]:+.5e}  ref={Z_ref[i1]:+.5e}  "
          f"|ratio|={ratio[i1]:.4f}  signs=({int(sign_ours[i1]):+d},{int(sign_ref[i1]):+d})")
    # sweep ratio over the small-angle regime where both are well-resolved
    small = g <= 60.0
    print(f"SWEEP |ratio| range (gamma<=60'): [{ratio[small].min():.3f}, "
          f"{ratio[small].max():.3f}]")
    print(f"SWEEP |ratio| range (all gamma):  [{ratio.min():.3f}, {ratio.max():.3f}]")

    np.savez(
        _OUT / "stage0_results.npz",
        gamma_arcmin=g, Z_ours=Z_ours, Z_ours_chi=Z_ours_chi, Z_ref=Z_ref,
        ratio_abs=ratio, sign_ours=sign_ours, sign_ref=sign_ref,
    )

    # --- plot: |Z| log-log + ratio panel --------------------------------------
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(7.0, 6.6), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1.4], hspace=0.08))

    ax0.plot(g, np.abs(Z_ours), "o-", color="C0", lw=1.8, ms=6,
             label=r"ours: $|\int d\lambda\,K^3\,\zeta_{TTT}|$")
    ax0.plot(g, np.abs(Z_ref), "s--", color="C3", lw=1.6, ms=5,
             label=r"ref: SPT tree Limber $|Z_\kappa|$")
    ax0.set_yscale("log")
    ax0.set_xscale("log")
    ax0.set_ylabel(r"$|Z_\kappa(\gamma)|$  (convergence 3PCF)")
    ax0.set_title(r"STAGE 0: scalar convergence 3PCF, squeezed, $z_s=5$ "
                  r"(both $Z_\kappa<0$)")
    ax0.grid(True, which="both", alpha=0.25)
    ax0.legend(frameon=False, fontsize=10)

    ax1.axhline(1.0, color="k", lw=0.8, ls=":")
    ax1.axhspan(0.5, 1.5, color="green", alpha=0.08,
                label=r"$|ratio-1|<0.5$ gate")
    ax1.plot(g, ratio, "o-", color="C2", lw=1.6, ms=5)
    ax1.set_xscale("log")
    ax1.set_ylim(0.0, 2.0)
    ax1.set_xlabel(r"$\gamma$ [arcmin]")
    ax1.set_ylabel(r"$|$ours$/$ref$|$")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(frameon=False, fontsize=9, loc="lower left")

    fig.savefig(_OUT / "stage0_compare.png", dpi=140, bbox_inches="tight")
    _FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG / "stage0_convergence_3pcf.pdf", bbox_inches="tight")
    print(f"\nsaved -> {_OUT/'stage0_compare.png'}")
    print(f"saved -> {_FIG/'stage0_convergence_3pcf.pdf'}")
    print(f"saved -> {_OUT/'stage0_results.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
