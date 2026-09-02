"""Angular power spectra of the FK variants, in the draft's own convention.

Reuses `analysis3/plot_analysis3_cl_decomposition.py` for the transform, so
this is the same curved-sky Wigner-d projection the draft's C_ell figure
uses (kappa-kappa goes through d^l_{0,0} = P_l), with the same apodisation
and DC subtraction. Only the input curves differ.

Usage (from ``driver_field_emulators/code``)::

    SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
    $SFT plot_fk_variants_cl.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSES = (HERE.parents[1] / "SFT-lensing-paper-analyses" / "sachs_sft"
            / "analyses" / "analysis3")
sys.path.insert(0, str(ANALYSES))

from _plot_style import MARKERS, PALETTE, apply_rcparams  # noqa: E402
from plot_analysis3_cl_decomposition import (  # noqa: E402
    ELL,
    build_curved_matrix,
    forward_curved,
)

from compare_fk_variants import load_curve  # noqa: E402
from plot_fk_variants import O0_NPZ, PRODUCTS, VARIANTS, FIGURES  # noqa: E402


def main() -> None:
    apply_rcparams()
    import matplotlib.pyplot as plt

    gamma, o0 = load_curve(O0_NPZ, order=0)
    # kappa-kappa is the spin-0 channel: d^l_{0,0} = P_l.
    setup = build_curved_matrix(gamma, ELL, 0, 0)
    cl_o0 = forward_curved(o0, setup)

    curves = []
    for label, path in VARIANTS:
        g, v = load_curve(path, order=2)
        if not np.allclose(g, gamma):
            raise ValueError(f"{path.name} gamma grid differs")
        curves.append((label, forward_curved(v, setup)))

    band = ELL * (ELL + 1) / (2.0 * np.pi)
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.2, 8.6), sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0]})

    ax_top.plot(ELL, band * np.abs(cl_o0), color="0.30", lw=3.2, alpha=0.35,
                label="Order-0 (reference)", zorder=1)
    for i, (label, cl) in enumerate(curves):
        pos = cl > 0
        ax_top.plot(ELL[pos], (band * np.abs(cl))[pos], ls="none",
                    marker=MARKERS[i], ms=5.0, color=PALETTE[i],
                    markerfacecolor=PALETTE[i], label=label)
        ax_top.plot(ELL[~pos], (band * np.abs(cl))[~pos], ls="none",
                    marker=MARKERS[i], ms=5.0, color=PALETTE[i],
                    markerfacecolor="none")
    ax_top.set_xscale("log"); ax_top.set_yscale("log")
    ax_top.set_ylabel(r"$\ell(\ell+1)\,|C_\ell^{\kappa\kappa}|\,/\,2\pi$")
    ax_top.legend(loc="lower left", frameon=False, fontsize=10.5)
    ax_top.text(0.98, 0.98, "filled: $+$   hollow: $-$", ha="right", va="top",
                transform=ax_top.transAxes, fontsize=10.5, color="0.35")

    for i, (label, cl) in enumerate(curves):
        ax_bot.plot(ELL, np.abs(cl / cl_o0), color=PALETTE[i],
                    marker=MARKERS[i], ms=4.5, lw=1.0, label=label)
    ax_bot.axhline(1.0, color="0.35", lw=1.0, ls=":")
    ax_bot.set_xscale("log"); ax_bot.set_yscale("log")
    ax_bot.set_ylim(1e-4, 1e2)
    ax_bot.set_xlabel(r"$\ell$")
    ax_bot.set_ylabel(r"$|C_\ell^{\rm FK}\,/\,C_\ell^{\rm Order\text{-}0}|$")

    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".pdf"):
        fig.savefig(FIGURES / f"fk_variants_cl{ext}", dpi=180)
    plt.close(fig)
    print(f"-> {FIGURES / 'fk_variants_cl.pdf'}")

    print("\n  ell    C_l O0      " + "  ".join(f"{l[:16]:>16s}" for l, _ in curves))
    for j in range(0, len(ELL), 4):
        row = "  ".join(f"{abs(c[j] / cl_o0[j]):16.4f}" for _, c in curves)
        print(f"  {int(ELL[j]):5d}  {cl_o0[j]:+.3e}  " + row)


if __name__ == "__main__":
    main()
