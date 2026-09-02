"""Figure and table: the deployed FK curve against the corrected one.

Everything except the cosine sampling of the kappa3 vertex is held fixed.
Same shells, same cosmology, same bispectrum, same multipole cutoff, same
quadrature windows, same fold. So the difference between the two curves is
the vertex sampling and nothing else.

Order-0 is drawn for scale because the draft's claims are about where FK
stands relative to it, not about FK in isolation.

Usage::

    python fig_grid_comparison.py --deployed <xi.npz> --corrected <xi.npz> \
        --o0 <xi.npz> --out ../../figures/corrected/grid_comparison.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle
from observables import OBSERVABLES, load_observables

#: Below this magnitude the FK combinations are cancellation noise.
FLOOR = 1e-15


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deployed", type=Path, required=True)
    ap.add_argument("--corrected", type=Path, required=True)
    ap.add_argument("--o0", type=Path, default=None)
    ap.add_argument("--ff", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--table-out", type=Path, default=None)
    args = ap.parse_args()

    g_dep, dep = load_observables(args.deployed, order=2)
    g_cor, cor = load_observables(args.corrected, order=2)
    if not np.allclose(g_dep, g_cor, rtol=1e-6):
        raise SystemExit("the two runs use different separation grids")
    gamma = g_dep
    o0 = load_observables(args.o0, order=0)[1] if args.o0 else None
    ff = load_observables(args.ff, order=2)[1] if args.ff else None

    plotstyle.apply()
    import matplotlib.pyplot as plt

    names = list(OBSERVABLES)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True)
    for i, name in enumerate(names):
        ax = axes[i // 2, i % 2]
        if o0 is not None and o0[name] is not None:
            plotstyle.signed_loglog(ax, gamma, o0[name], label="Order 0",
                                    color=plotstyle.PALETTE[5], marker="s",
                                    linestyle="--", floor=FLOOR)
        if ff is not None and ff[name] is not None:
            plotstyle.signed_loglog(ax, gamma, ff[name], label="FF",
                                    color=plotstyle.PALETTE[2], marker="^",
                                    linestyle=":", floor=FLOOR)
        plotstyle.signed_loglog(ax, gamma, dep[name], label="FK, deployed grid",
                                color=plotstyle.PALETTE[1], marker="o", floor=FLOOR)
        plotstyle.signed_loglog(ax, gamma, cor[name], label="FK, resolved grid",
                                color=plotstyle.PALETTE[0], marker="D", floor=FLOOR)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(plotstyle.OBSERVABLE_TITLES[name])
        if i // 2 == 1:
            ax.set_xlabel(r"separation $\gamma$ [arcmin]")
        if i % 2 == 0:
            ax.set_ylabel("magnitude")
        if i == 0:
            plotstyle.sign_note(ax, loc="upper right")
            ax.legend(loc="lower left")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[fig] -> {args.out}")

    lines = ["gamma[']      FK deployed    FK resolved      ratio"
             + ("        O0        FK_res/O0" if o0 else "")]
    for j, g in enumerate(gamma):
        d, c = dep["xi_kappa"][j], cor["xi_kappa"][j]
        row = f"{g:9.2f}   {d:+.5e}   {c:+.5e}   {d / c if c else np.nan:9.2f}"
        if o0 is not None:
            base = o0["xi_kappa"][j]
            row += f"   {base:+.4e}   {c / base if base else np.nan:+.5f}"
        lines.append(row)
    text = "\n".join(lines)
    print(text)
    if args.table_out:
        args.table_out.write_text(text + "\n")
        print(f"[fig] table -> {args.table_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
