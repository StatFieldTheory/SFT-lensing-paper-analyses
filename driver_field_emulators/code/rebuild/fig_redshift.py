"""The FK term against source redshift, and the cost of the line-of-sight floor.

Two questions share one figure because they are entangled and the draft
states one of them.

Left: the corrected FK contribution relative to Order-0, at each source
redshift. The draft says the leakage "strengthens steadily with source
redshift"; this is where that is tested on a resolved vertex.

Right: the same FK curve computed with the tabulated shells extended toward
the observer, divided by the standard one. That ratio is exactly what the
lambda_min floor costs, and it is shown against source redshift because the
kernel weight of the excised foreground grows as the source moves nearer.

Usage::

    python fig_redshift.py --fk <multiz_xi.npz> --o0 <order0_multiz_xi.npz> \
        --fk-extended <ext_multiz_xi.npz> --out fig.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle
from observables import load_observables, source_distances

FLOOR = 1e-15


def redshift_of(lam: np.ndarray) -> np.ndarray:
    """Source redshift for each affine distance, on the vertex table's mapping."""
    from scipy.interpolate import CubicSpline
    repo = Path(__file__).resolve().parents[3]
    table = (repo / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
             / "kappa3_vertex" / "equal_time_limber"
             / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz")
    d = np.load(table, allow_pickle=False)
    return CubicSpline(d["lambda_shells_Mpc"], d["z_shells"])(lam)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fk", type=Path, required=True)
    ap.add_argument("--o0", type=Path, default=None)
    ap.add_argument("--fk-extended", type=Path, default=None)
    ap.add_argument("--fk-standard", type=Path, default=None,
                    help="denominator for the extended/standard ratio. The "
                         "ratio is only meaningful between two folds that "
                         "share a callable, so this defaults to --fk but must "
                         "be set when --fk-extended was produced with a "
                         "different one.")
    ap.add_argument("--observable", default="xi_kappa")
    ap.add_argument("--gamma-max", type=float, default=200.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    lams = source_distances(args.fk)
    zs = redshift_of(lams)
    plotstyle.apply()
    import matplotlib.pyplot as plt

    n_panels = 2 if args.fk_extended else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(6.2 * n_panels, 4.6),
                             squeeze=False)
    left = axes[0][0]

    print(f"{'z_s':>6} {'lambda_s':>10} " +
          ("  FK/O0 at 0.5'   FK/O0 at 21.9'" if args.o0 else "  FK(0.5')"))
    for i, (lam, z) in enumerate(zip(lams, zs)):
        gamma, fk = load_observables(args.fk, order=2, t_final=float(lam))
        y = fk[args.observable]
        colour = plotstyle.PALETTE[i % len(plotstyle.PALETTE)]
        if args.o0:
            _, o0 = load_observables(args.o0, order=0, t_final=float(lam))
            base = o0[args.observable]
            keep = (np.abs(base) > FLOOR) & (gamma <= args.gamma_max)
            left.plot(gamma[keep], np.abs(y[keep] / base[keep]), "-o",
                      color=colour, markersize=3.6, label=f"$z_s={z:.1f}$")
            j0 = int(np.argmin(np.abs(gamma - 0.5)))
            j1 = int(np.argmin(np.abs(gamma - 21.88)))
            print(f"{z:6.2f} {lam:10.1f}  {y[j0]/base[j0]:14.5f} "
                  f"{y[j1]/base[j1]:16.5f}")
        else:
            plotstyle.signed_loglog(left, gamma, y, label=f"$z_s={z:.1f}$",
                                    color=colour, floor=FLOOR)
            print(f"{z:6.2f} {lam:10.1f}  {y[0]:+.5e}")
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_xlabel(r"separation $\gamma$ [arcmin]")
    left.set_ylabel(r"$|\xi^{\rm FK}_\kappa| / |\xi^{O0}_\kappa|$" if args.o0
                    else "magnitude")
    left.set_title("FK relative to Order-0, by source redshift")
    left.legend(loc="best", fontsize=8.5, ncol=2)

    if args.fk_extended:
        right = axes[0][1]
        print(f"\n{'z_s':>6}  median ratio (extended shells / standard), "
              f"gamma <= {args.gamma_max:.0f}'")
        standard = args.fk_standard or args.fk
        if standard != args.fk:
            print(f"  ratio panel uses {standard.name} as the standard-shell "
                  f"fold, to keep both sides on one callable")
        for i, (lam, z) in enumerate(zip(lams, zs)):
            gamma, a = load_observables(standard, order=2, t_final=float(lam))
            _, b = load_observables(args.fk_extended, order=2, t_final=float(lam))
            ya, yb = a[args.observable], b[args.observable]
            keep = (np.abs(ya) > FLOOR) & (gamma <= args.gamma_max)
            ratio = yb[keep] / ya[keep]
            right.plot(gamma[keep], ratio, "-o",
                       color=plotstyle.PALETTE[i % len(plotstyle.PALETTE)],
                       markersize=3.6, label=f"$z_s={z:.1f}$")
            print(f"{z:6.2f}  {np.median(ratio):.5f}   "
                  f"[{ratio.min():.5f}, {ratio.max():.5f}]")
        right.axhline(1.0, color="0.4", linewidth=0.9, linestyle="--")
        right.set_xscale("log")
        right.set_xlabel(r"separation $\gamma$ [arcmin]")
        right.set_ylabel("extended shells / standard")
        right.set_title(r"cost of the $\lambda_{\min}$ floor")
        right.legend(loc="best", fontsize=8.5, ncol=2)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"\n[fig] -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
