"""Figure: amplitude as a function of the multipole cutoff, and in k-space.

Three panels, because the cutoff question has three separate parts.

Left: the vertex itself, as a partial sum over octave bands, normalised to
the deployed cutoff. This says whether the integral has turned over.

Middle: the folded observable at the smallest plotted separation, for each
cutoff and each bispectrum model. This is what the paper would quote.

Right: the wavenumbers the cutoff reaches at the innermost and outermost
shells, against the range where the nonlinear bispectrum fit is calibrated
and where baryons stop being a small correction. This is what turns
"underconverged" into a bounded statement about model error.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle
from observables import load_observables

BIHALOFIT_KMAX = 10.0
BARYON_K = 1.0


def fold_amplitude(products: Path, model: str, cutoffs, gamma_index=0):
    out = []
    for cut in cutoffs:
        path = products / f"table_{model}_cut{cut}_r4_xi.npz"
        if not path.exists():
            out.append((cut, np.nan))
            continue
        _, obs = load_observables(path, order=2)
        out.append((cut, float(obs["xi_kappa"][gamma_index])))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vertex", type=Path, required=True,
                    help="npz from analyse_cutoff.py")
    ap.add_argument("--products", type=Path, required=True)
    ap.add_argument("--models", nargs="+", default=["tree", "bihalofit"])
    ap.add_argument("--cutoffs", nargs="+", type=int,
                    default=[960, 1920, 3840, 7680, 15360])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    d = np.load(args.vertex)
    plotstyle.apply()
    import matplotlib.pyplot as plt
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.4, 4.6))

    for mi, model in enumerate(args.models):
        key = f"{model}_edges"
        if key not in d:
            continue
        edges = d[key]
        curves = d[f"{model}_curves"]
        rows = d[f"{model}_rows"]
        gammas = d["gammas"]
        base = np.argmin(np.abs(edges - 960))
        for gi, g in enumerate(gammas[:3]):
            total = curves[:, rows[gi], -1]
            a.plot(edges, total / total[base],
                   marker="oDs^"[gi % 4], markersize=4.2,
                   linestyle="-" if mi == 0 else "--",
                   color=plotstyle.PALETTE[gi],
                   label=f"{model}, $\\gamma={g:.2g}'$" if mi < 2 else None)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel(r"cutoff $\ell_{\max}$")
    a.set_ylabel(r"$\zeta_{TTT}$ relative to $\ell_{\max}=960$")
    a.set_title("the vertex, partial sum over octaves")
    a.legend(loc="upper left", fontsize=8)

    for mi, model in enumerate(args.models):
        pairs = fold_amplitude(args.products, model, args.cutoffs)
        cut = np.array([p[0] for p in pairs], float)
        val = np.array([p[1] for p in pairs], float)
        keep = np.isfinite(val)
        if not keep.any():
            continue
        b.plot(cut[keep], np.abs(val[keep]), "o-",
               color=plotstyle.PALETTE[mi], markersize=5,
               label=f"{model}")
    b.set_xscale("log")
    b.set_yscale("log")
    b.set_xlabel(r"cutoff $\ell_{\max}$")
    b.set_ylabel(r"$|\xi^{\rm FK}_\kappa(0.5')|$")
    b.set_title("the folded observable")
    b.legend(loc="best")

    model = args.models[0]
    if f"{model}_chi_h" in d:
        chi_h = d[f"{model}_chi_h"]
        edges = d[f"{model}_edges"]
        for label, index, style in (("innermost shell", 0, "-"),
                                    ("source shell", -1, "--")):
            c.plot(edges, 2 * edges / chi_h[index], style, marker="o",
                   markersize=4.2,
                   color=plotstyle.PALETTE[0 if index == 0 else 2],
                   label=f"{label}, hard leg")
        c.axhspan(BIHALOFIT_KMAX, 1e3, color=plotstyle.PALETTE[1], alpha=0.12)
        c.axhline(BIHALOFIT_KMAX, color=plotstyle.PALETTE[1], linewidth=1.1)
        c.text(edges[0] * 1.1, BIHALOFIT_KMAX * 1.3,
               "beyond BiHalofit calibration", fontsize=8.5,
               color=plotstyle.PALETTE[1])
        c.axhline(BARYON_K, color=plotstyle.PALETTE[3], linewidth=1.1,
                  linestyle=":")
        c.text(edges[0] * 1.1, BARYON_K * 1.25, "baryons matter above",
               fontsize=8.5, color=plotstyle.PALETTE[3])
        c.set_ylim(1e-2, 3e2)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlabel(r"cutoff $\ell_{\max}$")
    c.set_ylabel(r"wavenumber $k$ [$h$/Mpc]")
    c.set_title("what the cutoff reaches")
    c.legend(loc="lower right", fontsize=8.5)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[fig] -> {args.out}")

    for model in args.models:
        print(f"\n{model}: folded |xi_FK(0.5')| against cutoff")
        for cut, val in fold_amplitude(args.products, model, args.cutoffs):
            print(f"  {cut:6d}  {val:+.5e}" if np.isfinite(val)
                  else f"  {cut:6d}  (not built)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
