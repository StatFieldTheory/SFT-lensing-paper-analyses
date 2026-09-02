"""Task D: tree against BiHalofit on the corrected grid, by cutoff.

Emits both a printed summary and a LaTeX table body, so the note's numbers
are generated rather than transcribed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from observables import load_observables

O0 = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
      "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")


def value(products: Path, model: str, cutoff: int):
    path = products / f"table_{model}_cut{cutoff}_r4_xi.npz"
    if not path.exists():
        return None
    g, o = load_observables(path, order=2)
    return g, o["xi_kappa"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--products", type=Path, required=True)
    ap.add_argument("--cutoffs", nargs="+", type=int,
                    default=[960, 1920, 3840, 7680, 15360])
    ap.add_argument("--tex-out", type=Path, default=None)
    args = ap.parse_args()

    _, base = load_observables(Path(O0), order=0)
    o0 = base["xi_kappa"][0]

    rows = []
    print(f"{'ell_max':>8} {'tree':>14} {'BiHalofit':>14} {'ratio':>8} "
          f"{'tree/O0':>9} {'BiHalofit/O0':>13}   k_hard at the near shell")
    for cut in args.cutoffs:
        t = value(args.products, "tree", cut)
        b = value(args.products, "bihalofit", cut)
        if t is None:
            continue
        tv = t[1][0]
        k_hard = 2 * cut / 292.6
        if b is None:
            print(f"{cut:8d} {tv:+14.5e} {'(pending)':>14} {'':>8} "
                  f"{100*tv/o0:8.2f}% {'':>13}   {k_hard:6.1f} h/Mpc")
            rows.append((cut, tv, None, k_hard))
            continue
        bv = b[1][0]
        print(f"{cut:8d} {tv:+14.5e} {bv:+14.5e} {bv/tv:8.3f} "
              f"{100*tv/o0:8.2f}% {100*bv/o0:12.2f}%   {k_hard:6.1f} h/Mpc")
        rows.append((cut, tv, bv, k_hard))

    if args.tex_out:
        lines = [r"\begin{table}[h]", r"\centering",
                 r"\begin{tabular}{rrrrrr}", r"\toprule",
                 r"$\elmax$ & tree & BiHalofit & ratio & tree/O0 & "
                 r"hard $k$, near shell \\",
                 r"\midrule"]
        for cut, tv, bv, k in rows:
            bcell = (f"${bv/1e-6:.3f}\\times10^{{-6}}$" if bv is not None
                     else "pending")
            rcell = f"{bv/tv:.3f}" if bv is not None else "---"
            lines.append(f"{cut} & ${tv/1e-6:.3f}\\times10^{{-6}}$ & {bcell} & "
                         f"{rcell} & {100*tv/o0:.2f}\\% & "
                         f"{k:.0f}\\,$h$/Mpc \\\\")
        lines += [r"\bottomrule", r"\end{tabular}",
                  r"\caption{The FK contribution to $\xi_\kappa$ at "
                  r"$\gamma=0.5'$, tree level against BiHalofit, on the "
                  r"corrected grid, as a function of the multipole cutoff. "
                  r"Order-0 there is $8.443\times10^{-4}$. The last column is "
                  r"the largest wavenumber the bispectrum is asked for at the "
                  r"innermost shell; BiHalofit is calibrated to "
                  r"$k\lesssim10\,h$/Mpc. Produced by "
                  r"\code{code/rebuild/table\_nonlinear.py}.}",
                  r"\label{tab:nonlinear}", r"\end{table}"]
        args.tex_out.write_text("\n".join(lines) + "\n")
        print(f"\n-> {args.tex_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
