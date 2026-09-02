"""What the callable corrections change, observable by observable.

Compares two folds that differ only in the vertex callable: the production
one, which sorts the query triple and symmetrises the coupling tensor, and
the permutation-aware one, which evaluates each tensor entry at the
geometry its own field assignment requires.

Both folds must use the same table, the same cutoff and the same geometry,
so the difference is the callable and nothing else.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rebuild"))
from observables import OBSERVABLES, load_observables  # noqa: E402

FLOOR = 1e-15
O0 = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
      "sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--legacy", type=Path, required=True,
                    help="fold with the production callable")
    ap.add_argument("--fixed", type=Path, required=True,
                    help="fold with the permutation-aware callable")
    ap.add_argument("--gamma-max", type=float, default=200.0)
    ap.add_argument("--tex-out", type=Path, default=None)
    args = ap.parse_args()

    g, legacy = load_observables(args.legacy, order=2)
    g2, fixed = load_observables(args.fixed, order=2)
    if not np.allclose(g, g2, rtol=1e-6):
        raise SystemExit("the two folds use different separation grids")
    _, base = load_observables(Path(O0), order=0)

    rows = []
    print("Effect of the callable corrections, at fixed table and geometry\n")
    for name in OBSERVABLES:
        a, b = legacy[name], fixed[name]
        if a is None or b is None:
            continue
        keep = (np.abs(a) > FLOOR) & (g <= args.gamma_max)
        peak = int(np.argmax(np.abs(b)))
        ratio = b[keep] / a[keep]
        print(f"{name}")
        print(f"   production callable, peak  {a[int(np.argmax(np.abs(a)))]:+.4e} "
              f"at {g[int(np.argmax(np.abs(a)))]:7.2f}'")
        print(f"   corrected  callable, peak  {b[peak]:+.4e} at {g[peak]:7.2f}'"
              f"   ({100*abs(b[peak]/base[name][peak]):.2f}% of Order-0)")
        if ratio.size:
            print(f"   corrected / production: median {np.median(ratio):+.4f}, "
                  f"range [{ratio.min():+.4f}, {ratio.max():+.4f}]")
        print(f"   at 0.5': {a[0]:+.5e} -> {b[0]:+.5e} "
              f"(change {100*(b[0]/a[0]-1):+.2f}%)" if a[0] else "")
        print()
        rows.append((name, a[0], b[0],
                     float(np.median(ratio)) if ratio.size else np.nan,
                     b[peak], g[peak]))

    if args.tex_out:
        label = {"xi_kappa": r"$\xi_\kappa$", "xi_plus": r"$\xi_+$",
                 "xi_minus": r"$\xi_-$",
                 "xi_kappa_gamma_t": r"$\xi_{\kappa\gamma_t}$"}
        lines = [r"\begin{table}[h]", r"\centering", r"\begin{tabular}{lrrrr}",
                 r"\toprule",
                 r"observable & at $0.5'$, production & at $0.5'$, corrected & "
                 r"median ratio & corrected peak \\", r"\midrule"]
        for name, a0, b0, med, pk, pg in rows:
            lines.append(
                f"{label.get(name, name)} & ${a0:+.4e}$ & ${b0:+.4e}$ & "
                f"{med:+.3f} & ${pk:+.2e}$ at ${pg:.0f}'$ \\\\".replace(
                    "e-0", r"\times10^{-").replace("e+0", r"\times10^{"))
        lines += [r"\bottomrule", r"\end{tabular}",
                  r"\caption{Effect of the two callable corrections at fixed "
                  r"table, cutoff and geometry. Produced by "
                  r"\code{code/callable\_fixed/compare\_callables.py}.}",
                  r"\label{tab:callablefix}", r"\end{table}"]
        args.tex_out.write_text("\n".join(lines) + "\n")
        print(f"-> {args.tex_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
