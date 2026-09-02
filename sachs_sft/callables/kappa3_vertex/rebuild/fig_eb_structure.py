"""Figure: the channel that carries E minus B, against separation.

The paper's equal-split claim, Delta C_EE = Delta C_BB, is the statement
that one vertex channel vanishes. That channel is zeta_TPP, the real part
of <Phi_00 Psi_0 Psi_0>, which in real components is
<Phi_00 (Psi_+^2 - Psi_x^2)>. Its partner in the sum, which carries
E plus B, is zeta_Bmod = <Phi_00 (Psi_+^2 + Psi_x^2)>.

Plotting their ratio shows where the claim holds and where it does not.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--shell", type=int, default=-1)
    ap.add_argument("--gamma-max", type=float, default=60.0,
                    help="beyond this the E+B channel passes through zero and "
                         "the ratio stops being meaningful")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    d = np.load(args.table, allow_pickle=False)
    t = np.asarray(d["cosine_triples"], float)
    m = np.isclose(t[:, 0], 1.0) & np.isclose(t[:, 1], t[:, 2])
    i = np.flatnonzero(m)
    gamma = np.degrees(np.arccos(np.clip(t[i, 1], -1, 1))) * 60.0
    order = np.argsort(gamma)
    i, gamma = i[order], gamma[order]
    tpp = np.asarray(d["zeta_TPP"], float)[i, args.shell]
    bmod = np.asarray(d["zeta_Bmod"], float)[i, args.shell]

    keep = (gamma > 0) & (gamma <= args.gamma_max) & (bmod != 0)
    g, ratio = gamma[keep], np.abs(tpp[keep] / bmod[keep])

    plotstyle.apply()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.plot(g, ratio, "-o", color=plotstyle.PALETTE[0], markersize=4.2,
            label=r"$|\zeta_{TPP}| / |\zeta_{Bmod}|$  (that is $|E-B|/|E+B|$)")

    anchor = np.argmin(np.abs(g - 2.0))
    reference = ratio[anchor] * (g / g[anchor]) ** 4
    ax.plot(g, reference, "--", color=plotstyle.PALETTE[1], linewidth=1.3,
            label=r"$\gamma^4$, matched at $\gamma=2'$")

    for level, style in ((0.01, ":"), (0.1, "-."), (1.0, "-")):
        hit = np.flatnonzero(ratio > level)
        if hit.size:
            ax.axvline(g[hit[0]], color="0.55", linestyle=style, linewidth=1.0)
            ax.text(g[hit[0]] * 1.03, 2e-6,
                    f"  {level:g} at {g[hit[0]]:.0f}'", rotation=90,
                    fontsize=8.5, color="0.35", va="bottom")

    ax.axhline(1.0, color="0.3", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"separation $\gamma$ [arcmin]")
    ax.set_ylabel("ratio of the difference channel to the sum channel")
    ax.set_title("where the equal-split claim holds")
    ax.legend(loc="upper left")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[fig] -> {args.out}")
    for level in (0.01, 0.1, 0.5, 1.0):
        hit = np.flatnonzero(ratio > level)
        if hit.size:
            print(f"  ratio first exceeds {level:5.2f} at gamma = {g[hit[0]]:.1f}'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
