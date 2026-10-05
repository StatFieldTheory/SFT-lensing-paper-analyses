"""Raw-channel ratio diagnostic, not an accepted FK E/B decomposition.

Plot |zeta_TPP|/|zeta_Bmod| for the explicitly selected collapsed slot0
family. The raw channel ratio is descriptive; it does not validate current
cyclic-placement tensor reconstruction, folded xi_plus/xi_minus or equal
E/B power. No accepted fold is evaluated by this script.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import plotstyle


def collapsed_slot0_rows(triples):
    """Select slot0 with exact 12-decimal keys, refusing duplicate rows."""
    t = np.asarray(triples, dtype=float)
    if t.ndim != 2 or t.shape[1] != 3 or not np.isfinite(t).all():
        raise ValueError("Expected finite cosine triples (n,3)")
    keys = np.round(t, 12)
    if len(np.unique(keys, axis=0)) != len(keys):
        raise ValueError("Duplicate 12-decimal geometry keys")
    mask = (keys[:, 0] == 1.0) & (keys[:, 1] == keys[:, 2])
    idx = np.flatnonzero(mask)
    if not idx.size:
        raise ValueError("No collapsed slot0 rows")
    if len(np.unique(keys[idx, 1])) != len(idx):
        raise ValueError("Duplicate slot0 separations")
    return idx


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
    i = collapsed_slot0_rows(t)
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
            label=r"$|\zeta_{TPP}| / |\zeta_{Bmod}|$  (raw-channel diagnostic)")

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
    ax.set_ylabel("raw-channel ratio |TPP| / |Bmod|")
    ax.set_title("raw-channel diagnostic, not folded E/B acceptance")
    ax.legend(loc="upper left")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight")
    print(f"[raw-channel diagnostic, not accepted FK E/B] -> {args.out}")
    for level in (0.01, 0.1, 0.5, 1.0):
        hit = np.flatnonzero(ratio > level)
        if hit.size:
            print(f"  ratio first exceeds {level:5.2f} at gamma = {g[hit[0]]:.1f}'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
