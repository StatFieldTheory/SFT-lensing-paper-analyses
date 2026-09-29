"""Paper figure: the FK amplitude as a function of the multipole cutoff.

One panel: the folded FK contribution to xi_kappa at gamma = 0.5', as a
percentage of Order-0, against the vertex-table multipole cutoff ell_max,
for the tree-level and BiHalofit bispectra. The selected products determine
the amplitudes and cutoff dependence. Active-manifest products are hash
verified. Uses the same house style as the other validation figure.

Usage::

    <PyCCL python> fig_cutoff_paper.py            # defaults: the cutoff_ladder run folder, outputs/
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[5]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
_LADDER = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "sftwick_outputs"
           / "2PCF" / "cutoff_ladder")
sys.path.insert(0, str(_REPO / "SFT-lensing-paper-analyses" / "reproduce"))
from product_paths import resolve_product  # noqa: E402

_LADDER = resolve_product("cutoff_ladder", _LADDER)

CUTS = (960, 1920, 3840, 7680, 15360)


def o0_at_half_arcmin() -> float:
    base = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "sftwick_outputs" / "2PCF"
    path = resolve_product("order0", base / "C_corr_op_O0" / "xi_C_corr_op_O0.npz")
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return float(np.interp(0.5, g[o], v[o]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--products", type=Path, default=_LADDER,
                    help="folder holding the folded cutoff ladder")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent / "outputs" / "fk_cutoff_convergence.pdf")
    a = ap.parse_args()

    o0 = o0_at_half_arcmin()
    series = {}
    for model in ("tree", "bihalofit"):
        vals = []
        for cut in CUTS:
            path = a.products / f"table_{model}_cut{cut}_r4_xi.npz"
            with np.load(path, allow_pickle=True) as data:
                mask = (data["a"] == 0) & (data["b"] == 0) & (data["order"] == 2)
                angles = []
                for x, y in zip(data["x"][mask], data["y"][mask], strict=True):
                    x, y = np.asarray(x, float), np.asarray(y, float)
                    cosine = np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y))
                    angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))) * 60)
                match = np.flatnonzero(np.isclose(angles, 0.5, rtol=0, atol=1e-6))
                if len(match) != 1:
                    raise ValueError(f"Expected one convergence sample at 0.5 arcmin: {path}")
                value = float(data["value"][mask][match[0]])
                if not np.isfinite(value):
                    raise ValueError(f"Nonfinite convergence value: {path}")
                vals.append(100.0 * value / o0)
        series[model] = np.array(vals)
        print(f"[{model}] " + "  ".join(f"{c}:{v:.2f}%"
                                        for c, v in zip(CUTS, vals)))

    sys.path.insert(0, str(_MC))
    from _plot_style import PALETTE, apply_rcparams
    apply_rcparams()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
    cuts = np.array(CUTS, float)
    ax.plot(cuts, series["tree"], "o-", color=PALETTE[2], lw=1.7, ms=6.5,
            label="tree-level bispectrum")
    ax.plot(cuts, series["bihalofit"], "s--", color=PALETTE[3], lw=1.5, ms=6.0,
            mfc="none", label="BiHalofit (uncontrolled)")
    ax.axhline(series["tree"][-1], color=PALETTE[2], lw=0.8, ls=":", alpha=0.6)
    ax.annotate(f"{series['tree'][-1]:.2f}%",
                xy=(cuts[0] * 1.05, series["tree"][-1] * 1.12),
                color=PALETTE[2], fontsize=10)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\ell_{\max}$")
    ax.set_ylabel(r"FK contribution to $\xi_\kappa(0.5')$  [% of Order-0]")
    ax.set_xticks(cuts)
    ax.set_xticklabels([str(c) for c in CUTS])
    ax.minorticks_off()
    ax.legend(loc="upper left")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out)
    fig.savefig(a.out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    print(f"[fig -> {a.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
