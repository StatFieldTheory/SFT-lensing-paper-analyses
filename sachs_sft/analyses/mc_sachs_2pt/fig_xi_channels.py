"""Convergence reference curves and matched stochastic-Sachs Monte Carlo checks.

Lines use the selected O0 and FK products. The FF line is the matched
local-covariance reference, as its legend states. The curves and markers are
read from the cache selected by reproduce/active_products.json; the isolated
Monte Carlo adapters that build it live in analyses/r1_sft061/mc_update/.
Markers show direct stochastic-Sachs estimates, including the FF full moment.
The FULL FF moment is compared directly (no connected/disconnected split).
O0 is shown as the independently calculated linear reference. The active
revision applies no amplitude anchor.

FK MARKERS: withdrawn 2026-08-26, restored 2026-08-28
=====================================================
They were withdrawn because an exact second-order expectation analysis of the
variance-reduced estimator (`simulate_fk_vr`) showed it measures a placement
SHARE of the vertex (0.211 at 1', crossing zero between 8' and 12') plus a
finite-sigma_lambda term; the apparent agreement at the adopted sigma_lambda=8
was those two errors summing through unity at the tuned value.  That diagnosis
stands, and `simulate_fk_vr` must never be used for FK markers again.

They are restored from a different estimator:
`analyses/mc_fk_complete/`, which is placement-complete (all
three legs of the deformation reach the vertex; proved as a symbolic identity)
and variance-controlled, and whose sigma_lambda -> 0 extrapolation reproduces
the folded FK channel to about a percent.  `_plot` draws FK markers only when
the cache supplies g_fk/fk_mc/fk_se. The cache is built by
`analyses/r1_sft061/mc_update/`; this script only draws it, and `--from-cache`
is accepted but no longer changes anything.

Output: figures/xi_kappa_channels.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "reproduce"))
from product_paths import resolve_product  # noqa: E402


def _cache_path() -> Path:
    return resolve_product("mc_cache")


def _plot(D: dict) -> None:
    """Plot the selected reference curves and their matched Monte Carlo checks."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _plot_style import (PALETTE, annotate_sign_legend, apply_rcparams,  # noqa: E402
                             clip_gamma_axis, plot_signed_line)
    apply_rcparams()
    import matplotlib.pyplot as plt

    g = np.asarray(D["g"], float); o0 = np.asarray(D["o0"], float)
    ff_mom = np.asarray(D["ff_mom"], float); fk3 = np.asarray(D["fk3"], float)
    g_ff = np.asarray(D.get("g_ff", g), float)
    ff_label = str(D.get("ff_reference_label", "FF (workflow)"))
    g_mc = np.asarray(D["g_mc"], float)
    ffc_mc = np.asarray(D["ffc_mc"], float); ffc_se = np.asarray(D["ffc_se"], float)
    C_O0, C_FF, C_FK = PALETTE[0], PALETTE[1], PALETTE[2]

    fig, ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)

    # Selected reference curves: sign-aware coloured lines. O0 changes sign at
    # large gamma (hollow sign markers there; see the filled/hollow note); the sign
    # markers are kept small so the positive FF/FK lines read cleanly.
    _smk = dict(sign_marker_size=3.2, sign_marker_alpha=0.45)
    plot_signed_line(ax, g, o0, color=C_O0, lw=1.7, alpha=0.95, label="O0 (workflow)", **_smk)
    plot_signed_line(ax, g_ff, ff_mom, color=C_FF, lw=1.7, alpha=0.95, label=ff_label, **_smk)
    # Faint beyond ~1 degree: the FK multipole sum stops converging there.
    plot_signed_line(ax, g, fk3, color=C_FK, lw=1.7, alpha=0.95,
                     label="FK (workflow)", faint_outside=(None, 60.0), **_smk)

    # direct Monte-Carlo overlay (positive in range -> filled markers, matching colour)
    ax.errorbar(g_mc, np.abs(ffc_mc), yerr=ffc_se, fmt="s", color=C_FF, ms=6.5,
                mfc=C_FF, mec=C_FF, capsize=2.5, lw=1.1, zorder=5,
                label="FF (local MC)")
    # Matched FK Monte Carlo estimates and their measured standard errors.
    if D.get("g_fk") is not None and len(np.atleast_1d(D["g_fk"])):
        ax.errorbar(np.asarray(D["g_fk"], float), np.abs(np.asarray(D["fk_mc"], float)),
                    yerr=np.asarray(D["fk_se"], float), fmt="o", color=C_FK, ms=6.5,
                    mfc=C_FK, mec=C_FK, capsize=2.5, lw=1.1, zorder=6, label="FK (MC)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\gamma\;[\mathrm{arcmin}]$")
    ax.set_ylabel(r"$\xi_\kappa(\gamma)\;:\;\langle\kappa\kappa\rangle$")
    clip_gamma_axis(ax, gamma_min=float(g.min()))  # default GAMMA_X_MAX_ARCMIN=2000', matches Fig. 11
    ax.legend(loc="lower left", ncol=1)
    annotate_sign_legend(ax)

    out_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "xi_kappa_channels.pdf")
    fig.savefig(out_dir / "xi_kappa_channels.png", dpi=200)
    plt.close(fig)
    print(f"[fig -> {out_dir / 'xi_kappa_channels.pdf'}]")


def main() -> None:
    D = dict(np.load(_cache_path(), allow_pickle=True))  # our own cache, safe
    print(f"[from cache {_cache_path()}]")
    _plot(D)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from-cache", action="store_true",
                    help="accepted for compatibility; the figure is always drawn "
                         "from the cache selected by reproduce/active_products.json")
    ap.parse_args()
    main()
