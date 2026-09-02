"""New figure: the full convergence 2-point decomposition xi_kk(gamma) into its
O0 / FF / FK channels, with the direct Monte-Carlo overlay for the FF channel
(replaces the archived apples/oranges 'mc_vs_analysis3_xi_kappa').

Lines  = analysis-3 channels (O0, FF full moment, FK).
Markers = Monte-Carlo FF full moment (simulate_ff_crn F-toggle).
The FULL FF moment is compared directly (no connected/disconnected split). O0 is the
dominant, separately-validated linear piece (anchor 0.881); its raw MC has a
large-gamma variance floor, so it is shown as the analytic line only.

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
the caller supplies g_fk/fk_mc/fk_se; `_compute` does not produce them, so the
paper figure is built through
`analyses/revision_2026-08/make_val_figure.py`.

Output: figures/xi_kappa_channels.pdf
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

import driver_stats as ds
import sachs_mc_core as mc

_BASE = ds._SACHS_SFT / "sftwick_outputs" / "2PCF"
LF = 2313.0288751857356
# FF Monte-Carlo sample points (on the analytic [0.5',5000'] grid nodes),
# spanning Figure 11's displayed range; the FF comparison stays clean over the
# whole grid.
GAMMA_MC = (1.0, 2.6, 6.7, 17.3, 44.4, 114.3, 293.9,
            372.2, 471.3, 596.9, 755.9, 957.2, 1212.2, 1535.2, 1944.1)


def _load_kk(name: str):
    d = np.load(_BASE / name / f"xi_{name}.npz", allow_pickle=True)  # trusted
    m = (d["a"] == 0) & (d["b"] == 0)
    xs, ys = d["x"][m], d["y"][m]
    vs = np.asarray(d["value"], float)[m]
    g = []
    for xi, yi in zip(xs, ys):
        xi = np.asarray(xi, float); yi = np.asarray(yi, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)), -1, 1)))) * 60)
    o = np.argsort(g)
    return np.asarray(g)[o], vs[o]


def _cache_path() -> Path:
    p = Path(__file__).resolve().parent / "outputs"
    p.mkdir(parents=True, exist_ok=True)
    return p / "appendix_mc_curve.npz"


def _compute(ff_real: int, ff_nl: int) -> dict:
    """Run the FF Monte-Carlo (deterministic given the fixed seed) plus the
    analysis-3 analytic channels, and cache the result so the figure can be
    restyled without recomputation (see ``--from-cache``)."""
    g, o0 = _load_kk("C_corr_op_O0")
    _, ff_mom = _load_kk("C_corr_op_K_limber_FF")
    _, fk3 = _load_kk("C_corr_op_K_limber_FK")

    ffc_mc = np.zeros(len(GAMMA_MC)); ffc_se = np.zeros(len(GAMMA_MC))
    for i, gv in enumerate(GAMMA_MC):
        cfg_ff = mc.MCConfig(n_real=ff_real, batch_size=3_000, n_lambda=ff_nl,
                             use_f_vertex=True, apply_anchor=False, seed=11)
        r = mc.simulate_ff_crn(cfg_ff, gv)
        ffc_mc[i] = r.ff_moment[0, 0]; ffc_se[i] = r.ff_moment_err[0, 0]
        print(f"  gamma={gv:7.2f}: FF_mc(full)={ffc_mc[i]:.3e} "
              f"+/-{ffc_se[i]:.1e}", flush=True)

    D = dict(g=g, o0=o0, ff_mom=ff_mom, fk3=fk3, g_mc=np.array(GAMMA_MC),
             ffc_mc=ffc_mc, ffc_se=ffc_se)
    np.savez(_cache_path(), **D)
    print(f"[cache -> {_cache_path()}]")
    return D


def _plot(D: dict) -> None:
    """Single-panel, house-style figure: the analysis-3 analytic channels
    (O0/FF/FK, sign-aware lines) overlaid with the direct stochastic-Sachs Monte
    Carlo of the FF channel (markers)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _plot_style import (PALETTE, annotate_sign_legend, apply_rcparams,  # noqa: E402
                             clip_gamma_axis, plot_signed_line)
    apply_rcparams()
    import matplotlib.pyplot as plt

    g = np.asarray(D["g"], float); o0 = np.asarray(D["o0"], float)
    ff_mom = np.asarray(D["ff_mom"], float); fk3 = np.asarray(D["fk3"], float)
    g_mc = np.asarray(D["g_mc"], float)
    ffc_mc = np.asarray(D["ffc_mc"], float); ffc_se = np.asarray(D["ffc_se"], float)
    C_O0, C_FF, C_FK = PALETTE[0], PALETTE[1], PALETTE[2]

    fig, ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)

    # analysis-3 analytic channels: sign-aware coloured lines.  O0 changes sign at
    # large gamma (hollow sign markers there; see the filled/hollow note); the sign
    # markers are kept small so the positive FF/FK lines read cleanly.
    _smk = dict(sign_marker_size=3.2, sign_marker_alpha=0.45)
    plot_signed_line(ax, g, o0, color=C_O0, lw=1.7, alpha=0.95, label="O0 (workflow)", **_smk)
    plot_signed_line(ax, g, ff_mom, color=C_FF, lw=1.7, alpha=0.95, label="FF (workflow)", **_smk)
    # Faint beyond ~1 degree: the FK multipole sum stops converging there.
    plot_signed_line(ax, g, fk3, color=C_FK, lw=1.7, alpha=0.95,
                     label="FK (workflow)", faint_outside=(None, 60.0), **_smk)

    # direct Monte-Carlo overlay (positive in range -> filled markers, matching colour)
    ax.errorbar(g_mc, np.abs(ffc_mc), yerr=ffc_se, fmt="s", color=C_FF, ms=6.5,
                mfc=C_FF, mec=C_FF, capsize=2.5, lw=1.1, zorder=5, label="FF (MC)")
    # FK overlay, only where the analytic FK multipole sum has converged.  The
    # error bars are sub-percent and therefore smaller than the symbols; the
    # quoted agreement lives in the appendix text, not in the plot.
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


def main(ff_real: int, ff_nl: int, from_cache: bool) -> None:
    if from_cache and _cache_path().exists():
        D = dict(np.load(_cache_path(), allow_pickle=True))  # our own cache, safe
        print(f"[from cache {_cache_path()}]")
    else:
        D = _compute(ff_real, ff_nl)
    _plot(D)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ff_real", type=int, default=48_000)
    ap.add_argument("--ff_nl", type=int, default=4000)
    ap.add_argument("--from-cache", action="store_true",
                    help="skip the MC and re-plot from outputs/appendix_mc_curve.npz "
                         "(for fast figure-style iteration)")
    a = ap.parse_args()
    main(a.ff_real, a.ff_nl, a.from_cache)
