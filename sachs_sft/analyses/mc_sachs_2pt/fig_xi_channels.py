"""New figure: the full convergence 2-point decomposition xi_kk(gamma) into its
O0 / FF / FK channels, with the CORRECT Monte-Carlo overlay (replaces the archived
apples/oranges 'mc_vs_analysis3_xi_kappa').

Lines  = analysis-3 channels (O0, FF full moment, FK).
Markers = Monte-Carlo: FK (simulate_fk_vr), FF full moment (simulate_ff_crn F-toggle).
The FULL FF moment is compared directly (no connected/disconnected split). O0 is the
dominant, separately-validated linear piece (anchor 0.881); its raw MC has a
large-gamma variance floor, so it is shown as the analytic line only.

FK MONTE-CARLO CONVERGENCE (2026-06-10)
=======================================
The FK MC injects the equal-time 3-cumulant through a finite-width AR(1)
colored-noise field of correlation length sigma_lambda (`simulate_fk_vr`).  The
analytic FK is the delta / white-noise (sigma_lambda -> 0) limit.  At a MODERATE,
well-resolved grid (n_lambda=1000, sigma_lambda/dlambda=4.7) decreasing
sigma_lambda drives the FK MC/analytic ratio monotonically toward 1.  Measured at
gamma=10' (6 seeds, n_real=9000/seed; `_fk_converge.py`-style sweep):

    sigma_lambda  18    14    11     9     7     5     4
    sigma/dlambda 9.4   7.3   5.8   4.7   3.7   2.6   2.1
    FK MC/analytic 1.99  1.69  1.35  1.08  0.79  0.52  0.42

i.e. a clean monotone descent through 1, crossing unity near sigma_lambda~9
(sigma/dlambda~4.7).  Below sigma_lambda~7 the ratio undershoots, partly the
genuine sigma->0 over-shoot and partly AR(1) under-resolution (sigma/dlambda<4).
IMPORTANT (counter-intuitive): the convergence knob is sigma_lambda, NOT
n_lambda.  At FIXED sigma_lambda a FINER grid makes the ratio WORSE (the
near-delta Levy zeta has skewness ~4-6, so a finer grid samples rare large
excursions more; n_lambda=2000 inflates gamma=10' to ~15-21x).  So n_lambda is
held at the moderate value and sigma_lambda alone is swept.

ADOPTED FIGURE SETTING: sigma_lambda=8, n_lambda=1000 (sigma/dlambda=4.2, still
AR(1)-resolved), 8 seeds, n_real=24000/seed -- the script default that builds the
published cache `outputs/appendix_mc_curve.npz`.  FK MC/analytic across the gamma
grid:

    gamma'       1.0   2.6   6.7  17.3  44.4 114.3 293.9   median
    FK MC/analytic 1.26 1.18 1.23 1.08 1.12 1.31 1.36    1.23

i.e. FK MC/analytic ~ 1.2-1.25 and roughly flat in gamma (range [1.08, 1.36] over
gamma<300'), down from ~2.4 at the old sigma_lambda=18 -- the variance-reduced
stochastic-Sachs MC reproduces the analytic FK workflow to within ~20-25% (the
finite-sigma_lambda overshoot, matching the paper caption's "~25%").  NOTE the
ratio itself drifts UP with sample count: a lower-statistics 9000/seed run gives
median ~1.12, because the near-delta Levy zeta is heavy-tailed (skewness ~4-6) so
more realisations sample rare large excursions.  The per-gamma inter-seed spread
is ~+/-40%.  This is an order-of-magnitude corroboration, not a converged
precision number; `fk_analytic` (0.992 vs analysis-3) remains THE precise FK
value.  The monotone sigma_lambda-trend above is the demonstration that the
analytic workflow IS the converged (sigma_lambda->0) limit of the direct MC.
See FK_NOTES.md for estimator details.

LARGE-GAMMA FK TAIL (gamma > 294', added 2026-06-13 to match Figure 11's range)
================================================================================
The FK MC/analytic agreement above is validated only on gamma <= 293.9'.  Beyond
that the GAMMA_MC grid now reaches ~1944', and two effects make the FK markers
there NOT a clean validation (FF stays clean, ratio ~0.88, to ~1944'):
  * The analytic FK kk channel is flat ~4.3e-6 only to ~150'; it then decays
    (3.4e-6 @597', 2.4e-6 @957', 9.2e-7 @1535') and CROSSES ZERO near 1944'.
  * The finite-sigma_lambda Levy-zeta systematic, ~1.1-1.35x at <=294', GROWS with
    gamma (measured ~1.6x @597', ~1.9x @1535'), so the FK markers sit above the
    decaying line.  Tightening it would need a per-node sigma_lambda convergence
    study, out of scope for this validation figure.
The FK markers past ~294' are therefore the honest (noisy, systematically high) MC
tail, not a precision check; the precision check lives at gamma <= 294'.

Output: figures/xi_kappa_channels.pdf
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

import driver_stats as ds
import fk_analytic as fka
import sachs_mc_core as mc

_BASE = ds._SACHS_SFT / "sftwick_outputs" / "2PCF"
LF = 2313.0288751857356
# MC sample points (on the analytic [0.5',5000'] grid nodes).  Extended past the
# original 293.9' out to ~1944' so the validation spans Figure 11's displayed
# range.  CAVEAT (see the FK-tail note in the docstring): beyond ~294' the FK MC
# carries a growing finite-sigma_lambda systematic (ratio ~1.6-1.9x at 600-1500')
# and the analytic FK decays through zero near 1944'; the FF MC stays clean (~0.88).
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


def _compute(n_real: int, n_lambda: int, sigma: float, seeds: int,
             ff_real: int, ff_nl: int) -> dict:
    """Run the stochastic-Sachs MC (deterministic given the fixed seeds) plus the
    analysis-3 analytic channels, and cache the result.  The heavy MC is cached so
    the figure can be restyled without recomputation (see ``--from-cache``)."""
    g, o0 = _load_kk("C_corr_op_O0")
    _, ff_mom = _load_kk("C_corr_op_K_limber_FF")
    _, fk3 = _load_kk("C_corr_op_K_limber_FK")

    fk_mc = np.zeros(len(GAMMA_MC)); fk_se = np.zeros(len(GAMMA_MC))
    ffc_mc = np.zeros(len(GAMMA_MC)); ffc_se = np.zeros(len(GAMMA_MC))
    for i, gv in enumerate(GAMMA_MC):
        vs = []
        for sd in range(1, seeds + 1):
            cfg = mc.MCConfig(n_real=n_real, batch_size=3_000, n_lambda=n_lambda,
                              use_f_vertex=True, skew_scale=1.0,
                              sigma_lambda=sigma, seed=sd)
            vs.append(mc.simulate_fk_vr(cfg, gv).fk[0, 0])
        vs = np.array(vs)
        fk_mc[i] = vs.mean(); fk_se[i] = vs.std() / np.sqrt(len(vs))
        cfg_ff = mc.MCConfig(n_real=ff_real, batch_size=3_000, n_lambda=ff_nl,
                             use_f_vertex=True, apply_anchor=False, seed=11)
        r = mc.simulate_ff_crn(cfg_ff, gv)
        ffc_mc[i] = r.ff_moment[0, 0]; ffc_se[i] = r.ff_moment_err[0, 0]
        print(f"  gamma={gv:7.2f}: FK_mc={fk_mc[i]:.3e} +/-{fk_se[i]:.1e}, "
              f"FF_mc(full)={ffc_mc[i]:.3e}", flush=True)

    D = dict(g=g, o0=o0, ff_mom=ff_mom, fk3=fk3, g_mc=np.array(GAMMA_MC),
             fk_mc=fk_mc, fk_se=fk_se, ffc_mc=ffc_mc, ffc_se=ffc_se,
             sigma=float(sigma), n_lambda=int(n_lambda))
    np.savez(_cache_path(), **D)
    print(f"[cache -> {_cache_path()}]")
    return D


def _plot(D: dict) -> None:
    """Single-panel, house-style figure: the analysis-3 analytic channels
    (O0/FF/FK, sign-aware lines) overlaid with the direct stochastic-Sachs Monte
    Carlo (markers).  No ratio panel; the convergence story lives in the caption."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _plot_style import (PALETTE, annotate_sign_legend, apply_rcparams,  # noqa: E402
                             clip_gamma_axis, plot_signed_line)
    apply_rcparams()
    import matplotlib.pyplot as plt

    g = np.asarray(D["g"], float); o0 = np.asarray(D["o0"], float)
    ff_mom = np.asarray(D["ff_mom"], float); fk3 = np.asarray(D["fk3"], float)
    g_mc = np.asarray(D["g_mc"], float)
    fk_mc = np.asarray(D["fk_mc"], float); fk_se = np.asarray(D["fk_se"], float)
    ffc_mc = np.asarray(D["ffc_mc"], float); ffc_se = np.asarray(D["ffc_se"], float)
    C_O0, C_FF, C_FK = PALETTE[0], PALETTE[1], PALETTE[2]

    fig, ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)

    # analysis-3 analytic channels: sign-aware coloured lines.  O0 changes sign at
    # large gamma (hollow sign markers there; see the filled/hollow note); the sign
    # markers are kept small so the positive FF/FK lines read cleanly.
    _smk = dict(sign_marker_size=3.2, sign_marker_alpha=0.45)
    plot_signed_line(ax, g, o0, color=C_O0, lw=1.7, alpha=0.95, label="O0 (workflow)", **_smk)
    plot_signed_line(ax, g, ff_mom, color=C_FF, lw=1.7, alpha=0.95, label="FF (workflow)", **_smk)
    plot_signed_line(ax, g, fk3, color=C_FK, lw=1.7, alpha=0.95, label="FK (workflow)", **_smk)

    # direct Monte-Carlo overlay (positive in range -> filled markers, matching colour)
    ax.errorbar(g_mc, np.abs(ffc_mc), yerr=ffc_se, fmt="s", color=C_FF, ms=6.5,
                mfc=C_FF, mec=C_FF, capsize=2.5, lw=1.1, zorder=5, label="FF (MC)")
    ax.errorbar(g_mc, np.abs(fk_mc), yerr=fk_se, fmt="o", color=C_FK, ms=6.5,
                mfc=C_FK, mec=C_FK, capsize=2.5, lw=1.1, zorder=5, label="FK (MC)")

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


def main(n_real: int, n_lambda: int, sigma: float, seeds: int,
         ff_real: int, ff_nl: int, from_cache: bool) -> None:
    if from_cache and _cache_path().exists():
        D = dict(np.load(_cache_path(), allow_pickle=True))  # our own cache, safe
        print(f"[from cache {_cache_path()}]")
    else:
        D = _compute(n_real, n_lambda, sigma, seeds, ff_real, ff_nl)
    _plot(D)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    # Converged FK-MC defaults (see the convergence study in the module docstring
    # and FK_NOTES.md): at a MODERATE grid n_lambda=1000 (sigma/dlambda ~ 4-5, the
    # AR(1) field is resolved and the fine-grid heavy-tail bias is absent),
    # decreasing the colored-noise correlation length sigma_lambda drives the FK
    # MC/analytic ratio monotonically toward 1 (the delta / white-noise limit the
    # analytic FK assumes).  sigma_lambda=8 (sigma/dlambda=4.2) is the converged,
    # still-AR(1)-resolved value: FK MC/analytic median ~1.23 at n_real=24000
    # (range [1.08, 1.36] over gamma<300'), down from ~2.4 at the old
    # sigma_lambda=18.  Do NOT
    # increase n_lambda to "converge": at fixed sigma a finer grid INFLATES the
    # ratio (heavy-tailed skewness ~4-6; n_lambda=2000 gives ~15-21x), so n_lambda
    # is held fixed and sigma_lambda is the convergence knob.
    ap.add_argument("--n_real", type=int, default=24_000)
    ap.add_argument("--n_lambda", type=int, default=1000)
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--ff_real", type=int, default=48_000)
    ap.add_argument("--ff_nl", type=int, default=600)
    ap.add_argument("--from-cache", action="store_true",
                    help="skip the MC and re-plot from outputs/appendix_mc_curve.npz "
                         "(for fast figure-style iteration)")
    a = ap.parse_args()
    main(a.n_real, a.n_lambda, a.sigma, a.seeds, a.ff_real, a.ff_nl, a.from_cache)
