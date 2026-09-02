"""The anchor integral with a CONVERGED quadrature, for all four Sigma2 sources.

driver_stats.order0_mc uses a 256-node Gauss-Legendre rule.  That rule is exact
for the smooth (corrected-diagonal) integrand but NOT for the production one,
whose Sigma2 is a sawtooth with 21 slope discontinuities.  This recomputes the
same integral on a fine uniform trapezoid grid so the four builders can be
compared without the rule's own error contaminating the difference.
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"
sys.path.insert(0, ROOT)
import _bootstrap  # noqa: E402

ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, ROOT + "/psd_fix")
import sigma2_repaired as R  # noqa: E402
from t5_order0_impact import CorrectedDiag  # noqa: E402

A = 406.01   # NOT exactly lam_lo: R0Reference.matrix has a degenerate
             # cell at x == lam_lo (searchsorted -> lo == hi) and returns an
             # extrapolated spike there; the omitted 0.01 Mpc sliver is ~1e-7
             # of the integral.


def main() -> None:
    bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
    lam_f = bgm.LAM_SOURCE_BASELINE
    n = 4001                     # Sigma2 evaluations per builder per gamma
    la = np.linspace(A, lam_f, n)
    G = ds.order0_window(bg, la, lam_f, n_gauss=256)
    builders = {
        "current": ds.Sigma2Builder(background=bg, apply_c0=False),
        "R1 knots": R.R1Knots(background=bg, apply_c0=False),
        "R0 (prod C)": R.R0Reference(background=bg, apply_c0=False),
        "CORR diag": CorrectedDiag(background=bg, apply_c0=False),
    }
    for gam in (1.0, 17.3):
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        print(f"\n=== gamma = {gam}'  (trapezoid, n = {n}) ===")
        vals = {}
        for k, b in builders.items():
            s = np.array([b.matrix(cosg, float(t))[0, 0] for t in la])
            vals[k] = float(np.trapezoid(s * G ** 2, la))
        base = vals["current"]
        r0 = vals["R0 (prod C)"]
        print(f"{'builder':>13} {'O0 (converged)':>16} {'vs current':>12} {'vs R0':>12}")
        for k, v in vals.items():
            print(f"{k:>13} {v:>16.6e} {v/base-1:>+11.3%} {v/r0-1:>+11.3%}")


if __name__ == "__main__":
    main()
