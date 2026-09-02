"""Effect of each repair on the validated Order-0 anchor.

driver_stats.ANCHOR_RATIO = 0.881 is the median of
    order0_mc(gamma; apply_c0=False) / O0_analysis3(gamma)
over the paper's gamma grid.  If a repair moves that ratio the constant would
have to be re-derived; if it does not, the repair is anchor-neutral.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/psd_fix")
import sigma2_repaired as R
from pathlib import Path

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3, o0_a3 = np.asarray(d["g"], float), np.asarray(d["o0"], float)

GAM = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 17.3, 44.4, 114.3, 293.9, 596.9])
CLS = [("current", ds.Sigma2Builder), ("R1 knots", R.R1Knots),
       ("R2 density", R.R2Density), ("R0 reference", R.R0Reference)]

for anchored in (False, True):
    tag = "anchored (apply_c0=True)" if anchored else "raw (apply_c0=False)"
    print(f"=== order0_mc / O0_analysis3, {tag} ===")
    builders = {k: cls(background=bg, apply_c0=anchored) for k, cls in CLS}
    print(f"{'gamma':>8} {'O0 analysis-3':>14} " + " ".join(f"{k:>13}" for k in builders))
    rows = {k: [] for k in builders}
    for gam in GAM:
        cosg = float(np.cos(np.deg2rad(gam / 60.0)))
        a3 = float(np.interp(gam, g_a3, o0_a3))
        vals = []
        for k, b in builders.items():
            v = float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                   builder=b, n_gauss=256))
            rows[k].append(v / a3); vals.append(v / a3)
        print(f"{gam:>8.2f} {a3:>14.5e} " + " ".join(f"{v:>13.5f}" for v in vals))
    print(f"{'':>23} " + " ".join(f"{'-'*13}" for _ in builders))
    med = {k: float(np.median(rows[k])) for k in rows}
    print(f"{'median':>23} " + " ".join(f"{med[k]:>13.5f}" for k in builders))
    base = med["current"]
    print(f"{'shift vs current':>23} " +
          " ".join(f"{med[k]/base-1:>+12.4%} " for k in builders))
    print()
