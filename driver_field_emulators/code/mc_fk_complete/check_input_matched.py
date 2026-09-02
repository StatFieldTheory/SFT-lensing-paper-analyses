"""Are the two sides actually being given the SAME driving-field statistics?

The comparison the paper wants is: hold the driving-field statistics fixed, and
ask whether the Monte-Carlo and the semi-analytic fold agree.  That requires the
MC's Sigma2 and zeta to BE the fold's Sigma2 and zeta.

zeta: shared by construction (same table, same callable).
Sigma2: the MC builds its own from corr_op by lambda-differentiation, and runs
with apply_anchor=False, i.e. WITHOUT the ANCHOR_C0 = 1/0.881 that exists
precisely to reconcile the two.  So this needs checking.
"""
import numpy as np
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from pathlib import Path

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3 = np.asarray(d["g"], float)
o0_a3, ff_a3 = np.asarray(d["o0"], float), np.asarray(d["ff_mom"], float)

print("=== 1. Order-0: does the MC's driving variance match the formalism's? ===")
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
print(f"{'gamma':>8} {'anchor OFF':>12} {'anchor ON':>12}   (order0_mc / O0_analysis3)")
for gam in (1.0, 2.6, 6.7, 17.3, 44.4):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    a3 = float(np.interp(gam, g_a3, o0_a3))
    r = []
    for anc in (False, True):
        b = ds.Sigma2Builder(background=bg, apply_c0=anc)
        r.append(float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE,
                                    builder=b, n_gauss=256)) / a3)
    print(f"{gam:>8.2f} {r[0]:>12.4f} {r[1]:>12.4f}")

print()
print("=== 2. FF: the markers, anchor OFF (as published) vs anchor ON ===")
print("   FF ~ Sigma2^2, so a mismatched variance shows up squared.")
print(f"{'gamma':>8} {'FF anchor OFF':>15} {'FF anchor ON':>14} {'analytic FF':>13} "
      f"{'OFF/ana':>9} {'ON/ana':>8}")
for gam in (1.0, 6.7):
    a3 = float(np.interp(gam, g_a3, ff_a3))
    vals = []
    for anc in (False, True):
        acc = []
        for s in range(4):
            cfg = mc.MCConfig(n_real=12000, batch_size=3000, n_lambda=2000,
                              use_f_vertex=True, apply_anchor=anc, seed=11 + 7 * s)
            acc.append(mc.simulate_ff_crn(cfg, gam).ff_moment[0, 0])
        vals.append(float(np.mean(acc)))
    print(f"{gam:>8.2f} {vals[0]:>15.4e} {vals[1]:>14.4e} {a3:>13.4e} "
          f"{vals[0]/a3:>9.4f} {vals[1]/a3:>8.4f}")
print()
print("  expected ratio ON/OFF from FF ~ Sigma2^2 :", f"{(1/0.881)**2:.4f}")
