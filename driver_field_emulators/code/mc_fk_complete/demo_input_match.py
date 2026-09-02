"""Are the driving-field statistics actually the SAME on both sides?

Two-point: NO. The MC uses a strict equal-time / white-noise collapse of the
corr_op kernel, which drops its cross-shell tail -- a 12% deficit
(driver_stats' ANCHOR_RATIO = 0.881).
Three-point: YES. Both sides read the same tabulated zeta6, and the deformation
injects its 6-fold symmetrisation exactly (ref_inj/ref_sym = 1.000000000).

So the channels are affected differently: O0 ~ Sigma2, FF ~ Sigma2^2, and
FK ~ zeta with NO Sigma2 dependence at all.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
from pathlib import Path

d = np.load(Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
                 "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz"),
            allow_pickle=True)
g_a, ff_a = np.asarray(d["g"], float), np.asarray(d["ff_mom"], float)

print("FF Monte-Carlo against the analytic FF, with and without the Sigma2 anchor")
print("(anchor multiplies Sigma2 by 1/0.881 = 1.135, so FF should move by 1.135^2 = 1.288)")
print(f"{'gamma':>7} {'anchor':>8} {'MC FF':>13} {'+/-':>10} {'analytic':>13} {'MC/analytic':>12}")
for gam in (1.0, 6.7):
    ana = float(np.interp(gam, g_a, ff_a))
    for anch in (False, True):
        vals = []
        for s in range(4):
            cfg = core.MCConfig(n_real=12_000, batch_size=3_000, n_lambda=2000,
                                use_f_vertex=True, apply_anchor=anch, seed=4400 + 77 * s)
            vals.append(core.simulate_ff_crn(cfg, gam).ff_moment[0, 0])
        v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(v.size)
        print(f"{gam:>7.2f} {str(anch):>8} {v.mean():>13.4e} {sem:>10.1e} "
              f"{ana:>13.4e} {v.mean()/ana:>12.4f}")
