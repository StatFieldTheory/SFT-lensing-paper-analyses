"""What ANCHOR_C0 is, what it fixes, and what it breaks -- channel by channel."""
import sys
import numpy as np
sys.path.insert(0, "..")
import _bootstrap
ds, bgm, mc = _bootstrap.wire()
from pathlib import Path
from ff_functional import ff_terms

cache = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
             "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(cache, allow_pickle=True)
g_a3 = np.asarray(d["g"], float)
o0_a3 = np.asarray(d["o0"], float)

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
print("ANCHOR_RATIO =", ds.ANCHOR_RATIO, "  ANCHOR_C0 =", ds.ANCHOR_C0)
print("\n=== A. is 0.881 a converged number? order0_mc vs n_gauss (anchor OFF) ===")
b = ds.Sigma2Builder(background=bg, apply_c0=False)
print(f"{'gamma':>7} " + " ".join(f"{'ng=%d' % n:>12}" for n in (48, 96, 192, 384, 768))
      + f"{'/O0_a3(384)':>13}")
for gam in (0.5, 1.0, 2.6, 6.7, 17.3, 44.4):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    v = [float(ds.order0_mc(cosg, lam_f=bgm.LAM_SOURCE_BASELINE, builder=b, n_gauss=n))
         for n in (48, 96, 192, 384, 768)]
    a3 = float(np.exp(np.interp(np.log(gam), np.log(g_a3), np.log(o0_a3))))
    print(f"{gam:>7.2f} " + " ".join(f"{x:12.5e}" for x in v) + f"{v[3]/a3:13.4f}")

print("\n  -> the default n_gauss=48 used when ANCHOR_RATIO was fixed is NOT converged.")
print("     0.881 vs the converged 0.890 is a quadrature artefact, not physics.")

print("\n=== B. FF term-by-term: which parts of the diagram feel the collapse? ===")
for gam in (1.0, 17.3):
    rf = ff_terms(gam, n=241, kind="full")
    rm = ff_terms(gam, n=241, kind="mc")
    print(f"  gamma = {gam}'")
    for k in ("conn22", "disc", "k1k3", "k3k1", "moment"):
        a, c = rf[k][0, 0], rm[k][0, 0]
        print(f"    {k:8s} full={a:12.5e}  mc={c:12.5e}  mc/full={c/a:8.5f}")

print("\n=== C. what the anchor does to each channel (it rescales Sigma2 by 1/0.881) ===")
c0 = ds.ANCHOR_C0
print(f"  O0 ~ Sigma2^1 : anchor multiplies the MC O0 by {c0:.4f}")
print(f"  FF ~ Sigma2^2 : anchor multiplies the MC FF by {c0**2:.4f}")
print(f"  FK ~ Sigma2^0 : anchor multiplies the MC FK by {1.0:.4f} (exactly)")
