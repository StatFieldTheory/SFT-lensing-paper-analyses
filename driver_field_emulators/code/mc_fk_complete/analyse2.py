"""Final assembly: paired-seed extrapolations and the honest error bars.

The MC/analytic ratios at different (gamma, sigma) share a seed list, so their
errors are strongly correlated.  Extrapolations are therefore done PER SEED and
the scatter of the extrapolated values is the error bar; the per-point SEMs are
not combined in quadrature.
"""
from __future__ import annotations
import numpy as np
from pathlib import Path


def load(p):
    d = np.load(p, allow_pickle=True)
    return {k: d[k] for k in d.files}


def paired_extrap(seeds_lo, seeds_hi, s_lo, s_hi, denom):
    """Linear-in-sigma extrapolation to sigma = 0, seed by seed."""
    w = s_hi / (s_hi - s_lo)
    ext = (np.asarray(seeds_lo) * w - np.asarray(seeds_hi) * (w - 1.0)) / denom
    return ext.mean(), ext.std(ddof=1) / np.sqrt(ext.size), ext.size


print("=" * 78)
print("A.  sigma_lambda -> 0 : the estimator's EXACT expectation vs analytic FK")
print("=" * 78)
d3 = load("_gate3_sigma.npz")
o = np.argsort(d3["sigma"])
s = d3["sigma"][o]; dev = 1.0 - (d3["exact"] / d3["ref_inj"])[o]
print(f"  {'sigma':>7} {'exact/ref_inj':>14} {'deficit':>10} {'deficit ratio':>14}")
for i in range(s.size):
    rr = dev[i] / dev[i - 1] if i else np.nan
    print(f"  {s[i]:>7.1f} {1-dev[i]:>14.5f} {dev[i]:>10.5f} {rr:>14.3f}")
p = np.polyfit(np.log(s), np.log(dev), 1)
print(f"  deficit ~ sigma^{p[0]:.3f}  (pure linear would be 1.000) -> "
      f"the sigma -> 0 intercept is 1.0000")

print()
print("=" * 78)
print("B.  Monte-Carlo at gamma = 1', high statistics (32 seeds x 24000)")
print("=" * 78)
h = load("_gate2_highstat.npz")
i8 = int(np.where(h["sigma"] == 8.0)[0][0]); i4 = int(np.where(h["sigma"] == 4.0)[0][0])
ana = float(h["analytic"][i8])
for i in (i8, i4):
    sd = np.asarray(h["seeds"][i], float)
    print(f"  sigma={h['sigma'][i]:>4.1f}  MC = {sd.mean():.5e} +/- "
          f"{sd.std(ddof=1)/np.sqrt(sd.size):.2e}   MC/analytic = "
          f"{sd.mean()/ana:.4f} +/- {sd.std(ddof=1)/np.sqrt(sd.size)/ana:.4f}"
          f"   MC/exact = {sd.mean()/h['exact'][i]:.4f}")
m, e, n = paired_extrap(h["seeds"][i4], h["seeds"][i8], 4.0, 8.0, ana)
print(f"  paired-seed extrapolation to sigma = 0 ({n} seeds):"
      f"  MC/analytic = {m:.4f} +/- {e:.4f}")

print()
print("=" * 78)
print("C.  gamma sweep, MC/analytic extrapolated to sigma = 0 (paired seeds)")
print("=" * 78)
g = load("_gate4_gamma.npz")
print(f"  {'gamma':>7} {'MC/ana(s->0)':>13} {'+/-':>8} {'exact/ana(s->0)':>16} "
      f"{'exact/ana s=4':>14}")
rows = []
for gam in sorted(set(g["gamma"])):
    m8 = (g["gamma"] == gam) & (g["sigma"] == 8.0)
    m4 = (g["gamma"] == gam) & (g["sigma"] == 4.0)
    a = float(g["analytic"][m8][0])
    mu, se, _ = paired_extrap(g["seeds"][m4][0], g["seeds"][m8][0], 4.0, 8.0, a)
    e0 = (2 * g["exact"][m4][0] - g["exact"][m8][0]) / a
    rows.append((gam, mu, se, e0))
    print(f"  {gam:>7.2f} {mu:>13.4f} {se:>8.4f} {e0:>16.4f} "
          f"{g['exact'][m4][0]/a:>14.4f}")
r = np.array(rows)
chi2 = np.sum(((r[:, 1] - 1.0) / r[:, 2]) ** 2)
print(f"\n  median MC/analytic at sigma->0 : {np.median(r[:,1]):.4f}")
print(f"  median exact/analytic at sigma->0: {np.median(r[:,3]):.4f}  "
      f"range [{r[:,3].min():.4f}, {r[:,3].max():.4f}]")
print(f"  chi2 of MC/analytic against 1 (correlated seeds, so indicative): "
      f"{chi2:.1f} for {r.shape[0]} points")

for f, lbl in (("_reseedA.npz", "seed block A"), ("_reseedB.npz", "seed block B")):
    if not Path(f).exists():
        continue
    print()
    print(f"  independent {lbl} at sigma = 8, 24 seeds:")
    x = load(f)
    print(f"  {'gamma':>7} {'MC/exact':>9} {'+/-':>8} {'MC/analytic':>12} "
          f"{'exact/ana':>10}")
    for i in np.argsort(x["gamma"]):
        sd = np.asarray(x["seeds"][i], float)
        se = sd.std(ddof=1) / np.sqrt(sd.size)
        print(f"  {x['gamma'][i]:>7.2f} {sd.mean()/x['exact'][i]:>9.4f} "
              f"{se/x['exact'][i]:>8.4f} {sd.mean()/x['analytic'][i]:>12.4f} "
              f"{x['exact'][i]/x['analytic'][i]:>10.4f}")
