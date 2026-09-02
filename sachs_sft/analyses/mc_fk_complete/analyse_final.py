"""Final numbers, pooling every independent seed block."""
from __future__ import annotations
import numpy as np
from pathlib import Path


def load(p):
    d = np.load(p, allow_pickle=True)
    return {k: d[k] for k in d.files}


def row(d, gam, sig):
    m = (d["gamma"] == gam) & (d["sigma"] == sig)
    i = int(np.where(m)[0][0])
    return np.asarray(d["seeds"][i], float), float(d["exact"][i]), float(d["analytic"][i])


print("=" * 76)
print("gamma = 1', pooled over two independent seed blocks (64 seeds x 24000)")
print("=" * 76)
A, B = load("_gate2_highstat.npz"), load("_gate2_highstat_B.npz")
pooled = {}
for sig in (8.0, 4.0):
    sa, ex, ana = row(A, 1.0, sig)
    sb, _, _ = row(B, 1.0, sig)
    s = np.concatenate([sa, sb]); pooled[sig] = s
    se = s.std(ddof=1) / np.sqrt(s.size)
    print(f"  sigma={sig:>4.1f}  n={s.size}  MC={s.mean():.5e} +/- {se:.2e}"
          f"   MC/exact={s.mean()/ex:.4f} +/- {se/ex:.4f}"
          f"   MC/analytic={s.mean()/ana:.4f} +/- {se/ana:.4f}")
# paired extrapolation, seed by seed (both blocks use the same seeds at both sigma)
_, _, ana = row(A, 1.0, 8.0)
ext = (2.0 * pooled[4.0] - pooled[8.0]) / ana
print(f"\n  linear sigma_lambda -> 0 extrapolation, per seed ({ext.size} seeds):")
print(f"    MC / analytic FK = {ext.mean():.4f} +/- "
      f"{ext.std(ddof=1)/np.sqrt(ext.size):.4f}")
exA = row(A, 1.0, 4.0)[1]; exB = row(A, 1.0, 8.0)[1]
print(f"    exact expectation, same extrapolation = {(2*exA-exB)/ana:.4f}")

print()
print("=" * 76)
print("gamma sweep at sigma_lambda = 8, pooled over the two independent blocks")
print("=" * 76)
RA, RB = load("_reseedA.npz"), load("_reseedB.npz")
G = load("_gate4_gamma.npz")
print(f"  {'gamma':>7} {'n':>4} {'MC/exact':>9} {'+/-':>7} {'sigma_dev':>10} "
      f"{'MC/analytic':>12} {'+/-':>7} {'exact/ana':>10}")
for gam in sorted(set(RA["gamma"])):
    sa, ex, ana = row(RA, gam, 8.0)
    sb, _, _ = row(RB, gam, 8.0)
    s = np.concatenate([sa, sb])
    se = s.std(ddof=1) / np.sqrt(s.size)
    print(f"  {gam:>7.2f} {s.size:>4d} {s.mean()/ex:>9.4f} {se/ex:>7.4f} "
          f"{(s.mean()-ex)/se:>10.1f} {s.mean()/ana:>12.4f} {se/ana:>7.4f} "
          f"{ex/ana:>10.4f}")
print("\n  (the original 8-seed block gave MC/exact = 1.029, 1.042, 1.062, 1.089,")
print("   1.130 at these gamma -- one correlated realisation, not a trend)")

print()
print("=" * 76)
print("small gamma, single 8-seed block (gate 4), for completeness")
print("=" * 76)
print(f"  {'gamma':>7} {'MC/ana s=8':>11} {'MC/ana s=4':>11} {'exact/ana s->0':>15}")
for gam in (0.5, 0.8, 1.3, 2.0):
    s8, e8, ana = row(G, gam, 8.0); s4, e4, _ = row(G, gam, 4.0)
    print(f"  {gam:>7.2f} {s8.mean()/ana:>11.4f} {s4.mean()/ana:>11.4f} "
          f"{(2*e4-e8)/ana:>15.4f}")
