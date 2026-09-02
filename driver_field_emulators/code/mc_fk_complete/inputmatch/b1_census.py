"""B1. Every independent gamma=1' block on disk, analysed with the right test.

The chi2-with-estimated-SEM test used previously is biased high when n is small
and unequal.  One-way ANOVA (and its Welch variant) is the correct test for a
per-block common offset.
"""
from __future__ import annotations
import os, numpy as np
from scipy import stats

os.chdir("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")

def grid_block(f, sig, lab, tag):
    d = np.load(f, allow_pickle=True)
    m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
    if not m.any(): return None
    i = int(np.where(m)[0][0])
    return (lab, tag, np.asarray(d["seeds"][i], float), float(d["exact"][i]))

def flat_block(f, lab, tag):
    d = np.load(f, allow_pickle=True)
    return (lab, tag, np.asarray(d["seeds"], float), float(d["exact"]))

def collect(sig):
    out = []
    for f, lab, tag in [("_gate2_highstat.npz", "A hs 20260827", "early"),
                        ("_gate2_highstat_B.npz", "B hsB 314159001", "early"),
                        ("_papergrid.npz", "C papergrid", "early"),
                        ("review/adv_sigma/_mcblockD.npz", "G advsigD", "mid"),
                        ("review/adjudication/_adj_mcblock.npz", "H adj 5150271", "mid"),
                        ("_blk6010001.npz", "I blk601", "late"),
                        ("_blk6020002.npz", "J blk602", "late"),
                        ("_blk6030003.npz", "K blk603", "late"),
                        ("_blk6040004.npz", "L blk604", "late")]:
        r = grid_block(f, sig, lab, tag)
        if r: out.append(r)
    s = "8" if sig == 8.0 else "4"
    for f, lab, tag in [(f"review/stats/_blockD_s{s}.npz", "D 8675309", "mid"),
                        (f"review/stats/_blockE_s{s}.npz", "E 271828183", "mid"),
                        (f"review/stats/_blockF_s{s}.npz", "F 161803399", "mid")]:
        if os.path.exists(f): out.append(flat_block(f, lab, tag))
    if sig == 8.0 and os.path.exists("review/stats/_bsize2000.npz"):
        out.append(flat_block("review/stats/_bsize2000.npz", "M bsize2000", "mid"))
    return out

for sig in (8.0, 4.0):
    B = collect(sig)
    ex = B[0][3]
    assert all(abs(b[3] - ex) < 1e-18 for b in B)
    print("=" * 92)
    print(f"gamma=1'  sigma_lambda={sig}  n_real=24000 batch=2000 N=1000   "
          f"exact = {ex:.7e}")
    print("=" * 92)
    print(f"{'block':>17} {'when':>6} {'n':>4} {'mean/exact':>11} {'sem':>8} "
          f"{'sd/exact':>9} {'skew':>6} {'exkurt':>7}")
    for lab, tag, s, _ in B:
        print(f"{lab:>17} {tag:>6} {s.size:>4d} {s.mean()/ex:>11.4f} "
              f"{s.std(ddof=1)/np.sqrt(s.size)/ex:>8.4f} {s.std(ddof=1)/ex:>9.4f} "
              f"{stats.skew(s):>6.2f} {stats.kurtosis(s):>7.2f}")
    grps = [b[2] for b in B]
    mu = np.array([g.mean() for g in grps]); n = np.array([g.size for g in grps])
    se = np.array([g.std(ddof=1) / np.sqrt(g.size) for g in grps])
    w = 1 / se ** 2; mm = (w * mu).sum() / w.sum()
    chi2 = float((w * (mu - mm) ** 2).sum()); dof = len(grps) - 1
    F, pF = stats.f_oneway(*grps)
    ag = stats.alexandergovern(*grps); W, pW = ag.statistic, ag.pvalue
    print(f"\n  naive chi2 (estimated SEMs) = {chi2:.2f}/{dof}  p={1-stats.chi2.cdf(chi2,dof):.4f}"
          f"   inflation {np.sqrt(chi2/dof):.2f}")
    print(f"  Alexander-Govern (unequal var) A = {W:.2f} p = {pW:.4f}")
    print(f"  one-way ANOVA  F = {F:.3f}  ({len(grps)}, {n.sum()-len(grps)}) dof  p = {pF:.4f}")
    print(f"  Levene (equal variance) p = {stats.levene(*grps).pvalue:.4f}   "
          f"Bartlett p = {stats.bartlett(*grps).pvalue:.4f}")
    # variance-components: sigma_between^2 estimate
    N_ = n.sum(); k = len(grps)
    MSW = sum(((g - g.mean())**2).sum() for g in grps) / (N_ - k)
    gm = np.concatenate(grps).mean()
    MSB = sum(g.size * (g.mean() - gm) ** 2 for g in grps) / (k - 1)
    n0 = (N_ - (n**2).sum() / N_) / (k - 1)
    s2b = (MSB - MSW) / n0
    print(f"  variance components: within sd = {np.sqrt(MSW)/ex:.4f},  "
          f"between-block sd = {np.sqrt(max(s2b,0))/ex:.4f} "
          f"({'<=' if s2b<=0 else ''}{np.sqrt(abs(s2b))/ex:.4f})")
    allv = np.concatenate(grps)
    print(f"  POOLED n={allv.size}: MC/exact = {allv.mean()/ex:.4f} "
          f"+/- {allv.std(ddof=1)/np.sqrt(allv.size)/ex:.4f}")
    # early vs late
    for t in ("early", "mid", "late"):
        g = np.concatenate([b[2] for b in B if b[1] == t])
        if g.size:
            print(f"    {t:>5}: n={g.size:>3d}  MC/exact = {g.mean()/ex:.4f} "
                  f"+/- {g.std(ddof=1)/np.sqrt(g.size)/ex:.4f}   "
                  f"sd = {g.std(ddof=1)/ex:.4f}")
    print()
