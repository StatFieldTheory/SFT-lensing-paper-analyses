"""B2. The statistic the heterogeneity claim is actually about:
ext = 2*MC(sigma=4) - MC(sigma=8), paired per seed, per block.
Twelve independent blocks on disk instead of four.
"""
from __future__ import annotations
import os, numpy as np
from scipy import stats
os.chdir("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")

def gb(f, sig):
    d = np.load(f, allow_pickle=True)
    m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
    if not m.any(): return None, None, None
    i = int(np.where(m)[0][0])
    return np.asarray(d["seeds"][i], float), float(d["exact"][i]), float(d["analytic"][i])

SRC = [("_gate2_highstat.npz", None, "A hs", "early"),
       ("_gate2_highstat_B.npz", None, "B hsB", "early"),
       ("_papergrid.npz", None, "C paper", "early"),
       ("review/adv_sigma/_mcblockD.npz", None, "G advsigD", "mid"),
       ("review/adjudication/_adj_mcblock.npz", None, "H adj", "mid"),
       ("_blk6010001.npz", None, "I blk601", "late"),
       ("_blk6020002.npz", None, "J blk602", "late"),
       ("_blk6030003.npz", None, "K blk603", "late"),
       ("_blk6040004.npz", None, "L blk604", "late"),
       ("review/stats/_blockD_s8.npz", "review/stats/_blockD_s4.npz", "D 8675309", "mid"),
       ("review/stats/_blockE_s8.npz", "review/stats/_blockE_s4.npz", "E 271828", "mid"),
       ("review/stats/_blockF_s8.npz", "review/stats/_blockF_s4.npz", "F 161803", "mid")]

rows = []
for f8, f4, lab, tag in SRC:
    if f4 is None:
        s8, e8, ana = gb(f8, 8.0); s4, e4, _ = gb(f8, 4.0)
    else:
        d8 = np.load(f8, allow_pickle=True); d4 = np.load(f4, allow_pickle=True)
        s8, e8, ana = np.asarray(d8["seeds"], float), float(d8["exact"]), float(d8["analytic"])
        s4, e4 = np.asarray(d4["seeds"], float), float(d4["exact"])
    if s8 is None or s4 is None: continue
    n = min(s8.size, s4.size)
    rows.append((lab, tag, (2 * s4[:n] - s8[:n]) / ana, ana, (2*e4 - e8)/ana))

ana = rows[0][3]; ex0 = rows[0][4]
print(f"paired sigma->0 extrapolation, gamma=1'.  analytic fold = {ana:.6e}")
print(f"exact expectation, same extrapolation / fold = {ex0:.4f}")
print(f"{'block':>12} {'when':>6} {'n':>4} {'ext/fold':>9} {'sem':>7} {'sd':>7} {'exkurt':>7}")
for lab, tag, e, _, _ in rows:
    print(f"{lab:>12} {tag:>6} {e.size:>4d} {e.mean():>9.4f} "
          f"{e.std(ddof=1)/np.sqrt(e.size):>7.4f} {e.std(ddof=1):>7.4f} "
          f"{stats.kurtosis(e):>7.2f}")
grps = [r[2] for r in rows]
mu = np.array([g.mean() for g in grps]); se = np.array([g.std(ddof=1)/np.sqrt(g.size) for g in grps])
w = 1/se**2; mm = (w*mu).sum()/w.sum(); chi2 = float((w*(mu-mm)**2).sum()); dof = len(grps)-1
F, pF = stats.f_oneway(*grps)
ag = stats.alexandergovern(*grps)
allv = np.concatenate(grps)
print(f"\n  naive chi2 = {chi2:.2f}/{dof}  p={1-stats.chi2.cdf(chi2,dof):.4f}  "
      f"inflation {np.sqrt(chi2/dof):.2f}")
print(f"  one-way ANOVA F = {F:.3f} p = {pF:.4f};  Alexander-Govern p = {ag.pvalue:.4f}")
print(f"  Levene p = {stats.levene(*grps).pvalue:.4f}  Bartlett p = {stats.bartlett(*grps).pvalue:.4f}")
n = np.array([g.size for g in grps]); N_=n.sum(); k=len(grps)
MSW = sum(((g-g.mean())**2).sum() for g in grps)/(N_-k)
MSB = sum(g.size*(g.mean()-allv.mean())**2 for g in grps)/(k-1)
n0 = (N_-(n**2).sum()/N_)/(k-1)
print(f"  variance components: within sd {np.sqrt(MSW):.4f}, between-block sd^2 "
      f"= {(MSB-MSW)/n0:+.2e}  -> sd {np.sqrt(abs((MSB-MSW)/n0)):.4f} "
      f"({'REAL' if MSB>MSW else 'consistent with ZERO'})")
print(f"  POOLED n={allv.size}: ext/fold = {allv.mean():.4f} +/- "
      f"{allv.std(ddof=1)/np.sqrt(allv.size):.4f}")
print(f"  MC/exact of the extrapolant = {allv.mean()/ex0:.4f} +/- "
      f"{allv.std(ddof=1)/np.sqrt(allv.size)/ex0:.4f}")
for t in ("early","mid","late"):
    g = np.concatenate([r[2] for r in rows if r[1]==t])
    print(f"    {t:>5}: n={g.size:>3d} ext/fold = {g.mean():.4f} +/- "
          f"{g.std(ddof=1)/np.sqrt(g.size):.4f}  sd={g.std(ddof=1):.4f}")
# time-ordered variance trend test
lab3 = [r[1] for r in rows]
g3 = [np.concatenate([r[2]-r[2].mean() for r in rows if r[1]==t]) for t in ("early","mid","late")]
print(f"  Levene on early/mid/late residuals: p = {stats.levene(*g3).pvalue:.4f}")
