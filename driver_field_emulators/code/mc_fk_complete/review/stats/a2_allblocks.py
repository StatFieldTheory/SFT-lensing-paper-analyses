"""Adversarial audit A2: pool every independent 32/24-seed block at gamma = 1'
(3 of the author's + N of mine) and ask whether the seed SEM is honest.
"""
from __future__ import annotations
import glob, os, sys
import numpy as np
from scipy import stats

os.chdir(os.path.dirname(os.path.abspath(__file__)))

AUTH = [("../../_gate2_highstat.npz", "A 20260827"),
        ("../../_gate2_highstat_B.npz", "B 314159001"),
        ("../../_papergrid.npz", "C 20260828")]
MINE = [("_blockD_s8.npz", "_blockD_s4.npz", "D 8675309"),
        ("_blockE_s8.npz", "_blockE_s4.npz", "E 271828183"),
        ("_blockF_s8.npz", "_blockF_s4.npz", "F 161803399"),
        ("_blockG_s8.npz", "_blockG_s4.npz", "G 424242424"),
        ("_blockH_s8.npz", "_blockH_s4.npz", "H 999331199")]

def collect(sig):
    out = []
    for f, l in AUTH:
        d = np.load(f, allow_pickle=True)
        m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
        if not m.any():
            continue
        i = int(np.where(m)[0][0])
        out.append((l, np.asarray(d["seeds"][i], float), None,
                    float(d["exact"][i]), float(d["analytic"][i])))
    for f8, f4, l in MINE:
        f = f8 if sig == 8.0 else f4
        if not os.path.exists(f):
            continue
        d = np.load(f, allow_pickle=True)
        out.append((l, np.asarray(d["seeds"], float),
                    np.asarray(d["batch_se"], float),
                    float(d["exact"]), float(d["analytic"])))
    return out

report = {}
for sig in (8.0, 4.0):
    B = collect(sig)
    ex, ana = B[0][3], B[0][4]
    assert all(abs(b[3] - ex) < 1e-16 for b in B), "exact expectation differs between blocks!"
    print("=" * 100)
    print(f"gamma = 1', sigma_lambda = {sig}, n_real = 24000, batch = 2000, N = 1000")
    print(f"exact expectation = {ex:.7e}   paper analytic fold = {ana:.7e}")
    print("=" * 100)
    print(f"{'block':>14} {'n':>4} {'MC/analytic':>12} {'seedSEM':>9} {'MC/exact':>9} "
          f"{'sd':>11} {'batchSE':>11} {'sd/bSE':>7} {'skew':>6} {'exk':>6} {'SW p':>6}")
    for l, s, b, _, _ in B:
        se = s.std(ddof=1) / np.sqrt(s.size)
        r = "" if b is None else f"{s.std(ddof=1)/ (b.mean()/0.9776):>7.3f}"
        bs = "" if b is None else f"{b.mean():>11.4e}"
        print(f"{l:>14} {s.size:>4d} {s.mean()/ana:>12.4f} {se/ana:>9.4f} "
              f"{s.mean()/ex:>9.4f} {s.std(ddof=1):>11.4e} {bs:>11} {r:>7} "
              f"{stats.skew(s):>6.2f} {stats.kurtosis(s):>6.2f} "
              f"{stats.shapiro(s).pvalue:>6.3f}")
    mu = np.array([b[1].mean() for b in B])
    se = np.array([b[1].std(ddof=1) / np.sqrt(b[1].size) for b in B])
    w = 1 / se**2
    mm = float((w * mu).sum() / w.sum())
    chi2 = float((w * (mu - mm) ** 2).sum()); dof = len(B) - 1
    pch = 1 - stats.chi2.cdf(chi2, dof)
    F, pan = stats.f_oneway(*[b[1] for b in B])
    print(f"\n  BETWEEN-BLOCK: chi2 = {chi2:.2f}/{dof} dof  p = {pch:.4f}   "
          f"inflation sqrt(chi2/dof) = {np.sqrt(chi2/dof):.2f}")
    print(f"  one-way ANOVA  F = {F:.3f}  p = {pan:.4f}   "
          f"Bartlett p = {stats.bartlett(*[b[1] for b in B]).pvalue:.4f}  "
          f"Levene p = {stats.levene(*[b[1] for b in B]).pvalue:.4f}")
    allv = np.concatenate([b[1] for b in B])
    sep = allv.std(ddof=1) / np.sqrt(allv.size)
    print(f"  POOLED n = {allv.size}: MC = {allv.mean():.6e} +/- {sep:.3e}")
    print(f"     MC/exact    = {allv.mean()/ex:.4f} +/- {sep/ex:.4f}")
    print(f"     MC/analytic = {allv.mean()/ana:.4f} +/- {sep/ana:.4f}")
    bm = np.array([b[1].mean() for b in B])
    print(f"  block-mean scatter estimate of the same error: "
          f"{bm.std(ddof=1)/np.sqrt(bm.size)/ana:.4f}  "
          f"(ratio to pooled SEM {bm.std(ddof=1)/np.sqrt(bm.size)/sep*ana/ana:.2f})")
    res = np.concatenate([b[1] - b[1].mean() for b in B])
    print(f"  pooled residual (n={res.size}): skew {stats.skew(res):+.3f}  "
          f"exkurt {stats.kurtosis(res):+.3f}  SW p {stats.shapiro(res).pvalue:.4f}  "
          f"JB p {stats.jarque_bera(res).pvalue:.4f}")
    rng = np.random.default_rng(11)
    bsd = np.array([rng.choice(allv, allv.size, replace=True).mean() for _ in range(30000)]).std(ddof=1)
    print(f"  bootstrap sd / seed SEM = {bsd/sep:.4f}")
    report[sig] = (B, allv, ex, ana)

# ---- paired sigma -> 0 extrapolation, block by block ----------------------
print()
print("=" * 100)
print("PAIRED sigma_lambda -> 0 EXTRAPOLATION  ext = (2*MC(4) - MC(8)) / analytic")
print("=" * 100)
B8 = {b[0]: b[1] for b in report[8.0][0]}
B4 = {b[0]: b[1] for b in report[4.0][0]}
ana = report[8.0][3]
labs = [l for l in B8 if l in B4]
exts = {}
print(f"{'block':>14} {'n':>4} {'ext mean':>10} {'sem':>8} {'sd':>8}")
for l in labs:
    n = min(B8[l].size, B4[l].size)
    e = (2 * B4[l][:n] - B8[l][:n]) / ana
    exts[l] = e
    print(f"{l:>14} {n:>4d} {e.mean():>10.4f} {e.std(ddof=1)/np.sqrt(n):>8.4f} {e.std(ddof=1):>8.4f}")
mu = np.array([exts[l].mean() for l in labs])
se = np.array([exts[l].std(ddof=1) / np.sqrt(exts[l].size) for l in labs])
w = 1 / se**2; mm = (w * mu).sum() / w.sum()
chi2 = float((w * (mu - mm) ** 2).sum()); dof = len(labs) - 1
alle = np.concatenate([exts[l] for l in labs])
sep = alle.std(ddof=1) / np.sqrt(alle.size)
print(f"\n  between-block chi2 = {chi2:.2f}/{dof}  p = {1-stats.chi2.cdf(chi2,dof):.4f}  "
      f"inflation {np.sqrt(chi2/dof):.2f}")
print(f"  POOLED n={alle.size}:  MC/analytic (sigma->0) = {alle.mean():.4f} +/- {sep:.4f}")
print(f"  inflated by sqrt(chi2/dof):                     {alle.mean():.4f} "
      f"+/- {sep*max(1.0,np.sqrt(chi2/dof)):.4f}")
print(f"  AUTHOR'S QUOTED VALUE (88 seeds, blocks A+B+C): 1.0022 +/- 0.0056")
abc = np.concatenate([exts[l] for l in labs if l[0] in "ABC"])
print(f"  reproduced here from A+B+C only (n={abc.size}): {abc.mean():.4f} +/- "
      f"{abc.std(ddof=1)/np.sqrt(abc.size):.4f}")
mine = np.concatenate([exts[l] for l in labs if l[0] not in "ABC"])
if mine.size:
    print(f"  MY independent blocks only (n={mine.size}):    {mine.mean():.4f} +/- "
          f"{mine.std(ddof=1)/np.sqrt(mine.size):.4f}")
    d = abc.mean() - mine.mean()
    sd = np.hypot(abc.std(ddof=1)/np.sqrt(abc.size), mine.std(ddof=1)/np.sqrt(mine.size))
    print(f"  author's blocks vs mine: difference {d:+.4f} +/- {sd:.4f}  ({d/sd:+.2f} sigma)")
