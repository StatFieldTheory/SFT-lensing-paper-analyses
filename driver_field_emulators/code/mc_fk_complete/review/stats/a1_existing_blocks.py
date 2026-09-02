"""Adversarial audit A1: are the three existing seed blocks mutually consistent,
and is the per-seed sampling distribution normal enough for a seed SEM at n=24-32?

Reads only the author's saved npz products.  No simulation.
"""
import numpy as np
from scipy import stats

BLOCKS = [("../../_gate2_highstat.npz", "A base 20260827"),
          ("../../_gate2_highstat_B.npz", "B base 314159001"),
          ("../../_papergrid.npz", "C base 20260828")]

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

data = {8.0: [], 4.0: []}
meta = {}
for f, lbl in BLOCKS:
    d = np.load(f, allow_pickle=True)
    for sig in (8.0, 4.0):
        m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
        if not m.any():
            continue
        i = int(np.where(m)[0][0])
        s = np.asarray(d["seeds"][i], float)
        data[sig].append((lbl, s))
        meta[sig] = (float(d["exact"][i]), float(d["analytic"][i]))

for sig in (8.0, 4.0):
    exact, ana = meta[sig]
    print("=" * 84)
    print(f"sigma_lambda = {sig}   gamma = 1'   exact = {exact:.6e}   analytic = {ana:.6e}")
    print("=" * 84)
    print(f"{'block':>18} {'n':>4} {'mean':>13} {'sd':>11} {'sem':>11} "
          f"{'skew':>7} {'exkurt':>8} {'SW p':>7} {'mean/exact':>11}")
    means, sems, ns = [], [], []
    allv = []
    for lbl, s in data[sig]:
        se = s.std(ddof=1) / np.sqrt(s.size)
        sw = stats.shapiro(s).pvalue
        print(f"{lbl:>18} {s.size:>4d} {s.mean():>13.6e} {s.std(ddof=1):>11.4e} "
              f"{se:>11.4e} {stats.skew(s):>7.2f} "
              f"{stats.kurtosis(s):>8.2f} {sw:>7.3f} {s.mean()/exact:>11.4f}")
        means.append(s.mean()); sems.append(se); ns.append(s.size)
        allv.append(s - s.mean())
    means = np.array(means); sems = np.array(sems); ns = np.array(ns)

    # --- inter-block consistency ------------------------------------------
    w = 1.0 / sems**2
    mu = float((w * means).sum() / w.sum())
    chi2 = float((w * (means - mu) ** 2).sum())
    dof = len(means) - 1
    p = 1.0 - stats.chi2.cdf(chi2, dof)
    print(f"\n  inverse-variance mean over blocks : {mu:.6e}")
    print(f"  chi2 between blocks               : {chi2:.2f} for {dof} dof   p = {p:.4f}")
    print(f"  -> scale factor on the error bar  : sqrt(chi2/dof) = {np.sqrt(chi2/dof):.2f}")
    # pairwise
    for i in range(len(means)):
        for j in range(i + 1, len(means)):
            dd = means[i] - means[j]
            sd = np.hypot(sems[i], sems[j])
            print(f"     {data[sig][i][0]:>18} vs {data[sig][j][0]:<18} "
                  f"diff = {dd:+.4e}  ({dd/sd:+.2f} sigma)")

    # --- pooled ------------------------------------------------------------
    pooled = np.concatenate([s for _, s in data[sig]])
    se_pool = pooled.std(ddof=1) / np.sqrt(pooled.size)
    print(f"\n  POOLED n={pooled.size}  mean={pooled.mean():.6e} "
          f"+/- {se_pool:.3e}   ratio/exact = {pooled.mean()/exact:.4f} "
          f"+/- {se_pool/exact:.4f}   ratio/analytic = {pooled.mean()/ana:.4f} "
          f"+/- {se_pool/ana:.4f}")

    # --- normality of the residual (block means removed) --------------------
    res = np.concatenate(allv)
    print(f"\n  pooled within-block residuals (n={res.size}):")
    print(f"    skew = {stats.skew(res):+.3f}   excess kurtosis = {stats.kurtosis(res):+.3f}")
    print(f"    Shapiro-Wilk  p = {stats.shapiro(res).pvalue:.4f}")
    print(f"    D'Agostino K2 p = {stats.normaltest(res).pvalue:.4f}")
    ad = stats.anderson(res / res.std(ddof=1), dist="norm")
    print(f"    Anderson-Darling A2 = {ad.statistic:.3f} "
          f"(5% crit {ad.critical_values[2]:.3f})")
    print(f"    Jarque-Bera   p = {stats.jarque_bera(res).pvalue:.4f}")
    q = np.percentile(res / res.std(ddof=1), [1, 5, 25, 50, 75, 95, 99])
    print(f"    standardised quantiles 1/5/25/50/75/95/99: "
          + " ".join(f"{x:+.2f}" for x in q))
    print(f"    max |z| = {np.abs(res/res.std(ddof=1)).max():.2f} "
          f"(Gaussian expectation for n={res.size}: ~{stats.norm.ppf(1-1/(2*res.size)):.2f})")

    # --- bootstrap of the pooled mean vs the SEM ---------------------------
    rng = np.random.default_rng(7)
    bs = np.array([rng.choice(pooled, pooled.size, replace=True).mean()
                   for _ in range(40000)])
    print(f"\n  bootstrap over seeds (40k):  sd = {bs.std(ddof=1):.3e} "
          f"vs seed SEM {se_pool:.3e}   ratio = {bs.std(ddof=1)/se_pool:.3f}")
    lo, hi = np.percentile(bs, [2.5, 97.5])
    print(f"    percentile 95% CI on MC/exact: "
          f"[{lo/exact:.4f}, {hi/exact:.4f}]   "
          f"(normal-theory: [{(pooled.mean()-1.96*se_pool)/exact:.4f}, "
          f"{(pooled.mean()+1.96*se_pool)/exact:.4f}])")
    # jackknife
    n = pooled.size
    jk = np.array([np.delete(pooled, i).mean() for i in range(n)])
    jse = np.sqrt((n - 1) / n * ((jk - jk.mean()) ** 2).sum())
    print(f"    jackknife SE = {jse:.3e}   ratio to seed SEM = {jse/se_pool:.3f}")
    print()
