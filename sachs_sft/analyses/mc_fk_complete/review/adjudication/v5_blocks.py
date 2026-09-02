import sys, glob
import numpy as np
base = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/"

def paired(npz, tag):
    d = np.load(npz, allow_pickle=True)
    g, s = d["gamma"], d["sigma"]; seeds = d["seeds"]; ana = d["analytic"]; ex = d["exact"]
    m8 = (np.abs(g-1.0) < 1e-9) & (np.abs(s-8.0) < 1e-9)
    m4 = (np.abs(g-1.0) < 1e-9) & (np.abs(s-4.0) < 1e-9)
    if not (m8.any() and m4.any()): return
    v8 = np.asarray(seeds[np.where(m8)[0][0]], float)
    v4 = np.asarray(seeds[np.where(m4)[0][0]], float)
    a = float(ana[m8][0]); e8 = float(ex[m8][0]); e4 = float(ex[m4][0])
    n = min(v8.size, v4.size); v8, v4 = v8[:n], v4[:n]
    r8, r4 = v8/a, v4/a
    xt = 2*r4 - r8
    def q(x): return x.mean(), x.std(ddof=1)/np.sqrt(x.size)
    print(f"{tag:<28} n={n:3d}  MC/ana(8)={q(r8)[0]:.4f}+-{q(r8)[1]:.4f}  "
          f"MC/ana(4)={q(r4)[0]:.4f}+-{q(r4)[1]:.4f}  "
          f"MC/exact(8)={v8.mean()/e8:.4f}  extrap={q(xt)[0]:.4f}+-{q(xt)[1]:.4f}")
    return xt

print("paired-seed sigma->0 extrapolation of MC/analytic at gamma=1', N=1000")
allx = []
for f, t in [(base+"review/adjudication/_adj_mcblock.npz", "ADJUDICATOR (base 5150271)"),
             (base+"_reseedA.npz", "existing _reseedA"),
             (base+"_reseedB.npz", "existing _reseedB"),
             (base+"_gate3_sigma.npz", "existing _gate3_sigma"),
             (base+"_papergrid.npz", "existing _papergrid"),
             (base+"review/adv_sigma/_mcblockD.npz", "adv_sigma block D")]:
    try:
        x = paired(f, t)
        if x is not None: allx.append(x)
    except Exception as e:
        print(f"{t:<28} [skip: {e}]")
if allx:
    pooled = np.concatenate(allx)
    print(f"\nPOOLED over available blocks: n={pooled.size}  "
          f"{pooled.mean():.4f} +- {pooled.std(ddof=1)/np.sqrt(pooled.size):.4f}")
    means = [x.mean() for x in allx]
    sems = [x.std(ddof=1)/np.sqrt(x.size) for x in allx]
    w = 1/np.array(sems)**2
    mu = (w*np.array(means)).sum()/w.sum()
    chi2 = (w*(np.array(means)-mu)**2).sum()
    print(f"block means: " + ", ".join(f"{m:.4f}+-{s:.4f}" for m, s in zip(means, sems)))
    print(f"inverse-variance mean {mu:.4f} +- {1/np.sqrt(w.sum()):.4f}   "
          f"chi2={chi2:.2f}/{len(means)-1}  scale={np.sqrt(chi2/max(len(means)-1,1)):.2f}")
