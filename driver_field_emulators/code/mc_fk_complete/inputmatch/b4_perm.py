"""B4. (i) Are any two 'independent' blocks secretly the same seeds?
   (ii) Exact permutation p-value for the block effect (no normality assumed).
   (iii) Same test on the SHUFFLED labels as a null calibration."""
import os, numpy as np
from scipy import stats
os.chdir("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
def gb(f, sig):
    d = np.load(f, allow_pickle=True); m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
    if not m.any(): return None, None, None
    i = int(np.where(m)[0][0])
    return np.asarray(d["seeds"][i], float), float(d["exact"][i]), float(d["analytic"][i])
SRC = [("_gate2_highstat.npz", None, "A hs"), ("_gate2_highstat_B.npz", None, "B hsB"),
       ("_papergrid.npz", None, "C paper"), ("review/adv_sigma/_mcblockD.npz", None, "G advsigD"),
       ("review/adjudication/_adj_mcblock.npz", None, "H adj"),
       ("_blk6010001.npz", None, "I blk601"), ("_blk6020002.npz", None, "J blk602"),
       ("_blk6030003.npz", None, "K blk603"), ("_blk6040004.npz", None, "L blk604"),
       ("review/stats/_blockD_s8.npz", "review/stats/_blockD_s4.npz", "D 8675309"),
       ("review/stats/_blockE_s8.npz", "review/stats/_blockE_s4.npz", "E 271828"),
       ("review/stats/_blockF_s8.npz", "review/stats/_blockF_s4.npz", "F 161803"),
       ("_gate3_sigma.npz", None, "* gate3 (base=A)")]
R = []
for f8, f4, lab in SRC:
    if f4 is None:
        s8, e8, ana = gb(f8, 8.0); s4, e4, _ = gb(f8, 4.0)
    else:
        d8 = np.load(f8, allow_pickle=True); d4 = np.load(f4, allow_pickle=True)
        s8, e8, ana = np.asarray(d8["seeds"], float), float(d8["exact"]), float(d8["analytic"])
        s4, e4 = np.asarray(d4["seeds"], float), float(d4["exact"])
    n = min(s8.size, s4.size)
    R.append((lab, s8[:n], s4[:n], (2*s4[:n]-s8[:n])/ana))
print("duplicate-seed audit (shared per-seed VALUES between 'independent' blocks):")
bad = 0
for i in range(len(R)):
    for j in range(i+1, len(R)):
        common = np.intersect1d(R[i][1], R[j][1])
        if common.size:
            bad += 1
            print(f"  !! {R[i][0]} and {R[j][0]} share {common.size} identical "
                  f"sigma=8 values (of {R[i][1].size}/{R[j][1].size})")
print(f"  {'no shared values' if bad==0 else str(bad)+' overlapping pairs'} "
      f"among {len(R)} blocks")
R = [r for r in R if not r[0].startswith("*")]
print()
rng = np.random.default_rng(7)
for name, idx in (("ext", 3), ("sigma8", 1), ("sigma4", 2)):
    grps = [r[idx] for r in R]
    F0 = stats.f_oneway(*grps).statistic
    allv = np.concatenate(grps); n = [g.size for g in grps]
    cnt = 0; NS = 20000
    for _ in range(NS):
        p = rng.permutation(allv); sp = np.split(p, np.cumsum(n)[:-1])
        cnt += stats.f_oneway(*sp).statistic >= F0
    print(f"{name:>7}: F={F0:.3f}  permutation p = {(cnt+1)/(NS+1):.4f}  "
          f"(parametric {stats.f_oneway(*grps).pvalue:.4f})")
