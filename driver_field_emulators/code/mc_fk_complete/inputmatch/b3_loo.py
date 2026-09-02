"""B3. Is the marginal block effect driven by one or two blocks?  And is the
per-seed sigma=8 / sigma=4 pair correlated the way CRN implies?"""
import os, numpy as np
from scipy import stats
os.chdir("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
exec(open("inputmatch/b2_ext.py").read().split("rows = []")[0].split("import os")[0]) if False else None
import importlib.util
spec = importlib.util.spec_from_file_location("b2", "inputmatch/b2_ext.py")

def gb(f, sig):
    d = np.load(f, allow_pickle=True)
    m = (d["gamma"] == 1.0) & (d["sigma"] == sig)
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
       ("review/stats/_blockF_s8.npz", "review/stats/_blockF_s4.npz", "F 161803")]
R = []
for f8, f4, lab in SRC:
    if f4 is None:
        s8, e8, ana = gb(f8, 8.0); s4, e4, _ = gb(f8, 4.0)
    else:
        d8 = np.load(f8, allow_pickle=True); d4 = np.load(f4, allow_pickle=True)
        s8, e8, ana = np.asarray(d8["seeds"], float), float(d8["exact"]), float(d8["analytic"])
        s4, e4 = np.asarray(d4["seeds"], float), float(d4["exact"])
    n = min(s8.size, s4.size)
    R.append((lab, s8[:n]/e8, s4[:n]/e4, (2*s4[:n]-s8[:n])/ana))

print("per-seed correlation of the sigma=8 and sigma=4 runs (common random numbers)")
for lab, a8, a4, e in R:
    print(f"  {lab:>10} n={a8.size:>3d}  corr(MC8,MC4) = {np.corrcoef(a8,a4)[0,1]:+.3f}")
allc = np.corrcoef(np.concatenate([r[1]-r[1].mean() for r in R]),
                   np.concatenate([r[2]-r[2].mean() for r in R]))[0,1]
print(f"  pooled (block means removed): {allc:+.3f}")

for name, idx in (("ext", 3), ("sigma8", 1), ("sigma4", 2)):
    grps = [r[idx] for r in R]
    F, p = stats.f_oneway(*grps)
    print(f"\n{name}: full ANOVA F={F:.3f} p={p:.4f}   leave-one-block-out:")
    for j, r in enumerate(R):
        g = [grps[i] for i in range(len(grps)) if i != j]
        Fj, pj = stats.f_oneway(*g)
        print(f"    drop {R[j][0]:>10}: F={Fj:.3f} p={pj:.4f}")
