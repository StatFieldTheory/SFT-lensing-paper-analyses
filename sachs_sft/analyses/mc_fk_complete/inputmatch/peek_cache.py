import numpy as np
from pathlib import Path
p = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses"
         "/sachs_sft/analyses/mc_sachs_2pt/outputs/appendix_mc_curve.npz")
d = np.load(p, allow_pickle=True)
for k in d.files:
    a = np.asarray(d[k])
    print(k, a.shape, a.dtype)
g = np.asarray(d["g"], float)
print("g:", g)
for k in ("o0","ff_mom"):
    print(k, np.asarray(d[k], float))
