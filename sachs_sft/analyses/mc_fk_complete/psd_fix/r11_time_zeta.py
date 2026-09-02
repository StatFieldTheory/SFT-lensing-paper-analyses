import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
c = float(np.cos(np.deg2rad(1.0/60.0)))
ds.zeta6(c, 1000.0)
t0=time.time()
for l in np.linspace(900,1100,10): ds.zeta6(c, float(l))
print("zeta6 per call (s):", (time.time()-t0)/10)
