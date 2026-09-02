import numpy as np
from scipy import stats
for f, lbl in (("_gate2_highstat.npz", "highstat g=1'"),
               ("_gate4_gamma.npz", "gamma sweep")):
    d = np.load(f, allow_pickle=True)
    print(f"=== {lbl} ({f}) ===")
    for i in range(len(d["gamma"])):
        s = np.asarray(d["seeds"][i], float); ex = d["exact"][i]
        if s.size < 8: continue
        print(f"  g={d['gamma'][i]:5.2f} s={d['sigma'][i]:4.1f} n={s.size:3d} "
              f"mean/exact={s.mean()/ex:7.4f} median/exact={np.median(s)/ex:7.4f} "
              f"skew={stats.skew(s):6.2f} min/ex={s.min()/ex:6.3f} "
              f"max/ex={s.max()/ex:6.3f} sem={s.std(ddof=1)/np.sqrt(s.size)/ex:.4f}")
