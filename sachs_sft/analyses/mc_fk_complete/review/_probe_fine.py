import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg = float(np.cos(1.0*np.pi/(180*60)))
for lo,hi in [(1400,1700),(1900,2100)]:
    l=np.arange(lo,hi+1,3.0)
    v=np.array([ds.sigma2_6x6(cg,float(x))[0,0] for x in l])
    print("=== %d-%d ==="%(lo,hi))
    for i in range(0,len(l),2):
        print("  %8.1f  %.5e"%(l[i],v[i]))
