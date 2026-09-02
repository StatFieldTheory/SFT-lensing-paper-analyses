import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg = float(np.cos(1.0*np.pi/(180*60)))
l = np.linspace(1290, 1330, 41)
v = np.array([ds.sigma2_6x6(cg,float(x))[0,0] for x in l])
for x,y in zip(l,v): print("%8.2f  %.6e"%(x,y))
print("--- coarse global profile ---")
l2 = np.linspace(406,2313,40)
v2 = np.array([ds.sigma2_6x6(cg,float(x))[0,0] for x in l2])
for x,y in zip(l2,v2): print("%8.2f  %.6e"%(x,y))
