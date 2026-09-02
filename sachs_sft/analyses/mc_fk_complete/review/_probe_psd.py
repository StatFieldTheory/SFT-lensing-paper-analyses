import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import _bootstrap
ds, bg, _ = _bootstrap.wire()
cg = float(np.cos(1.0*np.pi/(180*60)))
lam = 406.0 + (2313.029-406.0)/999*np.arange(1000)
bad=[]
mins=[]
diag0=[]
for k,l in enumerate(lam):
    S = ds.sigma2_6x6(cg, float(l))
    d = np.linalg.eigvalsh(S)
    mins.append(d.min()/d.max())
    diag0.append(S[0,0])
    if d.min() <= 0: bad.append(k)
mins=np.array(mins); diag0=np.array(diag0)
print("n nodes with a non-positive eig:", len(bad))
print("first/last bad k:", bad[:5], bad[-5:] if bad else None)
print("lam of bad:", lam[bad][:5], lam[bad][-5:] if bad else None)
print("min eig ratio over grid: %.3e  max: %.3e"%(mins.min(), mins.max()))
print("diag00 min %.3e"%diag0.min(), " any negative diag00:", (diag0<0).sum())
# where are the negatives concentrated?
import collections
print("fraction of grid bad: %.3f"%(len(bad)/1000))
print("mins percentiles:", np.percentile(mins,[0,1,5,50,95,100]))
# magnitude of the negative eig relative
S=ds.sigma2_6x6(cg,float(lam[bad[0]])) if bad else None
if S is not None:
    print("example bad matrix eigs:", np.linalg.eigvalsh(S))
