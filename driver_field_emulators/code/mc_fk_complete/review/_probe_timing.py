import sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import _bootstrap
ds, bg, core = _bootstrap.wire()
print("wired ok")
t0=time.time()
B = bg.Background(lam_min=120, lam_max=2313.029+5)
print("bg build", time.time()-t0)
gam = 1.0/60.0*np.pi/180.0
cg = float(np.cos(gam))
print("cos gamma", repr(cg), "1-cos", 1-cg)
t0=time.time()
S = ds.sigma2_6x6(cg, 1000.0)
print("first sigma2_6x6", time.time()-t0)
t0=time.time()
for l in np.linspace(406, 2313, 50):
    S = ds.sigma2_6x6(cg, float(l))
print("50 sigma2_6x6", time.time()-t0)
np.set_printoptions(precision=4, suppress=False, linewidth=180)
print("Sigma2 at 1000:\n", ds.sigma2_6x6(cg,1000.0))
ev = np.linalg.eigvalsh(ds.sigma2_6x6(cg,1000.0))
print("eigs:", ev, "ratio min/max", ev.min()/ev.max())
t0=time.time()
Z = ds.zeta6(cg, 1000.0)
print("first zeta6", time.time()-t0, "shape", Z.shape)
t0=time.time()
for l in np.linspace(406, 2313, 10):
    Z = ds.zeta6(cg, float(l))
print("10 zeta6", time.time()-t0)
