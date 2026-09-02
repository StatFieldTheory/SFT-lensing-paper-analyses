"""Is C(t,t) EXACTLY a quadratic Bezier on each table cell, from 3 corner values?

If yes, the density is available in CLOSED FORM and every interpolate-then-
differentiate scheme (R1 included) is an approximation to something exact.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _bootstrap
ds, bgmod, core = _bootstrap.wire()
import corr_op_C_callable as co

LAM = np.load(co.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
print("table lambda nodes:", LAM.size)
print(LAM)

N1 = ds._N1

def Cd(n2, t):
    return np.asarray(co.C_fn(N1, float(t), n2, float(t)), float)

def Cxy(n2, t1, t2):
    return np.asarray(co.C_fn(N1, float(t1), n2, float(t2)), float)

t0 = time.time(); _ = Cd(ds._n2_of_cos(1.0), 1000.0); print("first call %.2fs" % (time.time()-t0))
t0 = time.time()
for _ in range(50): Cd(ds._n2_of_cos(1.0), 1000.0)
print("per scalar C_fn call: %.4f ms" % ((time.time()-t0)/50*1e3))

for gam_arcmin in (1.0, 17.3):
    cg = np.cos(np.deg2rad(gam_arcmin/60.0))
    n2 = ds._n2_of_cos(cg)
    print(f"\n=== gamma = {gam_arcmin}' ===")
    worst = 0.0; worst_where = None
    for i in range(LAM.size-1):
        lo, hi = LAM[i], LAM[i+1]; h = hi-lo
        C00 = Cxy(n2, lo, lo); C01 = Cxy(n2, lo, hi); C10 = Cxy(n2, hi, lo); C11 = Cxy(n2, hi, hi)
        S = C01 + C10
        for u in (0.1, 0.25, 0.5, 0.73, 0.9):
            t = lo + u*h
            pred = (1-u)**2 * C00 + u*(1-u)*S + u**2 * C11
            act = Cd(n2, t)
            den = np.maximum(np.abs(act), 1e-300)
            rel = np.max(np.abs(pred-act)/den)
            if rel > worst:
                worst = rel; worst_where = (lo, hi, u)
    print("max rel |bezier - C_fn| over all 20 cells x 5 u x 9 entries: %.3e at %s" % (worst, worst_where))
