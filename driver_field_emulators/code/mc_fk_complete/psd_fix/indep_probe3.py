"""Probe 3: decisive one-sided derivative test at node 1300 using the SCALAR C_fn."""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

N1 = np.array([0.0, 0.0, 1.0])
def n2_of_cos(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1 - c*c)))
    return np.array([s, 0.0, c])
cosg = float(np.cos(np.deg2rad(1.0/60.0))); n2 = n2_of_cos(cosg)
def c00(t):
    return float(corr_op.C_fn(N1, float(t), n2, float(t))[0, 0])

x0 = 1300.0
print("SCALAR C_fn, entry (0,0), gamma=1'")
print("one-sided 2-point slopes approaching node 1300 from each side:")
for h in (10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1):
    left  = (c00(x0) - c00(x0-h))/h
    right = (c00(x0+h) - c00(x0))/h
    print(f"  h={h:6.2f}  slope- = {left:.8e}   slope+ = {right:.8e}   "
          f"ratio = {right/left:.5f}")
print()
print("value continuity:")
for d in (1e-6, 1e-9):
    print(f"  C(1300-{d:g}) = {c00(x0-d):.12e}   C(1300+{d:g}) = {c00(x0+d):.12e}")
print()
print("second differences inside each cell (should be constant per cell):")
for lo, hi in ((1150.0, 1300.0), (1300.0, 1447.0)):
    h = (hi-lo)/10
    for t in (lo+2*h, lo+5*h, lo+8*h):
        s = (c00(t+h) - 2*c00(t) + c00(t-h))/h**2
        print(f"  cell [{lo:.0f},{hi:.0f}] t={t:7.1f}: C'' = {s:.6e}")
