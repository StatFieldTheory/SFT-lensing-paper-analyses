"""Probe 2: is C_ab(t,t) exactly piecewise quadratic in every cell?  Is it C^1?"""
import sys
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

NODES = np.asarray(corr_op._TABLE.cl_table.lambda_grid, float)
N1 = np.array([0.0, 0.0, 1.0])
def n2_of_cos(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1 - c*c)))
    return np.array([s, 0.0, c])

def diag(tarr, cosg):
    ta = np.asarray(tarr, float)
    return corr_op.C_fn_batch(N1, ta, n2_of_cos(cosg), ta)

for gam in (1.0, 17.3):
    cosg = float(np.cos(np.deg2rad(gam/60.0)))
    print(f"\n===== gamma = {gam}'  (cos = {cosg:.12f}) =====")
    worst = 0.0
    for i in range(NODES.size - 1):
        lo, hi = NODES[i], NODES[i+1]
        eps = 0.02*(hi-lo)
        t = np.linspace(lo+eps, hi-eps, 25)
        C = diag(t, cosg)
        x = (t - 0.5*(lo+hi)) / (hi-lo)
        msg = []
        for a in range(3):
            for b in range(3):
                y = C[:, a, b]
                sc = np.abs(y).max()
                if sc < 1e-300:
                    continue
                r2 = np.abs(y - np.polyval(np.polyfit(x, y, 2), x)).max()/sc
                r3 = np.abs(y - np.polyval(np.polyfit(x, y, 3), x)).max()/sc
                worst = max(worst, r2)
                if r2 > 1e-12:
                    msg.append(f"({a}{b}) r2={r2:.2e} r3={r3:.2e}")
        if msg:
            print(f"  cell [{lo:.0f},{hi:.0f}]  NONQUADRATIC: " + "; ".join(msg))
    print(f"  worst deg-2 relative residual over all 20 cells x 9 entries: {worst:.3e}")

# C^1 test at interior nodes: fit quadratic each side, compare value/slope/curv
print("\n===== continuity at interior nodes (entry 00, gamma=1') =====")
cosg = float(np.cos(np.deg2rad(1.0/60.0)))
print(f"{'node':>8} {'jump C/C':>12} {'jump C1/C1':>12} {'C2 left':>12} {'C2 right':>12} {'ratio':>8}")
for i in range(1, NODES.size-1):
    x0 = NODES[i]
    out = []
    for lo, hi in ((NODES[i-1], NODES[i]), (NODES[i], NODES[i+1])):
        eps = 0.02*(hi-lo)
        t = np.linspace(lo+eps, hi-eps, 15)
        y = diag(t, cosg)[:, 0, 0]
        co = np.polyfit(t - x0, y, 2)          # y = co0 (t-x0)^2 + co1 (t-x0) + co2
        out.append((co[2], co[1], 2*co[0]))
    (v0,d0,s0),(v1,d1,s1) = out
    print(f"{x0:>8.0f} {abs(v1-v0)/abs(v0):>12.2e} {abs(d1-d0)/abs(d0):>12.2e} "
          f"{s0:>12.4e} {s1:>12.4e} {s1/s0:>8.3f}")
