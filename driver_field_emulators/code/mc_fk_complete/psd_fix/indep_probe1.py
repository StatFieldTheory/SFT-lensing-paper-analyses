"""Probe 1: table lambda grid, C_fn diagonal structure, timing."""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op

tab = corr_op._TABLE
lam_grid = np.asarray(tab.cl_table.lambda_grid, dtype=float)
print("lambda_grid n =", lam_grid.size)
print(np.array2string(lam_grid, precision=4, max_line_width=200))

N1 = np.array([0.0, 0.0, 1.0])
def n2_of_cos(c):
    c = float(np.clip(c, -1, 1)); s = float(np.sqrt(max(0.0, 1 - c*c)))
    return np.array([s, 0.0, c])

cos1 = float(np.cos(np.deg2rad(1.0/60.0)))
n2 = n2_of_cos(cos1)

t0 = time.time()
for _ in range(20):
    corr_op.C_fn(N1, 1250.0, n2, 1250.0)
print("scalar C_fn time per call: %.4f s" % ((time.time()-t0)/20))

# batch timing
ts = np.full(200, 1250.0)
t0 = time.time()
b = corr_op.C_fn_batch(N1, ts, n2, ts)
print("batch(200) same-t time: %.4f s, shape" % (time.time()-t0), b.shape)
ts = np.linspace(1250.0, 1290.0, 200)
t0 = time.time()
b = corr_op.C_fn_batch(N1, ts, n2, ts)
print("batch(200) varying-t time: %.4f s" % (time.time()-t0))

# structure along the diagonal inside one cell: nodes ... 1300 ...
# fit polynomials of increasing degree to C_00(t,t) inside [1150,1300] and across
def diagC(tarr, cosg):
    nn2 = n2_of_cos(cosg)
    ta = np.asarray(tarr, float)
    return corr_op.C_fn_batch(N1, ta, nn2, ta)

for lo, hi, tag in [(1160.0, 1290.0, "inside cell [1150,1300]"),
                    (1240.0, 1360.0, "straddling node 1300")]:
    t = np.linspace(lo, hi, 41)
    C = diagC(t, cos1)[:, 0, 0]
    print(f"\n{tag}: C_00 range [{C.min():.6e},{C.max():.6e}]")
    x = (t - t.mean())/ (t[-1]-t[0])
    for deg in (2,3,4,5,6,8):
        co = np.polyfit(x, C, deg)
        res = C - np.polyval(co, x)
        print(f"   deg {deg}: max|resid|/max|C| = {np.abs(res).max()/np.abs(C).max():.3e}")
