"""The bigger defect: corr_op interpolates its 21x21 table BILINEARLY, but C has
a ridge along the equal-time diagonal.  At a cell midpoint the bilinear value
averages two on-diagonal corners with two off-diagonal ones and undershoots.

The table stores C_ab(lam_i, lam_i) EXACTLY at all 21 nodes, so a 1-D
interpolation through those values is a clean check of the 2-D interpolant.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
from scipy.interpolate import PchipInterpolator
N1 = ds._N1
g21 = np.load(corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"]

# exact equal-time diagonal AT the nodes (bilinear is exact on grid points)
Cn = np.array([corr_op.C_fn(N1, float(t), N1, float(t))[0, 0] for t in g21])
p = PchipInterpolator(g21, np.log(np.abs(Cn)))
print("C_00 on the equal-time diagonal: bilinear 2-D vs 1-D PCHIP through the")
print("21 exact nodal values, evaluated at cell midpoints")
print(f"{'cell':>18} {'width':>7} {'mid lam':>9} {'bilinear':>13} {'1-D nodal':>13} {'bilin/nodal':>12}")
worst = 0.0
for a, b in zip(g21[:-1], g21[1:]):
    m = 0.5 * (a + b)
    bl = corr_op.C_fn(N1, float(m), N1, float(m))[0, 0]
    nd = float(np.sign(Cn[0]) * np.exp(p(m)))
    r = bl / nd
    worst = min(worst, r) if b - a > 20 else worst
    if b - a > 20:
        print(f"{f'{a:.0f}-{b:.0f}':>18} {b-a:>7.0f} {m:>9.1f} {bl:>13.5e} "
              f"{nd:>13.5e} {r:>12.4f}")
print(f"\nworst mid-cell ratio over the wide cells: {worst:.4f} "
      f"({(worst-1)*100:+.1f}%)")
print()
print("=== why: the four bilinear corners at a mid-cell diagonal point ===")
a, b = 1300.0, 1447.0; m = 0.5 * (a + b)
for (x, y, lbl) in ((a, a, "(lo,lo) on-diagonal"), (b, b, "(hi,hi) on-diagonal"),
                    (a, b, "(lo,hi) OFF-diagonal"), (b, a, "(hi,lo) OFF-diagonal")):
    v = corr_op.C_fn(N1, float(x), N1, float(y))[0, 0]
    print(f"  C({x:.0f},{y:.0f}) {lbl:<22} = {v:.5e}")
print(f"  bilinear average at the midpoint      = "
      f"{corr_op.C_fn(N1, m, N1, m)[0,0]:.5e}")
print(f"  1-D nodal interpolation there         = "
      f"{float(np.sign(Cn[0])*np.exp(p(m))):.5e}")
