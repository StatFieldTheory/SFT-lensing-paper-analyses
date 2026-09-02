"""Is the sampled corr_op equal-time diagonal smooth, or does it carry noise
that the spline derivative amplifies?"""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
from scipy.interpolate import CubicSpline

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
N1 = ds._N1

print("=== A. sample C_ab(t,t) densely and measure its smoothness ===")
lam = np.linspace(1240.0, 1380.0, 281)          # 0.5 Mpc spacing around the defect
C = np.array([corr_op.C_fn(N1, float(t), N1, float(t)) for t in lam])
D4 = np.asarray(bg.D(lam), float) ** 4
f = D4[:, None, None] * C
for (a, b) in ((0, 0), (1, 1), (2, 2), (0, 1)):
    y = f[:, a, b]
    # 4th difference is a clean noise probe: ~0 for a smooth function
    d4 = np.diff(y, 4)
    smooth = np.abs(y).max()
    print(f"  entry ({a},{b}): |f|max={smooth:.4e}  rms 4th-diff/|f|max = "
          f"{d4.std()/smooth:.3e}  max|4th diff|/|f|max = {np.abs(d4).max()/smooth:.3e}")

print()
print("=== B. is C(t,t) itself PSD, and are its INCREMENTS PSD? ===")
eC = np.array([np.linalg.eigvalsh(f[i]) for i in range(len(lam))])
print(f"  min eig of D^4 C over the window: {eC.min():+.4e} (should be >= 0)")
inc = f[1:] - f[:-1]
einc = np.array([np.linalg.eigvalsh(inc[i]) for i in range(len(inc))])
neg = (einc.min(axis=1) < 0).sum()
print(f"  increments over 0.5 Mpc: {neg}/{len(inc)} have a negative eigenvalue")
print(f"  most negative increment eigenvalue / median positive: "
      f"{einc.min()/np.median(einc[einc>0]):+.4f}")
for h in (1, 2, 4, 8, 24):          # coarser increments, h*0.5 Mpc
    inc = f[h:] - f[:-h]
    e = np.array([np.linalg.eigvalsh(inc[i]) for i in range(len(inc))])
    print(f"    increment over {h*0.5:5.1f} Mpc: "
          f"{(e.min(axis=1) < 0).sum():3d}/{len(inc)} negative, "
          f"worst {e.min()/np.median(e[e>0]):+.4f} of median")

print()
print("=== C. spline derivative vs direct finite differences of the dense data ===")
nodes = np.linspace(406.0, 2328.0, 160)
Cn = np.array([corr_op.C_fn(N1, float(t), N1, float(t)) for t in nodes])
fn = np.asarray(bg.D(nodes), float)[:, None, None] ** 4 * Cn
sp = [[CubicSpline(nodes, fn[:, a, b]) for b in range(3)] for a in range(3)]
print(f"{'lam':>8} {'spline deriv minEig':>20} {'FD 0.5Mpc minEig':>18} "
      f"{'FD 4Mpc minEig':>16} {'FD 12Mpc minEig':>16}")
fd = {}
for h in (0.5, 4.0, 12.0):
    fd[h] = None
for t in (1295.6, 1303.2, 1307.0, 1308.9, 1312.7, 1320.0):
    d4t = float(bg.D(t)) ** 4
    S = np.array([[sp[a][b](t, 1) for b in range(3)] for a in range(3)]) / d4t
    S = 0.5 * (S + S.T)
    row = [np.linalg.eigvalsh(S).min()]
    for h in (0.5, 4.0, 12.0):
        tp, tm = t + h / 2, t - h / 2
        fp = float(bg.D(tp)) ** 4 * corr_op.C_fn(N1, tp, N1, tp)
        fm = float(bg.D(tm)) ** 4 * corr_op.C_fn(N1, tm, N1, tm)
        G = (fp - fm) / h / d4t
        G = 0.5 * (G + G.T)
        row.append(np.linalg.eigvalsh(G).min())
    print(f"{t:>8.1f} {row[0]:>20.4e} {row[1]:>18.4e} {row[2]:>16.4e} {row[3]:>16.4e}")
