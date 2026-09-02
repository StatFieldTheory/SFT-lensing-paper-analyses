"""Confirm: corr_op C(t,t) has curvature kinks at its 21 table nodes, and the
160-node interpolating cubic spline rings just past them."""
import sys, numpy as np
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as corr_op
from scipy.interpolate import CubicSpline

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
N1 = ds._N1
grid21 = np.load(corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"]

print("=== A. second derivative of C_00(t,t) across a table node (lam = 1300) ===")
t = np.linspace(1270.0, 1330.0, 241)          # 0.25 Mpc
c = np.array([corr_op.C_fn(N1, float(x), N1, float(x))[0, 0] for x in t])
d2 = np.gradient(np.gradient(c, t), t)
print(f"{'lam':>8} {'C_00':>14} {'d2C/dlam2':>14}")
for i in range(0, 241, 20):
    print(f"{t[i]:>8.2f} {c[i]:>14.6e} {d2[i]:>14.3e}")
left = d2[(t > 1280) & (t < 1298)].mean(); right = d2[(t > 1302) & (t < 1320)].mean()
print(f"  mean d2C left of 1300  = {left:.4e}")
print(f"  mean d2C right of 1300 = {right:.4e}   jump factor {right/left:.2f}")

print()
print("=== B. distance from each pathological lambda to the nearest table node ===")
b = ds.Sigma2Builder(background=bg, apply_c0=False)
lam = np.linspace(406.0, 2313.0, 1000)
mn = np.array([np.linalg.eigvalsh(b.matrix(1.0, float(l))).min() for l in lam])
loc = np.convolve(mn, np.ones(21) / 21, mode="same")
bad = np.where(mn / loc < 0.5)[0]
print(f"  {bad.size} of 1000 lambda have min-eig below half the local mean")
runs, cur = [], [bad[0]]
for i in bad[1:]:
    if i == cur[-1] + 1: cur.append(i)
    else: runs.append(cur); cur = [i]
runs.append(cur)
print(f"{'window [Mpc]':>22} {'nearest table node':>20} {'offset':>9}")
for r in runs:
    lo, hi = lam[r[0]], lam[r[-1]]
    j = int(np.argmin(np.abs(grid21 - lo)))
    print(f"{f'{lo:.1f} - {hi:.1f}':>22} {grid21[j]:>20.1f} {lo-grid21[j]:>9.1f}")

print()
print("=== C. the ringing: spline derivative vs truth, cell by cell near 1300 ===")
nodes = np.linspace(406.0, 2328.0, 160)
fn = np.asarray(bg.D(nodes), float)[:, None, None] ** 4 * np.array(
    [corr_op.C_fn(N1, float(x), N1, float(x)) for x in nodes])
sp = [[CubicSpline(nodes, fn[:, a, bb]) for bb in range(3)] for a in range(3)]
print(f"  the 160-node spline cells bracketing lam=1300: "
      f"{nodes[(nodes>1280)&(nodes<1330)].round(1)}")
print(f"{'lam':>8} {'spline d/dlam minEig':>21} {'FD(4 Mpc) minEig':>18} {'error':>10}")
for x in (1292.0, 1298.0, 1301.0, 1305.0, 1307.0, 1309.0, 1313.0, 1317.0, 1325.0):
    d4x = float(bg.D(x)) ** 4
    S = np.array([[sp[a][bb](x, 1) for bb in range(3)] for a in range(3)]) / d4x
    S = 0.5 * (S + S.T)
    fp = float(bg.D(x + 2)) ** 4 * corr_op.C_fn(N1, x + 2, N1, x + 2)
    fm = float(bg.D(x - 2)) ** 4 * corr_op.C_fn(N1, x - 2, N1, x - 2)
    G = 0.5 * ((fp - fm) / 4.0 / d4x + ((fp - fm) / 4.0 / d4x).T)
    a1, a2 = np.linalg.eigvalsh(S).min(), np.linalg.eigvalsh(G).min()
    print(f"{x:>8.1f} {a1:>21.4e} {a2:>18.4e} {a1/a2-1:>+10.2%}")
