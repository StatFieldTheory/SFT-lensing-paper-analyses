import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import cleanroom_exact as CR
from _env import ds

N, s = 400, 8.0
cosg = float(np.cos(np.radians(1.0 / 60.0)))
grid = CR.build_grid(N)
S, Z = CR.driver_tables(cosg, grid["lam"])
rho, V, Cov, rp = CR.ar1_stats(S, grid["dlam"], s)
B = CR.smeared_B(Cov, rp, rho, grid["dlam"])
Qn = CR.solve_Q_all(V, Z / (2 * s) ** 2)
Qs = CR.solve_Q_all(B, Z)

fn = CR.fk_exact(grid, Cov, rho, rp, Qn)["total"]
fs = CR.fk_exact(grid, Cov, rho, rp, Qs)["total"]
alpha = float(np.sum(Qs * Qn) / np.sum(Qn * Qn))       # best scalar fit Qs ~ alpha Qn
dQ = Qs - alpha * Qn
print(f"best scalar alpha (Qs ~ alpha Qn) = {alpha:.6f}")
print(f"||dQ||/||Qs||                     = {np.linalg.norm(dQ)/np.linalg.norm(Qs):.4e}")
fd = CR.fk_exact(grid, Cov, rho, rp, dQ)["total"]
print(f"FK(Qn)={fn:+.6e}  alpha*FK(Qn)={alpha*fn:+.6e}  FK(dQ)={fd:+.6e}  "
      f"sum={alpha*fn+fd:+.6e}  FK(Qs)={fs:+.6e}   (linearity check)")
print(f"=> a {np.linalg.norm(dQ)/np.linalg.norm(Qs):.2e} relative change in Q moves FK by "
      f"{abs(fd/(alpha*fn)):.2f} x  -> amplification ~"
      f"{abs(fd/(alpha*fn))/(np.linalg.norm(dQ)/np.linalg.norm(Qs)):.3g}")

# --- cancellation depth inside the T1/T2 contractions ----------------------
def depth(Q):
    N_ = grid["N"]; dlam = grid["dlam"]; D2 = grid["D2"]; Wd = grid["Wd"]
    A = np.empty_like(Cov); A[0] = D2[0]*Cov[0]
    for k in range(1, N_): A[k] = rho*A[k-1] + D2[k]*Cov[k]
    h = np.zeros(N_)
    for k in range(N_-2, -1, -1): h[k] = rho*(h[k+1] + D2[k+1])
    m = np.arange(N_-1); WdJ = Wd[1:]; pref = dlam/D2[m]
    num = 0.0; den = 0.0
    for l in range(N_):
        lo = m < l; hi = ~lo
        U = np.empty((N_-1, 6, 6))
        if lo.any():
            ml = m[lo]; U[lo] = (pref[lo]*rp[l-ml])[:, None, None]*A[ml]
        if hi.any():
            mh = m[hi]; coef = h[l] - rp[mh-l]*h[mh]
            U[hi] = pref[hi][:, None, None]*(A[l][None] + coef[:, None, None]*Cov[l][None])
        U0 = U[:, :, 0:3]
        MA = np.einsum("j,jam,jbm->ab", WdJ, U0, U0, optimize=True)
        t = Wd[l]*Q[l, 3]*MA
        num += t.sum(); den += np.abs(t).sum()
    return num, den
num, den = depth(Qn)
print(f"T1_A: net={num:+.6e}  sum|terms|={den:.6e}  cancellation depth={den/abs(num):.3e}"
      f"  -> float64 rel. error ~ {1e-16*den/abs(num):.2e}")
