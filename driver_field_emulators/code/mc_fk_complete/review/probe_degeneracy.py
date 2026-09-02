import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import cleanroom_exact as CR
from _env import ds

N = 1000
s = 8.0
cosg = float(np.cos(np.radians(1.0 / 60.0)))
grid = CR.build_grid(N)
S, Z = CR.driver_tables(cosg, grid["lam"])
rho, V, Cov, rp = CR.ar1_stats(S, grid["dlam"], s)
B = CR.smeared_B(Cov, rp, rho, grid["dlam"])
k = N // 2
r = float(np.median(B[:, 0, 0] / S[:, 0, 0]))
print("median B/S2 (00) =", r)
dS = np.linalg.eigvalsh(S[k]); dB = np.linalg.eigvalsh(B[k])
print("Sigma2 eig :", np.array2string(dS, precision=6))
print("B      eig :", np.array2string(dB, precision=6))
print("B/S2 eig ratios:", np.array2string(dB/dS, precision=6))
print("elementwise (B/S2) at k:", np.array2string(B[k]/S[k], precision=5))

Qn = ds.solve_Q(S[k], Z[k])                      # == solve_Q(V, Z/(2s)^2)
Qs = ds.solve_Q(B[k], Z[k])
print("max|Qn| =", np.abs(Qn).max(), " max|Qs| =", np.abs(Qs).max(),
      " ratio =", np.abs(Qs).max()/np.abs(Qn).max())
print("Qn[3] =\n", np.array2string(Qn[3], precision=3))
print("Qs[3] =\n", np.array2string(Qs[3], precision=3))

# how much of the FK comes from the small-eigenvalue directions? vary rcond
for rc in (1e-8, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
    Qn_all = np.array([ds.solve_Q(V[i], Z[i] / (2*s)**2, rcond=rc) for i in range(N)])
    Qs_all = np.array([ds.solve_Q(B[i], Z[i], rcond=rc) for i in range(N)])
    dn = CR.fk_exact(grid, Cov, rho, rp, Qn_all)
    dsm = CR.fk_exact(grid, Cov, rho, rp, Qs_all)
    print(f"rcond={rc:8.1e}  nominal tot={dn['total']:+.5e}  smeared tot={dsm['total']:+.5e}")
