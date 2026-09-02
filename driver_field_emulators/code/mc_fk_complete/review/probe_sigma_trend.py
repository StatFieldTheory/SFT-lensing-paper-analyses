import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import cleanroom_exact as CR

N = 2000
cosg = float(np.cos(np.radians(1.0 / 60.0)))
grid = CR.build_grid(N)
S, Z = CR.driver_tables(cosg, grid["lam"])
loc = CR.fk_local(grid, Z)
print(f"N={N} dlam={grid['dlam']:.4f}   LOCAL ref total = {loc['total']:+.6e}")
print(f"{'sigma':>6} {'dlam/sig':>8} {'B/S2':>8} | {'nom T1':>12} {'nom tot':>12} "
      f"| {'sme T1':>12} {'sme tot':>12} {'sme/loc':>8}")
for s in (8.0, 4.0, 2.0, 1.0, 0.5):
    rho, V, Cov, rp = CR.ar1_stats(S, grid["dlam"], s)
    B = CR.smeared_B(Cov, rp, rho, grid["dlam"])
    Qn = CR.solve_Q_all(V, Z / (2.0 * s) ** 2)
    Qs = CR.solve_Q_all(B, Z)
    dn = CR.fk_exact(grid, Cov, rho, rp, Qn)
    dsm = CR.fk_exact(grid, Cov, rho, rp, Qs)
    r = float(np.median(B[:, 0, 0] / S[:, 0, 0]))
    print(f"{s:6.2f} {grid['dlam']/s:8.3f} {r:8.4f} | {dn['T1']:+12.5e} {dn['total']:+12.5e} "
          f"| {dsm['T1']:+12.5e} {dsm['total']:+12.5e} {dsm['total']/loc['total']:8.4f}")

# eigen-structure of Sigma2 / B at a mid shell
k = N // 2
rho, V, Cov, rp = CR.ar1_stats(S, grid["dlam"], 8.0)
B = CR.smeared_B(Cov, rp, rho, grid["dlam"])
for name, M in (("Sigma2", S[k]), ("V(s=8)", V[k]), ("B(s=8)", B[k]), ("Cov(s=8)", Cov[k])):
    d = np.linalg.eigvalsh(M)
    print(f"{name:9s} eig/max = {np.array2string(d/d.max(), precision=3)}")
