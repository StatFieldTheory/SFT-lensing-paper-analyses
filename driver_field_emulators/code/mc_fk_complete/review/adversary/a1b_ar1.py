"""A1.2 redone: the analytic R[k,p] used by fk_expect_exact vs a direct simulation
of the chain fk_complete_core actually runs.  Deviations are normalised by the
LOCAL variance scale |A[min(k,p)]| so that exponentially small far-lag entries do
not produce meaningless ratios."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex

rng = np.random.default_rng(11)
N = 60
grid = ex.build_grid(n_lambda=N)
print(f"N={N}  dlam={grid.dlam:.3f}")
V = np.empty((N, 6, 6))
for k in range(N):
    B = rng.standard_normal((6, 6)) * (1.0 + 0.5 * np.sin(0.3 * k))
    V[k] = B @ B.T + 0.4 * np.eye(6)

for sig in (100.0, 300.0):
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V); A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    Lf = np.array([np.linalg.cholesky(V[k]) for k in range(N)])
    m = 600_000
    s1mr = np.sqrt(1 - rho * rho)
    zs = np.empty((N, 6, m)); z = np.zeros((6, m))
    for k in range(N):
        innov = Lf[k] @ rng.standard_normal((6, m))
        z = innov if k == 0 else rho * z + s1mr * innov
        zs[k] = z
    print(f"\n  sigma={sig:6.1f}  rho={rho:.5f}   (MC noise ~ {1/np.sqrt(m):.4f})")
    print(f"    {'(k,p)':>10} {'|emp-ana|/|A_min|':>20} {'ctrl A->V':>12} {'ctrl rho^2lag':>14}")
    for (k, p) in ((0, 0), (5, 5), (20, 18), (40, 37), (59, 59), (3, 6), (30, 33), (50, 45)):
        emp = zs[k] @ zs[p].T / m
        scale = np.max(np.abs(A[min(k, p)]))
        ana = rho ** abs(k - p) * A[min(k, p)]
        bad1 = rho ** abs(k - p) * V[min(k, p)]
        bad2 = rho ** (2 * abs(k - p)) * A[min(k, p)]
        print(f"    {str((k,p)):>10} {np.max(np.abs(emp-ana))/scale:20.5f}"
              f" {np.max(np.abs(emp-bad1))/scale:12.5f} {np.max(np.abs(emp-bad2))/scale:14.5f}")
