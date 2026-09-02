"""ADVERSARIAL TEST A1.

Three conventions the attack names, each tested so that a WRONG convention fails.

(1) F vertex and its linearisation.  fk_complete_core hand-codes _f_vertex6 and
    _f_vertex6_lin; fk_expect_exact uses tensors F6 and FS6 built from
    sachs_mc_core._F.  These are two SEPARATE implementations, so a mismatch is
    a bug in one of them.  In addition the linearisation is checked against a
    central finite difference of the vertex itself, which would catch a missing
    or spurious factor 2 / a missing symmetrisation even if BOTH sides shared it.

(2) The AR(1) chain.  fk_expect_exact models the simulated chain by
    A[k] = rho^2 A[k-1] + (1-rho^2) V[k],  R[k,p] = rho^|k-p| A[min(k,p)].
    That is an ANALYTIC MODEL of what fk_complete_core actually samples, so a
    direct simulation of the very same recursion is an independent check of it.

(3) The linear propagator.  s1[k] = sum_{p<=k} dlam (D_p/D_k)^2 z[p] is asserted
    by both sides; checked by running the recursion and reading the sensitivity.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc

rng = np.random.default_rng(7)
print("=" * 78)
print("A1.1  F vertex: MC hand-code vs exact-expectation tensor")
print("=" * 78)
s = rng.standard_normal((6, 5))
mc_F = fc._f_vertex6(s)
ten_F = np.einsum("abc,bm,cm->am", ex.F6, s, s, optimize=True)
e = np.max(np.abs(mc_F - ten_F)) / np.max(np.abs(ten_F))
print(f"  _f_vertex6(s)  vs  F6_abc s_b s_c        max rel dev = {e:.3e}")

d = rng.standard_normal((6, 5))
mc_L = fc._f_vertex6_lin(s, d)
ten_L = np.einsum("abc,bm,cm->am", ex.FS6, s, d, optimize=True)
e2 = np.max(np.abs(mc_L - ten_L)) / np.max(np.abs(ten_L))
print(f"  _f_vertex6_lin(s,d) vs FS6_abc s_b d_c   max rel dev = {e2:.3e}")

# central finite difference of the vertex -- independent of both encodings
for h in (1e-3, 1e-4, 1e-5):
    fd = (fc._f_vertex6(s + h * d) - fc._f_vertex6(s - h * d)) / (2 * h)
    e3 = np.max(np.abs(fd - mc_L)) / np.max(np.abs(mc_L))
    print(f"  d/dh F(s+h d)|_0  vs _f_vertex6_lin   h={h:.0e}  max rel dev = {e3:.3e}")

# FALSIFICATION CONTROLS -- these must NOT be ~0
bad1 = np.einsum("abc,bm,cm->am", ex.F6, s, d, optimize=True)          # forgot F+F^T
bad2 = 2.0 * ten_L                                                     # spurious 2
print(f"  control: FS -> F (no symmetrisation)     max rel dev = "
      f"{np.max(np.abs(bad1 - mc_L)) / np.max(np.abs(mc_L)):.3e}   (must be O(1))")
print(f"  control: FS -> 2 FS                      max rel dev = "
      f"{np.max(np.abs(bad2 - mc_L)) / np.max(np.abs(mc_L)):.3e}   (must be O(1))")

print()
print("=" * 78)
print("A1.2  AR(1) chain: direct simulation vs the analytic R used by the exact side")
print("=" * 78)
N = 60
grid = ex.build_grid(n_lambda=N)
# synthetic but non-stationary V[k] (the point is the recursion, not the table)
V = np.empty((N, 6, 6))
for k in range(N):
    B = rng.standard_normal((6, 6)) * (1.0 + 0.5 * np.sin(0.3 * k))
    V[k] = B @ B.T + 0.4 * np.eye(6)
for sig in (8.0, 30.0):
    rho = float(np.exp(-grid.dlam / sig))
    A = np.empty_like(V)
    A[0] = V[0]
    for k in range(1, N):
        A[k] = rho * rho * A[k - 1] + (1 - rho * rho) * V[k]
    Lf = np.array([np.linalg.cholesky(V[k]) for k in range(N)])
    m = 400_000
    s1mr = np.sqrt(1 - rho * rho)
    zs = np.empty((N, 6, m))
    z = np.zeros((6, m))
    for k in range(N):
        innov = Lf[k] @ rng.standard_normal((6, m))
        z = innov if k == 0 else rho * z + s1mr * innov
        zs[k] = z
    # analytic R for a few (k,p)
    worst = 0.0
    for (k, p) in ((0, 0), (5, 5), (17, 3), (40, 12), (59, 59), (3, 41), (30, 31)):
        emp = zs[k] @ zs[p].T / m
        ana = rho ** abs(k - p) * A[min(k, p)]
        r = np.max(np.abs(emp - ana)) / np.max(np.abs(ana))
        worst = max(worst, r)
    print(f"  sigma={sig:5.1f}  rho={rho:.6f}  worst rel dev over 7 (k,p) pairs = {worst:.4f}"
          f"   [MC noise ~ {1/np.sqrt(m):.4f}]")
    # falsification control: the stationary guess A -> V
    ctrl = 0.0
    for (k, p) in ((17, 3), (40, 12), (30, 31)):
        emp = zs[k] @ zs[p].T / m
        bad = rho ** abs(k - p) * V[min(k, p)]
        ctrl = max(ctrl, np.max(np.abs(emp - bad)) / np.max(np.abs(emp)))
    print(f"           control A -> V (stationary guess): rel dev = {ctrl:.4f}  (must be O(1))")

print()
print("=" * 78)
print("A1.3  linear propagator s1[k] = sum_p dlam (D_p/D_k)^2 z[p]")
print("=" * 78)
N = 200
grid = ex.build_grid(n_lambda=N)
D, resp, dlam = grid.D, grid.resp, grid.dlam
# run the recursion with a unit impulse at node p and read s1[k]
for p in (0, 37, 150):
    s1 = np.zeros(N)
    cur = 0.0
    for k in range(N):
        cur = resp[k] * cur + (dlam if k == p else 0.0)
        s1[k] = cur
    ks = np.array([p, p + 10, N - 1])
    ana = dlam * (D[p] / D[ks]) ** 2
    print(f"  impulse at p={p:4d}: max rel dev at k={list(ks)} = "
          f"{np.max(np.abs(s1[ks] - ana) / np.abs(ana)):.3e}")
    bad = dlam * (D[p] / D[ks]) ** 4
    print(f"      control exponent 4: rel dev = {np.max(np.abs(s1[ks]-bad)/np.abs(s1[ks])):.3e}")
