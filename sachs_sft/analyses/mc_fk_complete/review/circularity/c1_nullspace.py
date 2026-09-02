"""What EXACTLY cancels?  The null space of the comparison, measured.

The attack says the comparison is circular because both sides read the same
tabulated zeta6.  Granting that, the question is the DIMENSION of the null
space.  Measured here:

 1. the estimator is EXACTLY linear in the tabulated cumulant  ->  any overall
    rescaling of the table is invisible (both sides scale, ratio fixed);
 2. it is linear node-by-node, so an arbitrary lambda-reweighting g(lambda) of
    the table changes the estimator by exactly the same reweighting of its own
    kernel -- the ratio to the fold then moves only by the MISMATCH between the
    two kernels;
 3. so the comparison constrains  sum_m [K_MC(m) - K_fold(m)] zeta(m)  and
    nothing else.  Item 3 is what the cross-table and multi-source tests probe.
"""
import _env, numpy as np, fk_expect_exact as ex
from run_gate2 import analytic_fk

N = 1000
grid = ex.build_grid(n_lambda=N)
gam = 1.0
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
ana = float(analytic_fk([gam])[0])
_z0 = ex._DS.zeta6
zt = np.array([_z0(cosg, float(l)) for l in grid.lam])

print("=== 1. exact linearity of the ESTIMATOR in the tabulated cumulant ===")
base = None
for c in (1.0, 1.3, -2.7):
    ex._DS.zeta6 = (lambda cg, l, c=c: c * _z0(cg, l))
    V, A, Q, rho = ex.node_stats(grid, cosg, 8.0)
    T1, T2 = ex.expectation(grid, cosg, 8.0, V=V, A=A, Q=Q, rho=rho)
    v = float(((T1 + T2) + (T1 + T2).T)[0, 3])
    if base is None:
        base = v
    print(f"  zeta -> {c:+.2f} zeta : exact = {v: .10e}   v/(c*v0) - 1 = "
          f"{v / (c * base) - 1: .3e}")
ex._DS.zeta6 = _z0

print("\n  => the estimator is linear in zeta to machine precision.")
print("     The order-2 sft-wick fold carries exactly ONE K insertion, so it is")
print("     linear too.  A common rescaling of the table therefore cancels in")
print("     the ratio EXACTLY: the amplitude of zeta6 is 100% unconstrained.")

print("\n=== 2. the estimator's lambda kernel (what a table error is weighted by) ===")
kern = grid.dlam * grid.Wd * grid.Hd
integ = kern * np.einsum("Abc,mbcB->m", ex.F6, zt)[:]  # A=0,B=3 slice below
M = np.einsum("m,Abc,mbcB->mAB", kern, ex.F6, zt, optimize=True)
w = 2.0 * M[:, 0, 3]
w = w / w.sum()
c = np.cumsum(w)
q = [grid.lam[np.searchsorted(c, f)] for f in (0.05, 0.25, 0.5, 0.75, 0.95)]
print(f"  FK integrand quantiles in lambda [Mpc]: 5% {q[0]:.0f}  25% {q[1]:.0f}"
      f"  50% {q[2]:.0f}  75% {q[3]:.0f}  95% {q[4]:.0f}")
print(f"  grid spans [{grid.lam[0]:.0f}, {grid.lam[-1]:.0f}] Mpc")

print("\n=== 3. how much does the lambda weighting move with gamma? ===")
print("  (if it barely moves, the 9-point gamma sweep is ~1 constraint, not 9)")
qs = {}
for g in (0.5, 1.0, 5.0, 17.3):
    cg = float(np.cos(np.deg2rad(g / 60.0)))
    z = np.array([_z0(cg, float(l)) for l in grid.lam])
    Mg = np.einsum("m,Abc,mbcB->mAB", kern, ex.F6, z, optimize=True)
    wg = 2.0 * Mg[:, 0, 3]; wg = wg / wg.sum()
    cc = np.cumsum(wg)
    qs[g] = [grid.lam[np.searchsorted(cc, f)] for f in (0.25, 0.5, 0.75)]
    print(f"  gamma={g:5.1f}'  quartiles {qs[g][0]:.0f} / {qs[g][1]:.0f} / "
          f"{qs[g][2]:.0f} Mpc")
