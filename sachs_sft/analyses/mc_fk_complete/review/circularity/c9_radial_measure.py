"""The historically REAL bug, run through the comparison: the (1+z)^-4 measure.

On 2026-06-09 the equal-shell kappa3 vertex was found to be missing a
(1+z)^-4 radial-measure Jacobian; fixing it changed the FK kappa-kappa
baseline by a factor 0.035 (+1.219e-4 -> +4.29e-6).  driver_stats' own
docstring says the per-leg (1+z)^4 Sachs prefactor is "already baked in" to the
tabulated vertex -- i.e. that factor lives INSIDE the shared table.

So: apply it to the estimator side only, to measure how large the error is;
then note that because it is baked into the object BOTH sides read, the ratio
is invariant and the comparison is blind to it.
"""
import _env, numpy as np
import fk_expect_exact as ex
from run_gate2 import analytic_fk

N, GAM = 1000, 1.0
grid = ex.build_grid(n_lambda=N)
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
ana = float(analytic_fk([GAM])[0])
zz = ex._DS.zeta6
opz = 1.0 + np.asarray(grid.bg.z(grid.lam), float)
print(f"(1+z) over the grid: {opz.min():.3f} .. {opz.max():.3f}")


def run(fac_fn, label):
    ex._DS.zeta6 = (lambda cg, l: fac_fn(l) * zz(cg, l))
    r = {}
    for sig in (8.0, 4.0):
        V, A, Q, rho = ex.node_stats(grid, cosg, sig)
        T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
        r[sig] = float(((T1 + T2) + (T1 + T2).T)[0, 3]) / ana
    ex._DS.zeta6 = zz
    print(f"{label:<52}{r[8.0]:>9.4f}{r[4.0]:>9.4f}{2*r[4.0]-r[8.0]:>10.4f}")


zof = grid.bg.z
print(f"\n{'zeta modification (ESTIMATOR SIDE ONLY)':<52}{'s=8':>9}{'s=4':>9}"
      f"{'r(s->0)':>10}")
run(lambda l: 1.0, "none")
run(lambda l: (1.0 + float(zof(l))) ** -4, "zeta *= (1+z)^-4   (the 2026-06-09 bug)")
run(lambda l: (1.0 + float(zof(l))) ** -2, "zeta *= (1+z)^-2")
run(lambda l: (1.0 + float(zof(l))) ** +2, "zeta *= (1+z)^+2")
print("\nBUT: driver_stats.zeta_tensor -- 'lambda-density, per-leg (1+z)^4 Sachs")
print("prefactor already baked in' -- so this factor is INSIDE the shared table.")
print("Applied to both sides the ratio is invariant EXACTLY (estimator is linear")
print("in zeta node-by-node; the order-2 fold carries one K insertion, also")
print("linear).  A factor of this size is therefore invisible to the test.")
