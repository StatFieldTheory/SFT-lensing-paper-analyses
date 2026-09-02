"""The SECOND shared input the NOTES caveat does not mention: D(lambda).

``background.py`` (Monte-Carlo) is a thin wrapper over exactly the
``scripts/D_callable.py`` that the sft-wick config names as its
``R_time_module``.  The F vertex is shared too (``scripts/inputs/F_tensor.npy``
is numerically identical to the hard-coded ``sachs_mc_core._F``).  So the
comparison is blind to errors in D and in F as well as in zeta6.

This measures the size of that blind spot: how far would the ratio move if the
warp were applied to ONE side only?  That is exactly how much a COMMON error of
the same size hides.
"""
import _env, numpy as np
import fk_expect_exact as ex, inj_exact as ij
from run_gate2 import analytic_fk

N = 1000
GAM = 1.0
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
ana = float(analytic_fk([GAM])[0])
base = ex.build_grid(n_lambda=N)
print(f"{'D warp (estimator side only)':<40}{'s=8':>9}{'s=4':>9}{'r(s->0)':>10}")
for mode, eps in (("none", 0.0), ("const", 0.02), ("const", 0.05),
                  ("linear", 0.01), ("linear", 0.02), ("linear", 0.05),
                  ("linear", -0.05)):
    g = base if mode == "none" else ij.build_grid_warped(N, eps=eps, mode=mode)
    r = {}
    for sig in (8.0, 4.0):
        V, A, Q, rho = ex.node_stats(base, cosg, sig)
        T1, T2 = ij.expectation(g, cosg, sig, V, A, Q, rho)
        r[sig] = ij.fk_kk(T1, T2) / ana
    lab = "none" if mode == "none" else f"D -> D (1 + {eps:+.2f} u), u={mode}"
    print(f"{lab:<40}{r[8.0]:>9.4f}{r[4.0]:>9.4f}{2*r[4.0]-r[8.0]:>10.4f}")
