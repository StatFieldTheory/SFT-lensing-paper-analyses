"""Why "MC / exact expectation" is a gate that cannot fail.

`exact` is the closed-form Wick expectation of the estimator, computed from the
SAME discrete model the Monte-Carlo simulates: same propagator, same F vertex,
same Sigma2, same zeta.  So MC/exact asks only "does my sampler sample the model
I wrote down?"  It cannot notice that the model is wrong.

Demonstrated by breaking the model on purpose.  The Sachs response is
R(u,v) = [D(v)/D(u)]^2.  Substituting D -> sqrt(D) everywhere rebuilds resp, the
observable window Wd and the F-leg response Hd consistently with exponent 1
instead of 2 -- a different, WRONG physical model, but a perfectly
self-consistent one.  The driving-field statistics are left untouched.
"""
import numpy as np
import fk_expect_exact as ex
import fk_complete_core as fc
from run_gate2 import analytic_fk

GAM, SIG, N, SEEDS = 1.0, 8.0, 1000, 8
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))
ana = float(analytic_fk([GAM])[0])


def make_world(p):
    """Grid whose response exponent is p (production physics is p = 2)."""
    g = ex.build_grid(n_lambda=N)
    D = np.asarray(g.bg.D(g.lam), float) ** (p / 2.0)   # D -> D^(p/2)
    w = np.full(N, g.dlam); w[0] = w[-1] = 0.5 * g.dlam
    resp = np.empty(N); resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    Wd = D**2 * np.cumsum((w / D**2)[::-1])[::-1]
    gj = g.dlam * Wd[1:] / D[:-1] ** 4
    t4 = np.zeros(N); t4[:-1] = np.cumsum(gj[::-1])[::-1]
    g.D, g.resp, g.Wd, g.Hd = D, resp, Wd, D**4 * t4
    return g


print(f"gamma = {GAM}'   sigma_lambda = {SIG}   n_lambda = {N}   "
      f"{SEEDS} seeds x 24000")
print(f"sft-wick fold (the independent answer): {ana:.6e}\n")
print(f"{'world':>26} {'exact expectation':>19} {'exact/fold':>11} "
      f"{'MC':>13} {'MC/fold':>9} {'MC/exact':>10} {'+/-':>8}")
for p, label in ((2.0, "production, p = 2"), (1.0, "BROKEN, p = 1")):
    g = make_world(p)
    nodes = ex.node_stats(g, cosg, SIG, calibrate="smeared")
    V, A, Q, rho = nodes
    T1, T2 = ex.expectation(g, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
    exact = ((T1 + T2) + (T1 + T2).T)[0, 3]
    vals = [fc.simulate_fk_complete(GAM, sigma_lambda=SIG, n_lambda=N,
                                    n_real=24000, batch_size=2000,
                                    seed=90210 + 991 * s, grid=g, nodes=nodes).kk
            for s in range(SEEDS)]
    v = np.array(vals); sem = v.std(ddof=1) / np.sqrt(v.size)
    print(f"{label:>26} {exact:>19.6e} {exact/ana:>11.4f} {v.mean():>13.6e} "
          f"{v.mean()/ana:>9.4f} {v.mean()/exact:>10.4f} {sem/exact:>8.4f}")
