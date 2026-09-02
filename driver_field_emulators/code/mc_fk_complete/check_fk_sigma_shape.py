"""Is FK immune to the SHAPE of Sigma2(lambda), not merely to its overall scale?

This is the sharp form of "are the two sides given the same driving-field
statistics".  The FK estimator uses Sigma2 only through Q = solve_Q(B, zeta),
where Q ~ zeta/Sigma2^2, while each of the two propagator legs carries one power
of Sigma2.  If the cancellation is POINTWISE in lambda, then FK depends on zeta
and the background alone, and an arbitrary lambda-dependent distortion of Sigma2
leaves it unchanged.

That matters because the Monte-Carlo drives the Sachs equation with a
lambda-LOCAL field of density Sigma2, whereas the formalism consumes the full
non-local corr_op C(lam1, lam2).  Those differ (the equal-time collapse costs
12% in Order-0).  If FK is shape-invariant, that difference cannot touch it.
"""
import numpy as np
import _bootstrap
ds, _bgm, _core = _bootstrap.wire()
import fk_expect_exact as ex

GAM, SIG, N = 1.0, 8.0, 400
cosg = float(np.cos(np.deg2rad(GAM / 60.0)))


class Distorted(ds.Sigma2Builder):
    """Sigma2 -> g(lambda) * Sigma2, for an arbitrary smooth positive g."""
    g = staticmethod(lambda lam: 1.0)

    def matrix(self, cos_gamma, lam):
        return float(self.g(float(lam))) * super().matrix(cos_gamma, lam)


DISTORTIONS = [
    ("undistorted",              lambda l: 1.0),
    ("x 2 (global)",             lambda l: 2.0),
    ("linear ramp 0.5 -> 2.0",   lambda l: 0.5 + 1.5 * (l - 406.0) / 1907.0),
    ("tilt (lam/1200)^1.5",      lambda l: (l / 1200.0) ** 1.5),
    ("oscillatory 1+0.5 sin",    lambda l: 1.0 + 0.5 * np.sin(6.0 * np.pi
                                                              * (l - 406.0) / 1907.0)),
    ("step-like tanh 0.4->2.5",  lambda l: 0.4 + 2.1 / (1.0 + np.exp(-(l - 1300.0) / 60.0))),
]

print(f"gamma = {GAM}'   sigma_lambda = {SIG}   n_lambda = {N}")
print(f"{'distortion of Sigma2(lambda)':<28} {'Order-0':>13} {'O0 ratio':>9} "
      f"{'FK exact':>14} {'FK ratio':>10}")
base_fk = base_o0 = None
for name, gfun in DISTORTIONS:
    g = ex.build_grid(n_lambda=N)
    b = Distorted(background=g.bg, apply_c0=False)
    b.g = staticmethod(gfun)
    g.builder = b
    V, A, Q, rho = ex.node_stats(g, cosg, SIG, calibrate="smeared")
    T1, T2 = ex.expectation(g, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
    fk = ((T1 + T2) + (T1 + T2).T)[0, 3]
    o0 = float(ds.order0_mc(cosg, lam_f=g.lam[-1], builder=b, n_gauss=256))
    if base_fk is None:
        base_fk, base_o0 = fk, o0
    print(f"{name:<28} {o0:>13.5e} {o0/base_o0:>9.4f} {fk:>14.7e} "
          f"{fk/base_fk:>10.7f}")
