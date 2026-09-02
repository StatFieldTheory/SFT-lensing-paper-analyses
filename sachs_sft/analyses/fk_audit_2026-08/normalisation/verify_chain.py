"""Step 1: does the pipeline's factor chain reduce to standard Limber?

Six independent numerical checks, each printed with a PASS/FAIL and the
worst relative residual:

  T1  D = a chi solves the background Jacobi equation in the affine
      parameter, with the paper's background Ricci focusing
      Phi_00^(0) = -(1+z)^4 (calH^2 - calH').
  T2  the affine lensing kernel built from the paper's response propagator,
      g = -D^2 int_lambda^{lambda_s} d lambda / D^2, evaluated by brute-force
      quadrature, equals the closed form -a^2 chi (chi_s-chi)/chi_s.
  T3  the assembled comoving weight q = (d lambda/d chi) g A (1+z)^4 equals
      the textbook lensing efficiency W, sign included.
  T4  the per-leg screen-Hessian multiplier cancels the Poisson 1/k^2 under
      Limber, leaving exactly A(a)(1+z)^4 (and the size of the l(l+1) vs
      l^2 residual that the code's k = l/chi choice leaves behind).
  T5  the equal-shell collapse Jacobian is (d lambda/d chi)^2 = (1+z)^-4,
      i.e. the n=3 analogue of the PyCCL-validated n=2 factor (1+z)^-2.
  T6  the full assembly of B_kappa from pipeline factors equals
      int d chi W^3 / chi^4 B_delta.

Run:
  CAN=/Users/zzhang/projects/angular_statistics/canoes
  PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python normalisation/verify_chain.py
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from canoes.cosmo.background import H_DEFAULT, OMEGA_M_DEFAULT
from canoes.cosmo.lightcone import build_lightcone
from canoes.sachs.kappa3 import (
    _H0_H_PER_MPC,
    _kappa3_radial_density_conversion,
)

from normalisation.chain import (
    SimpleCosmo,
    affine_lensing_kernel,
    affine_lensing_kernel_quadrature,
    b_zeta_lambda_density,
    background,
    d_lambda_d_chi,
    lensing_efficiency_from_chain,
    lensing_efficiency_standard,
    poisson_amp,
)

COSMO = SimpleCosmo(Omega_m=OMEGA_M_DEFAULT, h=H_DEFAULT)
RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str) -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def relerr(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)


# ---------------------------------------------------------------- T1
def t1_background_focusing() -> None:
    """Phi_00^(0) = A(a)(1+z)^4: the same amplitude the vertex uses per leg.

    The paper's background Ricci focusing is
    ``Phi_00^(0) = -(E_0^2 (1+z)^2/a^2)(calH^2 - calH')`` [Eq. driving bg].
    For flat LCDM ``calH^2 - calH' = 3/2 Omega_m H0^2 (1+z) = -A(a)``, so the
    background focusing is exactly the per-leg Limber response evaluated on a
    unit density contrast.  Checked against the canoes background splines.
    """
    chi = np.linspace(50.0, 4000.0, 2001)
    lc, opz = build_lightcone(COSMO, chi, distance_unit="Mpc/h")
    lhs = lc.calH**2 - lc.calHp                      # calHp = -d calH/d chi = calH'
    rhs = 1.5 * COSMO.Omega_m * _H0_H_PER_MPC**2 * opz
    err = relerr(lhs, rhs)
    phi00_bg = -(opz**4) * lhs
    leg = poisson_amp(opz, COSMO) * opz**4
    err_leg = relerr(phi00_bg, leg)
    record(
        "T1 Phi00^(0) = A(a)(1+z)^4 (background focusing = per-leg response)",
        float(max(err.max(), err_leg.max())) < 2e-4,
        f"calH^2-calH' vs 3/2 Om H0^2 (1+z): max rel diff {err.max():.2e} "
        f"(canoes spline accuracy); Phi00^(0) vs A(1+z)^4: {err_leg.max():.2e}",
    )


def t1b_jacobi_ode() -> None:
    """Integrate D'' = Phi_00^(0) D in lambda and recover D = a chi."""
    from scipy.integrate import solve_ivp

    chi_max = 4000.0

    def a_of_chi(c):
        c = np.atleast_1d(np.asarray(c, float))
        _lc, opz = build_lightcone(COSMO, np.clip(c, 1e-6, None),
                                   distance_unit="Mpc/h")
        return 1.0 / opz

    def rhs(c, y):
        a = float(a_of_chi(c)[0])
        opz = 1.0 / a
        phi00_bg = float(poisson_amp(np.array([opz]), COSMO)[0]) * opz**4
        # d/d chi = (d lambda/d chi) d/d lambda = a^2 d/d lambda
        return [a**2 * y[1], a**2 * phi00_bg * y[0]]

    chi_eval = np.linspace(1.0, chi_max, 200)
    sol = solve_ivp(rhs, (1e-6, chi_max), [0.0, 1.0], t_eval=chi_eval,
                    rtol=1e-10, atol=1e-14, method="DOP853")
    d_num = sol.y[0]
    d_exact = a_of_chi(chi_eval) * chi_eval
    err = relerr(d_num, d_exact)
    record(
        "T1b Jacobi ODE solution equals D = a chi",
        float(err.max()) < 1e-4,
        f"max rel diff {err.max():.2e} over chi in [1, {chi_max:.0f}] Mpc/h",
    )


# ---------------------------------------------------------------- T2
def t2_kernel() -> None:
    chi_s = 2000.0
    chi = np.array([50.0, 200.0, 600.0, 1000.0, 1500.0, 1900.0])
    g_closed = affine_lensing_kernel(COSMO, chi, chi_s)
    g_quad = affine_lensing_kernel_quadrature(COSMO, chi, chi_s, n_quad=40001)
    opz, _lam, _d = background(COSMO, chi)
    g_textbook = -(1.0 / opz**2) * chi * (chi_s - chi) / chi_s
    e1 = relerr(g_closed, g_quad)
    e2 = relerr(g_closed, g_textbook)
    record(
        "T2 affine lensing kernel g from the response propagator",
        float(max(e1.max(), e2.max())) < 1e-6,
        f"quadrature vs closed form {e1.max():.2e}; "
        f"closed form vs -a^2 chi(chi_s-chi)/chi_s {e2.max():.2e}",
    )


# ---------------------------------------------------------------- T3
def t3_efficiency() -> None:
    chi_s = 2000.0
    chi = np.linspace(20.0, 1980.0, 400)
    q = lensing_efficiency_from_chain(COSMO, chi, chi_s)
    w = lensing_efficiency_standard(COSMO, chi, chi_s)
    err = relerr(q, w)
    same_sign = bool(np.all(np.sign(q) == np.sign(w)))
    record(
        "T3 assembled q(chi) = (d lam/d chi) g A (1+z)^4 equals W(chi)",
        float(err.max()) < 1e-10 and same_sign,
        f"max rel diff {err.max():.2e}, signs agree: {same_sign}; "
        f"q/W median {np.median(q / w):.12f}",
    )


# ---------------------------------------------------------------- T4
def t4_limber_cancellation() -> None:
    """Screen-Hessian multiplier times Poisson 1/k^2, at k = l/chi."""
    chi = 1500.0
    ell = np.array([2.0, 10.0, 60.0, 100.0, 500.0, 1000.0, 4000.0])
    opz = np.array([2.0])
    amp = float(poisson_amp(opz, COSMO)[0])
    geom = float(opz[0]) ** 4

    # what the code does: k = l/chi, multiplier l^2/chi^2, potential 2A/k^2
    k_code = ell / chi
    leg_code = geom * (ell**2 / chi**2) * (2.0 * amp / k_code**2) / 2.0
    leg_target = np.full_like(ell, amp * geom)
    err_code = relerr(leg_code, leg_target)

    # what the exact spherical multiplier would give at the same k = l/chi:
    # Phi_00 multiplier L^2/chi^2 with L^2 = l(l+1) (appendix Eq. screen hessian)
    l2 = ell * (ell + 1.0)
    leg_exact = geom * (l2 / chi**2) * (2.0 * amp / k_code**2) / 2.0
    err_exact = relerr(leg_exact, leg_target)

    # spin-2 leg: sqrt(L^2 (L^2-2)) / chi^2
    leg_psi = geom * (np.sqrt(l2 * (l2 - 2.0)) / chi**2) * (2.0 * amp / k_code**2) / 2.0
    err_psi = relerr(np.abs(leg_psi), np.abs(leg_target))

    record(
        "T4 screen Hessian cancels Poisson 1/k^2 per leg -> A(a)(1+z)^4",
        float(err_code.max()) < 1e-14,
        f"exact in the code's own k = l/chi convention (max {err_code.max():.1e}); "
        f"the l(l+1) vs l^2 ambiguity is {err_exact.max()*100:.0f}% at l=2, "
        f"{err_exact[3]*100:.1f}% at l=100, {err_exact[-1]*100:.3f}% at l=4000; "
        f"spin-2 sqrt(L^2(L^2-2)) deficit {err_psi[3]*100:.2f}% at l=100",
    )


# ---------------------------------------------------------------- T5
def t5_collapse_jacobian() -> None:
    chi = np.array([200.0, 800.0, 2000.0])
    opz, _lam, _d = background(COSMO, chi)
    jac = d_lambda_d_chi(opz)
    fac3, name3, _d3 = _kappa3_radial_density_conversion(jac, "lambda")
    target3 = opz**-4
    err3 = relerr(fac3, target3)
    # n = 2 analogue used by the PyCCL-validated 2-point path
    target2 = opz**-2
    err2 = relerr(jac, target2)
    record(
        "T5 equal-shell collapse Jacobian (d lam/d chi)^2 = (1+z)^-4",
        float(max(err3.max(), err2.max())) < 1e-14,
        f"3-point factor '{name3}' matches (1+z)^-4 to {err3.max():.1e}; "
        f"2-point analogue matches (1+z)^-2 to {err2.max():.1e}",
    )


# ---------------------------------------------------------------- T6
def t6_bkappa_assembly() -> None:
    """Full B_kappa assembly from pipeline factors vs int W^3/chi^4 B_delta."""
    chi_s = 2000.0
    chi_grid = np.linspace(20.0, 1980.0, 601)
    rng = np.random.default_rng(7)
    # a stand-in B_delta(chi): the assembly is linear in it, so any
    # positive smooth profile exercises the measure bookkeeping.
    b_delta = 1e4 * (1.0 + 0.5 * np.sin(chi_grid / 300.0)) + rng.uniform(0, 1, chi_grid.size)

    opz, lam_grid, _d = background(COSMO, chi_grid)
    g = affine_lensing_kernel(COSMO, chi_grid, chi_s)

    # pipeline route: kernel g^3 folded against zeta as a LAMBDA density,
    # with the affine measure written out explicitly so both routes share
    # the same quadrature grid (the change of variable is not hidden).
    b_zeta_lam = np.array([
        b_zeta_lambda_density(COSMO, float(c), np.array([bd]))[0]
        for c, bd in zip(chi_grid, b_delta)
    ])
    pipeline = np.trapezoid(d_lambda_d_chi(opz) * g**3 * b_zeta_lam, chi_grid)

    # same thing but integrated on the pipeline's own lambda grid, to show
    # the change of variable is numerically consistent too
    pipeline_lam = np.trapezoid(g**3 * b_zeta_lam, lam_grid)

    # standard Limber route
    w = lensing_efficiency_standard(COSMO, chi_grid, chi_s)
    standard = np.trapezoid(w**3 / chi_grid**4 * b_delta, chi_grid)

    err = abs(pipeline - standard) / abs(standard)
    err_lam = abs(pipeline_lam - standard) / abs(standard)
    record(
        "T6 B_kappa from pipeline factors = int d chi W^3/chi^4 B_delta",
        err < 1e-12,
        f"pipeline/standard = {pipeline / standard:.14f} (rel diff {err:.2e}); "
        f"same fold on the lambda grid {pipeline_lam / standard:.8f} "
        f"(rel diff {err_lam:.1e}, quadrature only)",
    )


def main() -> None:
    print(f"cosmology: Omega_m={COSMO.Omega_m:.6f}, h={COSMO.h:.4f}, "
          f"H0={_H0_H_PER_MPC:.6e} h/Mpc\n")
    t1_background_focusing()
    t1b_jacobi_ode()
    t2_kernel()
    t3_efficiency()
    t4_limber_cancellation()
    t5_collapse_jacobian()
    t6_bkappa_assembly()
    n_fail = sum(1 for _n, ok, _d in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks pass")


if __name__ == "__main__":
    main()
