"""FLRW background quantities for the Monte-Carlo Sachs integrator.

Thin wrapper over the project's closed-form ``D = a * chi`` background
(``scripts/sachs_sft/scripts/D_callable.py``).  The Sachs fluctuation ODE
integrated by the MC is

    ds_a/dlambda = -2 theta_sa(lambda) s_a + F_abc s_b s_c + f_a(n, lambda)

whose linear (F=0) Green's function is identically the project's

    response_R(t, lambda') = Theta(t - lambda') [D(lambda') / D(t)]^2,

because d/dt [D(lam')/D(t)]^2 = -2 (D'/D) [.]^2 = -2 theta_sa(t) [.]^2.
So this module exposes theta_sa(lambda) = d ln D / dlambda = D'/D, plus
pass-throughs to D, chi, z, a and the closed-form response for validation.

Cosmology is hardcoded in D_callable (flat LCDM, Omega_m=0.31609, h=0.6711).
All lengths are physical Mpc; lambda = 0 at the observer, increasing toward
the source.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

# --- locate and import the project's closed-form background --------------
_D_CALLABLE_DIR = (
    Path(__file__).resolve().parents[1].parent / "scripts"  # .../sachs_sft/scripts
)
if str(_D_CALLABLE_DIR) not in sys.path:
    sys.path.insert(0, str(_D_CALLABLE_DIR))

import D_callable as _Dc  # noqa: E402  (D_at, chi_of_lambda, z_of_lambda, response_R, ...)

# lambda at the fiducial source redshift z_s = 5 (physical Mpc).
# D_callable gives 2317.96; the sft-wick production baseline npz used
# t_final = 2313.029.  Default to the D_callable value; callers that want
# bit-for-bit comparison with the baseline can pass LAM_SOURCE_BASELINE.
LAM_SOURCE = float(_Dc.lambda_of_z(5.0)) if hasattr(_Dc, "lambda_of_z") else 2317.961
LAM_SOURCE_BASELINE = 2313.029


class Background:
    """Tabulated FLRW background with a smooth theta_sa(lambda) = D'/D.

    Parameters
    ----------
    lam_min, lam_max : float
        Spline domain in physical Mpc.  Defaults span the union of the
        driving-field cumulant tables (limber C: 210..2318, zeta: 397..2327).
    n_grid : int
        Number of spline nodes (log-spaced away from the lam=0 singularity).
    """

    def __init__(self, lam_min: float = 150.0, lam_max: float = 2330.0,
                 n_grid: int = 2048):
        self.lam_min = float(lam_min)
        self.lam_max = float(lam_max)
        # log-spaced: theta_sa ~ 1/lambda diverges near 0, so resolve small lam.
        lam = np.geomspace(lam_min, lam_max, n_grid)
        D = np.asarray(_Dc.D_at(lam), dtype=float)
        lnD = np.log(D)
        self._lam = lam
        self._D_spline = CubicSpline(lam, D)
        self._lnD_spline = CubicSpline(lam, lnD)
        self._theta_spline = self._lnD_spline.derivative(1)  # d lnD/dlam = D'/D

    # --- background fields -------------------------------------------------
    def D(self, lam):
        """Background Jacobi amplitude D = a * chi (physical Mpc)."""
        return self._D_spline(lam)

    def theta_sa(self, lam):
        """Saddle-point expansion theta_sa = D'/D = d ln D / dlambda (1/Mpc)."""
        return self._theta_spline(lam)

    def z(self, lam):
        return _Dc.z_of_lambda(lam)

    def chi(self, lam):
        return _Dc.chi_of_lambda(lam)

    def a(self, lam):
        return 1.0 / (1.0 + np.asarray(_Dc.z_of_lambda(lam), dtype=float))

    def response_R(self, t, lam_prime):
        """Closed-form response [D(lam')/D(t)]^2 * Theta(t-lam') (validation)."""
        return _Dc.response_R(t, lam_prime)

    def response_R_from_spline(self, t, lam_prime):
        """Same response built from this object's D spline (self-consistency)."""
        t = np.asarray(t, dtype=float)
        lp = np.asarray(lam_prime, dtype=float)
        ratio = self._D_spline(lp) / self._D_spline(t)
        return np.where(lp <= t, ratio * ratio, 0.0)


if __name__ == "__main__":
    # Quick sanity check: theta_sa = D'/D should match response_R structure,
    # and the spline response should agree with D_callable.response_R.
    bg = Background()
    lam = np.array([300.0, 800.0, 1500.0, 2200.0])
    print("lambda      :", lam)
    print("D = a*chi    :", bg.D(lam))
    print("z(lambda)    :", np.asarray(bg.z(lam)))
    print("theta_sa=D'/D:", bg.theta_sa(lam))
    print("1/lambda ref :", 1.0 / lam, "(theta_sa ~ 1/lam near observer)")
    # response self-consistency
    t = LAM_SOURCE
    r_spline = bg.response_R_from_spline(t, lam)
    r_exact = np.asarray([bg.response_R(t, x) for x in lam])
    print("R spline     :", r_spline)
    print("R D_callable :", r_exact)
    print("max rel diff :", np.max(np.abs(r_spline - r_exact) / np.abs(r_exact)))
    print("LAM_SOURCE(z=5):", LAM_SOURCE)
