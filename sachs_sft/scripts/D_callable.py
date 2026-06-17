"""Canoes-pipeline Jacobi response for sft-wick FF runs.

This module intentionally has no ``canoes`` dependency because the FF runner is
executed in the PyCCL environment.  It mirrors the canoes-pipeline fiducial
cosmology and manuscript light-cone convention:

    d lambda / d chi = a(chi)^2 = 1 / (1 + z)^2,
    D(lambda) = a(lambda) * chi(lambda),

with lambda and chi in physical Mpc and observer normalization E_0 = 1.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline

C_LIGHT_KM_S: float = 299792.458
OMEGA_M: float = 0.3160919980475834
H_REDUCED: float = 0.6711
Z_MAX: float = 10.0
N_Z: int = 5001


def _tables() -> dict[str, np.ndarray | CubicSpline]:
    z = np.linspace(0.0, Z_MAX, N_Z)
    e_z = np.sqrt(OMEGA_M * (1.0 + z) ** 3 + (1.0 - OMEGA_M))
    hubble_radius_mpc = C_LIGHT_KM_S / (100.0 * H_REDUCED)

    chi = np.zeros_like(z)
    chi[1:] = cumulative_trapezoid(hubble_radius_mpc / e_z, z)

    lam = np.zeros_like(z)
    lam[1:] = cumulative_trapezoid(hubble_radius_mpc / ((1.0 + z) ** 2 * e_z), z)

    z_end = float(z[-1])
    one_plus_z_end_sq = (1.0 + z_end) ** 2
    return {
        "z": z,
        "chi": chi,
        "lambda": lam,
        "z_of_lambda": CubicSpline(
            lam,
            z,
            bc_type=((1, 1.0 / hubble_radius_mpc), (1, one_plus_z_end_sq * e_z[-1] / hubble_radius_mpc)),
        ),
        "chi_of_lambda": CubicSpline(
            lam,
            chi,
            bc_type=((1, 1.0), (1, one_plus_z_end_sq)),
        ),
    }


_CACHE = _tables()


def z_of_lambda(lam: np.ndarray | float) -> np.ndarray | float:
    """Return background redshift for physical affine parameter lambda in Mpc."""
    lam_arr = np.asarray(lam, dtype=float)
    z = _CACHE["z_of_lambda"](lam_arr)
    return float(z) if np.ndim(lam) == 0 else z


def chi_of_lambda(lam: np.ndarray | float) -> np.ndarray | float:
    """Return comoving distance for physical affine parameter lambda in Mpc."""
    lam_arr = np.asarray(lam, dtype=float)
    chi = _CACHE["chi_of_lambda"](lam_arr)
    return float(chi) if np.ndim(lam) == 0 else chi


def D_at(lam: np.ndarray | float) -> np.ndarray | float:
    """Jacobi amplitude D(lambda) = a(lambda) chi(lambda) in physical Mpc."""
    lam_arr = np.asarray(lam, dtype=float)
    z = _CACHE["z_of_lambda"](lam_arr)
    chi = _CACHE["chi_of_lambda"](lam_arr)
    out = chi / (1.0 + z)
    return float(out) if np.ndim(lam) == 0 else out


@lru_cache(maxsize=8192)
def _D_scalar_cached(lam: float) -> float:
    return float(D_at(lam))


def _D_scalar(lam: float) -> float:
    return _D_scalar_cached(round(lam, 8))


def theta_sa_integral(
    lam_upper: np.ndarray | float,
    lam_lower: np.ndarray | float,
) -> np.ndarray | float:
    """Return integral theta^(sa) d lambda via log D endpoints."""
    return np.log(D_at(lam_upper)) - np.log(D_at(lam_lower))


def response_R(
    t: np.ndarray | float,
    lam_prime: np.ndarray | float,
) -> np.ndarray | float:
    """Response R(t, lambda') = Theta(t-lambda') [D(lambda') / D(t)]^2."""
    if isinstance(t, (int, float, np.floating, np.integer)) and isinstance(
        lam_prime, (int, float, np.floating, np.integer)
    ):
        t_f = float(t)
        lp_f = float(lam_prime)
        if t_f < lp_f:
            return 0.0
        D_t = _D_scalar(t_f)
        if D_t == 0.0:
            return 0.0
        ratio = _D_scalar(lp_f) / D_t
        return ratio * ratio

    t_arr = np.asarray(t, dtype=float)
    lp_arr = np.asarray(lam_prime, dtype=float)
    D_t = D_at(t_arr)
    D_lp = D_at(lp_arr)
    safe_D_t = np.where(D_t > 0.0, D_t, 1.0)
    ratio2 = np.where(D_t > 0.0, (D_lp / safe_D_t) ** 2, 0.0)
    return np.where(t_arr >= lp_arr, ratio2, 0.0)


__all__ = [
    "D_at",
    "theta_sa_integral",
    "response_R",
    "z_of_lambda",
    "chi_of_lambda",
]
