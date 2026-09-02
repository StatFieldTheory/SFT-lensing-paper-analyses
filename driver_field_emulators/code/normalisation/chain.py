"""The pipeline's per-leg transfer chain, assembled into a lensing kernel.

Every factor here is taken from the canoes source, not re-derived:

* ``poisson_amp``  : ``kappa3.compute_kappa3_sigma3_high`` line
                     ``-1.5 * Omega_m * _H0_H_PER_MPC**2 * (1+z)``
* Sachs weight     : ``kappa3._kappa3_limber_response_product_h``
                     (``A(a) (1+z)^4`` per Phi00 leg, minus that per Psi0 leg)
* affine measure   : ``canoes.cosmo.lightcone._lambda_integrand`` /
                     ``kappa2._lambda_jacobian`` -> ``d lambda/d chi = (1+z)^-2``
* collapse jacobian: ``kappa3._kappa3_radial_density_conversion``
* background        : ``canoes.cosmo.lightcone.build_lightcone``

The one factor the vertex code does NOT contain is the affine lensing kernel
``g(lambda)``, i.e. the weight with which a driving field at affine distance
lambda enters the observed convergence.  That kernel belongs to the sft-wick
fold, and is built here from the paper's own response propagator
[Eq. "explicit resp op"], ``R(lambda' -> lambda) = (D(lambda')/D(lambda))^2``.

Units: everything in this module is h-native (chi, lambda in Mpc/h, H0 in
h/Mpc, P in (Mpc/h)^3, B in (Mpc/h)^6), exactly like the canoes internals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from canoes.cosmo.lightcone import build_lightcone, lambda_of_chi
from canoes.sachs.kappa2 import _lambda_jacobian
from canoes.sachs.kappa3 import (
    _H0_H_PER_MPC,
    _kappa3_limber_response_product_h,
    _kappa3_radial_density_conversion,
)

LAMBDA_CONVENTION = "project"


@dataclass(frozen=True)
class SimpleCosmo:
    """Minimal ``CosmoLike``: canoes only reads ``Omega_m`` and ``h``."""

    Omega_m: float
    h: float


def poisson_amp(one_plus_z: NDArray[np.float64], cosmo: SimpleCosmo):
    """``A(a) = -3/2 Omega_m H0^2 (1+z)`` in (h/Mpc)^2, canoes' expression."""
    return -1.5 * float(cosmo.Omega_m) * (_H0_H_PER_MPC**2) * np.asarray(
        one_plus_z, dtype=np.float64
    )


def d_lambda_d_chi(one_plus_z: NDArray[np.float64]) -> NDArray[np.float64]:
    """``d lambda / d chi = (1+z)^-2`` (E_0 = 1 observer normalisation)."""
    return _lambda_jacobian(np.asarray(one_plus_z, np.float64),
                            convention=LAMBDA_CONVENTION)


def background(cosmo: SimpleCosmo, chi_h: NDArray[np.float64]):
    """Return ``(one_plus_z, lambda_h, D_jacobi)`` on a chi grid in Mpc/h."""
    chi_h = np.asarray(chi_h, dtype=np.float64)
    _lc, one_plus_z = build_lightcone(cosmo, chi_h, distance_unit="Mpc/h")
    lam_h = lambda_of_chi(cosmo, chi_h, convention=LAMBDA_CONVENTION,
                          distance_unit="Mpc/h")
    d_jacobi = chi_h / one_plus_z          # D = a chi, the flat-FLRW Jacobi field
    return one_plus_z, lam_h, d_jacobi


def affine_lensing_kernel(
    cosmo: SimpleCosmo,
    chi_h: NDArray[np.float64],
    chi_s_h: float,
) -> NDArray[np.float64]:
    """``g(lambda)`` with ``kappa = int d lambda g(lambda) Phi_00(lambda)``.

    Built from the paper's response propagator only:

        delta theta(lambda) = int_0^lambda d lambda' R(lambda',lambda) Phi_00
        R(lambda', lambda) = (D(lambda')/D(lambda))^2
        kappa = - int_0^{lambda_s} d lambda  delta theta(lambda)

    so that, after exchanging the order of integration,

        g(lambda') = - D(lambda')^2 int_{lambda'}^{lambda_s} d lambda / D^2 .

    The integral is done exactly: with D = a chi and d lambda = a^2 d chi it
    collapses to int d chi / chi^2.
    """
    chi_h = np.asarray(chi_h, dtype=np.float64)
    one_plus_z, _lam, d_jac = background(cosmo, chi_h)
    tail = 1.0 / chi_h - 1.0 / float(chi_s_h)     # int_chi^chi_s d chi / chi^2
    return -(d_jac**2) * tail


def affine_lensing_kernel_quadrature(
    cosmo: SimpleCosmo,
    chi_h: NDArray[np.float64],
    chi_s_h: float,
    n_quad: int = 20001,
) -> NDArray[np.float64]:
    """Same kernel, but with the lambda integral done by brute-force quadrature.

    Independent of the analytic ``int d chi / chi^2`` collapse, so the two
    routes together check that collapse rather than assuming it.
    """
    chi_h = np.asarray(chi_h, dtype=np.float64)
    out = np.empty_like(chi_h)
    for i, chi_i in enumerate(chi_h):
        grid = np.linspace(float(chi_i), float(chi_s_h), n_quad)
        opz, _lam, d_jac = background(cosmo, grid)
        integrand = (1.0 / opz**2) / d_jac**2      # d lambda/d chi / D^2
        tail = np.trapezoid(integrand, grid)
        _opz_i, _l_i, d_i = background(cosmo, np.array([chi_i]))
        out[i] = -(d_i[0] ** 2) * tail
    return out


def lensing_efficiency_from_chain(
    cosmo: SimpleCosmo,
    chi_h: NDArray[np.float64],
    chi_s_h: float,
) -> NDArray[np.float64]:
    """``q(chi)``: the comoving-density weight of delta in kappa.

        kappa = int d lambda g Phi_00
              = int d chi (d lambda/d chi) g [A (1+z)^4] delta
              = int d chi q(chi) delta .

    Every factor comes from canoes: ``A`` from ``poisson_amp``, the
    ``(1+z)^4`` Sachs photon weight and the per-leg sign from
    ``_kappa3_limber_response_product_h`` (one Phi_00 leg), and the measure
    from ``_lambda_jacobian``.
    """
    chi_h = np.asarray(chi_h, dtype=np.float64)
    one_plus_z, _lam, _d = background(cosmo, chi_h)
    amp = poisson_amp(one_plus_z, cosmo)
    # Single Phi_00 leg, read out of the canoes helper itself: the TTT
    # response product is the cube of one Phi_00 leg, so cbrt recovers it
    # (sign included; A < 0 makes the product negative).
    leg = np.cbrt(np.array([
        _kappa3_limber_response_product_h("TTT", float(a), float(o))
        for a, o in zip(amp, one_plus_z)
    ]))
    g = affine_lensing_kernel(cosmo, chi_h, chi_s_h)
    return d_lambda_d_chi(one_plus_z) * g * leg


def lensing_efficiency_standard(
    cosmo: SimpleCosmo,
    chi_h: NDArray[np.float64],
    chi_s_h: float,
) -> NDArray[np.float64]:
    """Textbook ``W(chi) = 3/2 Omega_m H0^2 (1+z) chi (chi_s-chi)/chi_s``."""
    chi_h = np.asarray(chi_h, dtype=np.float64)
    one_plus_z, _lam, _d = background(cosmo, chi_h)
    return (1.5 * cosmo.Omega_m * _H0_H_PER_MPC**2 * one_plus_z
            * chi_h * (float(chi_s_h) - chi_h) / float(chi_s_h))


def b_zeta_lambda_density(
    cosmo: SimpleCosmo,
    chi_h: float,
    b_delta: NDArray[np.float64],
    channel: str = "TTT",
) -> NDArray[np.float64]:
    """The pipeline's flat-sky driving-field bispectrum, as a lambda density.

    ``compute_kappa3_sigma3_high`` builds, per shell,

        zeta_chi_fourier = R_X(z) * chi^-4 * B_delta(l_i/chi)

    then multiplies by ``_kappa3_radial_density_conversion`` (``(1+z)^-4``
    for ``radial_measure='lambda'``) and by ``h^4`` for physical units.
    This helper returns the h-native lambda-density value.
    """
    chi_arr = np.atleast_1d(np.asarray([chi_h], dtype=np.float64))
    one_plus_z, _lam, _d = background(cosmo, chi_arr)
    amp = poisson_amp(one_plus_z, cosmo)
    response = _kappa3_limber_response_product_h(
        channel, float(amp[0]), float(one_plus_z[0])
    )
    jac, _name, _desc = _kappa3_radial_density_conversion(
        d_lambda_d_chi(one_plus_z), "lambda"
    )
    return response * float(chi_h) ** -4 * np.asarray(b_delta) * float(jac[0])
