"""Spin-2 angular correlation kernels for the (Re Psi_0, Im Psi_0) sector.

A spin-2 field on the sphere transforms as ``(Q + i U) -> e^{-2 i psi}
(Q + i U)`` under a passive rotation of the local tangent basis by angle
``psi``. For two points (n_1, n_2) at angular separation ``gamma``, working
in the *great-circle* frame (local x-axis at each point along the connecting
great circle) the two-point correlations are isotropic and given by

    xi_plus(gamma)  = (1/4 pi) sum_l (2 l + 1) [C_l^EE + C_l^BB] d^l_{2,2}(gamma)
    xi_minus(gamma) = (1/4 pi) sum_l (2 l + 1) [C_l^EE - C_l^BB] d^l_{2,-2}(gamma)
    xi_E_scalar(g)  = (1/4 pi) sum_l (2 l + 1) C_l^{T E}        d^l_{0,2}(gamma)

where d^l_{m m'} is the Wigner small-d. In linear scalar Phi=Psi theory
B-modes vanish, so xi_plus and xi_minus both reduce to a single E-only
spectrum and ``xi_E_scalar`` carries the (Phi_00, Psi_0) cross-correlation.

In a *global* tangent basis (the canonical (e_theta, e_phi) at each
n_i), the correlations are obtained by rotating each leg by its own
angle ``psi_i`` relative to the great-circle frame. The closed-form
relations used by ``kappa2_callable.py`` are (see component_to_psi0.wl
for the symbolic derivation):

    <Q_1 Q_2> = (1/2)[ cos(2(p1 - p2)) xi_p + cos(2(p1 + p2)) xi_m ]
    <U_1 U_2> = (1/2)[ cos(2(p1 - p2)) xi_p - cos(2(p1 + p2)) xi_m ]
    <Q_1 U_2> = (1/2)[ sin(2(p1 + p2)) xi_m - sin(2(p1 - p2)) xi_p ]
    <U_1 Q_2> = (1/2)[ sin(2(p1 + p2)) xi_m + sin(2(p1 - p2)) xi_p ]

    <T_1 Q_2> = cos(2 p2) xi_TE,    <T_1 U_2> = sin(2 p2) xi_TE
    <Q_1 T_2> = cos(2 p1) xi_TE,    <U_1 T_2> = sin(2 p1) xi_TE

with ``T = Phi_00`` (spin-0), ``Q = Re Psi_0`` and ``U = Im Psi_0`` in the
global basis at each point. See Zaldarriaga & Seljak 1997 (ApJ 482 28),
Ng & Liu 1999, and Smith et al. 2007 for the conventions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import lpmv


@dataclass(frozen=True)
class GreatCircleAngles:
    """Per-pair spin-rotation angles relative to the global (e_theta, e_phi) basis."""

    cos_gamma: float
    sin_gamma: float
    psi_1: float  # angle at n_1 from e_theta_1 to t_12 (CCW from outside)
    psi_2: float  # angle at n_2 from e_theta_2 to t_21


def _to_unit(n: np.ndarray) -> np.ndarray:
    n = np.asarray(n, dtype=float)
    norm = np.linalg.norm(n)
    if norm == 0.0:
        raise ValueError("Direction vector has zero norm")
    return n / norm


def _theta_phi_basis(n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (e_theta, e_phi) at the unit vector n in the standard convention.

    e_theta points towards increasing polar angle theta (south pole when at
    equator); e_phi points towards increasing azimuth phi. At the poles, phi
    is degenerate; we pick a convention where e_phi -> (0, 1, 0) at z = +1
    and (0, -1, 0) at z = -1 to keep the basis right-handed.
    """
    x, y, z = n
    rho = np.hypot(x, y)
    if rho < 1e-12:
        # Pole: pick a deterministic frame.
        sign = 1.0 if z > 0 else -1.0
        e_theta = np.array([sign, 0.0, 0.0])
        e_phi = np.array([0.0, sign, 0.0])
        return e_theta, e_phi
    cos_phi = x / rho
    sin_phi = y / rho
    sin_theta = rho
    cos_theta = z
    e_theta = np.array([cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta])
    e_phi = np.array([-sin_phi, cos_phi, 0.0])
    return e_theta, e_phi


def great_circle_angles(n1: np.ndarray, n2: np.ndarray) -> GreatCircleAngles:
    """Compute (cos gamma, sin gamma, psi_1, psi_2) for the spin-2 rotation.

    psi_i is the angle, measured CCW as seen from outside the sphere, from
    the global ``e_theta_i`` axis to the unit tangent ``t_{ij}`` (the
    tangent at n_i to the great circle pointing towards n_j).
    """
    n1 = _to_unit(n1)
    n2 = _to_unit(n2)
    cos_gamma = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
    sin_gamma = float(np.sqrt(max(0.0, 1.0 - cos_gamma ** 2)))

    if sin_gamma < 1e-12:
        # Coincident or antipodal: psi is undefined; great-circle frame is
        # degenerate. Return zeros and let downstream code use the scalar
        # limit (xi_+ at gamma = 0, etc.).
        return GreatCircleAngles(cos_gamma, sin_gamma, 0.0, 0.0)

    # Tangent at n_1 to great circle towards n_2:
    t12 = (n2 - cos_gamma * n1) / sin_gamma
    t21 = (n1 - cos_gamma * n2) / sin_gamma

    e_th1, e_ph1 = _theta_phi_basis(n1)
    e_th2, e_ph2 = _theta_phi_basis(n2)
    psi_1 = float(np.arctan2(np.dot(t12, e_ph1), np.dot(t12, e_th1)))
    psi_2 = float(np.arctan2(np.dot(t21, e_ph2), np.dot(t21, e_th2)))
    return GreatCircleAngles(cos_gamma, sin_gamma, psi_1, psi_2)


def _wigner_d(l: int, m: int, mp: int, cos_beta: float) -> float:
    """Wigner small-d ``d^l_{m m'}(beta)`` via the spin-spherical-harmonic
    relation ``d^l_{m m'}(beta) = (-)^{m'} sqrt((l+m)!(l-m)!/((l+m')!(l-m')!))
    P_l^{m, m'}`` is awkward; we use the closed Jacobi-polynomial form.

    For our needs only ``(m, m') in {(0, 2), (2, 2), (2, -2)}`` appear, so
    we hard-code those cases via lpmv (associated Legendre) plus the
    inversion formula. Specifically, with ``x = cos(beta)``:
        d^l_{2, 2}(beta) = ((l - 2)! / (l + 2)!) * P_l^2(x) * <prefactor>
    and similarly for d^l_{2, -2}, d^l_{0, 2}.

    A more numerically stable approach (used here) is via the explicit
    formula in terms of associated Legendre functions:

        d^l_{m m'}(beta) = (-1)^{m - m'} sqrt(((l+m')! (l-m')!)/((l+m)! (l-m)!)) *
                           [(1 - x)/2]^{(m - m')/2} *
                           [(1 + x)/2]^{(m + m')/2} *
                           P_{l-m}^{(m - m', m + m')}(x)

    For our restricted (m, m') we expand directly. We require l >= max(|m|, |m'|).
    """
    # Reduce to the canonical (m, m') with m >= m' >= 0 if possible via:
    #   d^l_{m, m'}(beta) = (-1)^{m - m'} d^l_{m', m}(beta)
    #   d^l_{-m, -m'}(beta) = (-1)^{m - m'} d^l_{m', m}(beta)
    # (These are the standard symmetry relations.)
    sign = 1.0
    if abs(mp) > abs(m):
        m, mp = mp, m
        sign *= (-1.0) ** (m - mp)
    if m < 0:
        m, mp = -m, -mp
        sign *= (-1.0) ** (m - mp)

    # Now |m| >= |m'|, m >= 0.
    if l < max(abs(m), abs(mp)):
        return 0.0

    # Use the Jacobi-polynomial form, which is numerically stable for large l.
    # Closed-form sum with t-cancellations overflows once individual terms
    # exceed float range (around l ~ 170 for naive factorials, and around
    # l ~ 600 even with log-space evaluation due to term cancellation).
    #
    # Wikipedia ("Wigner D-matrix"), with M = max(|m|, |mp|) and
    # N = min(|m|, |mp|):
    #   d^l_{m, mp}(beta) = sqrt[(l + M)! (l - M)! / ((l + N)! (l - N)!)]
    #                       * sin(beta/2)^{|m - mp|}
    #                       * cos(beta/2)^{|m + mp|}
    #                       * P^{(|m - mp|, |m + mp|)}_{l - M}(cos beta).
    # In our canonical form (m >= |mp|, m >= 0): M = m, N = |mp|, so
    # (l + m)! (l - m)! sits in the numerator. The factorial sum is symmetric
    # under mp -> -mp, so lgamma(l + mp + 1) + lgamma(l - mp + 1) is
    # equivalent to lgamma(l + |mp| + 1) + lgamma(l - |mp| + 1).
    from math import lgamma, exp, log, cos, sin, acos
    from scipy.special import eval_jacobi

    beta = acos(max(-1.0, min(1.0, float(cos_beta))))
    cos_half = cos(0.5 * beta)
    sin_half = sin(0.5 * beta)

    log_pref = 0.5 * (
        lgamma(l + m + 1) + lgamma(l - m + 1)
        - lgamma(l + mp + 1) - lgamma(l - mp + 1)
    )
    pref = exp(log_pref)

    # Both exponents are non-negative since m >= |mp|.
    cos_pow = cos_half ** (m + mp) if (m + mp) > 0 else 1.0
    sin_pow = sin_half ** (m - mp) if (m - mp) > 0 else 1.0

    P = float(eval_jacobi(l - m, m - mp, m + mp, float(cos_beta)))
    return sign * pref * cos_pow * sin_pow * P


def _wigner_d_array(l_arr: np.ndarray, m: int, mp: int, cos_beta: float) -> np.ndarray:
    """Vectorised Wigner-d for a 1-D array of l values, fixed (m, mp, cos_beta).

    Replaces a Python loop over ``_wigner_d`` with one ``scipy.special.eval_jacobi``
    call on the whole ell array. Same canonical-form symmetry reductions as
    ``_wigner_d``; output is zero for any ell below ``max(|m|, |mp|)``.
    """
    from math import lgamma, exp, cos, sin, acos
    from scipy.special import eval_jacobi

    l_arr = np.asarray(l_arr, dtype=int)
    m_in, mp_in = int(m), int(mp)
    sign = 1.0
    if abs(mp_in) > abs(m_in):
        m_in, mp_in = mp_in, m_in
        sign *= (-1.0) ** (m_in - mp_in)
    if m_in < 0:
        m_in, mp_in = -m_in, -mp_in
        sign *= (-1.0) ** (m_in - mp_in)

    valid = l_arr >= max(abs(m_in), abs(mp_in))
    out = np.zeros(l_arr.shape, dtype=float)
    if not np.any(valid):
        return out
    l_v = l_arr[valid]

    beta = acos(max(-1.0, min(1.0, float(cos_beta))))
    cos_half = cos(0.5 * beta)
    sin_half = sin(0.5 * beta)
    cos_pow = cos_half ** (m_in + mp_in) if (m_in + mp_in) > 0 else 1.0
    sin_pow = sin_half ** (m_in - mp_in) if (m_in - mp_in) > 0 else 1.0

    # log-prefactor sqrt[(l+m)!(l-m)! / ((l+mp)!(l-mp)!)] per ell, vectorised
    # via lgamma of arrays. (Canonical form m >= |mp|, m >= 0; sum is symmetric
    # under mp -> -mp, so the explicit sign of mp does not matter.)
    from scipy.special import gammaln
    log_pref = 0.5 * (
        gammaln(l_v + m_in + 1) + gammaln(l_v - m_in + 1)
        - gammaln(l_v + mp_in + 1) - gammaln(l_v - mp_in + 1)
    )
    pref = np.exp(log_pref)

    P = eval_jacobi(l_v - m_in, m_in - mp_in, m_in + mp_in, float(cos_beta))
    out[valid] = sign * pref * cos_pow * sin_pow * P
    return out


def xi_plus_minus(
    cos_gamma: float, ells: np.ndarray, Cl_EE: np.ndarray, Cl_BB: np.ndarray | None = None
) -> tuple[float, float]:
    """Return (xi_plus(gamma), xi_minus(gamma)) for the spin-2 power spectra.

    ``ells`` is sorted; ``Cl_EE`` and ``Cl_BB`` (optional) have the same
    shape. For Phi=Psi linear theory, omit ``Cl_BB`` (taken as zero).
    """
    ells = np.asarray(ells, dtype=int)
    Cl_EE = np.asarray(Cl_EE, dtype=float)
    if Cl_BB is None:
        Cl_BB = np.zeros_like(Cl_EE)
    else:
        Cl_BB = np.asarray(Cl_BB, dtype=float)
    weight = (2.0 * ells.astype(float) + 1.0) / (4.0 * np.pi)
    d22 = _wigner_d_array(ells, 2, 2, cos_gamma)
    d2m2 = _wigner_d_array(ells, 2, -2, cos_gamma)
    xi_p = float(np.sum(weight * (Cl_EE + Cl_BB) * d22))
    xi_m = float(np.sum(weight * (Cl_EE - Cl_BB) * d2m2))
    return xi_p, xi_m


def xi_TE_scalar(cos_gamma: float, ells: np.ndarray, Cl_TE: np.ndarray) -> float:
    """Return the scalar-spin-2 cross-correlation in the great-circle frame.

    xi_TE(gamma) = (1/4 pi) sum_l (2 l + 1) C_l^{T E} d^l_{0, 2}(gamma).
    """
    ells = np.asarray(ells, dtype=int)
    Cl_TE = np.asarray(Cl_TE, dtype=float)
    weight = (2.0 * ells.astype(float) + 1.0) / (4.0 * np.pi)
    d02 = _wigner_d_array(ells, 0, 2, cos_gamma)
    return float(np.sum(weight * Cl_TE * d02))


def spin2_real_block(
    xi_plus: float,
    xi_minus: float,
    psi_1: float,
    psi_2: float,
) -> np.ndarray:
    """Real-space (Q, U)-(Q, U) correlation 2x2 block in the global frame.

    Returns ``M[a, b] = <field_a(n_1) field_b(n_2)>`` with
    ``a, b in {Q, U}``.
    """
    cp_diff = np.cos(2.0 * (psi_1 - psi_2))
    cp_sum = np.cos(2.0 * (psi_1 + psi_2))
    sp_diff = np.sin(2.0 * (psi_1 - psi_2))
    sp_sum = np.sin(2.0 * (psi_1 + psi_2))
    QQ = 0.5 * (cp_diff * xi_plus + cp_sum * xi_minus)
    UU = 0.5 * (cp_diff * xi_plus - cp_sum * xi_minus)
    QU = 0.5 * (sp_sum * xi_minus - sp_diff * xi_plus)
    UQ = 0.5 * (sp_sum * xi_minus + sp_diff * xi_plus)
    return np.array([[QQ, QU], [UQ, UU]])


def scalar_spin2_pair(xi_TE: float, psi_other: float) -> tuple[float, float]:
    """Cross-correlation <T(n_T) Q(n_other)> and <T(n_T) U(n_other)>."""
    return float(np.cos(2.0 * psi_other) * xi_TE), float(np.sin(2.0 * psi_other) * xi_TE)


__all__ = [
    "GreatCircleAngles",
    "great_circle_angles",
    "scalar_spin2_pair",
    "spin2_real_block",
    "xi_TE_scalar",
    "xi_plus_minus",
]
