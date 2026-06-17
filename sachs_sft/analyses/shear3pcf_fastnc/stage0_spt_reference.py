"""STAGE 0 (reference side): independent hand-rolled SPT tree-level flat-sky
Limber convergence 3PCF on the COLLAPSED (squeezed) config at z_s = 5.

Interpreter (ABSOLUTE):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
(only needs numpy/scipy + canoes background; PyCCL env has canoes too).

PURPOSE
-------
A SECOND, fully independent implementation of the SAME tree+Limber+flat-sky
family as the canoes pipeline, with a DIFFERENT quadrature and an explicitly
hand-coded normalization, to cross-check the absolute normalization / radial
measure / spin-0 projection of the driving-field cumulant zeta_TTT.

It is an INTERNAL consistency check (same physics family + same P(k)), NOT an
approximation-independent test, so ~30% agreement is the expectation.

CONSTRUCTION (mirrors appendix `eq: appendix limber reduced bispectrum` +
`eq: appendix zeta squeezed` for the spin-0 TTT channel, hand-coded):

For each source shell lambda(chi), the squeezed (two coincident points, third at
gamma) convergence-3PCF integrand is the flat-sky Limber 3PCF of the screen-space
spin-0 Hessian of the lightcone potential:

  zeta_TTT(gamma, chi)
    = (1/(2pi)^4) * (D(z)^4 / chi^4) * [A(a) (1+z)^4]^3
      * INT du dv dphi  [u v]  B_delta(k1, k2, k3)  * 2pi J0(v * gamma),

with  k1 = w/chi,  k2 = u/chi,  k3 = v/chi,
      w = sqrt(u^2 + v^2 + 2 u v cos phi)   (the third triangle side),
      A(a) = -3/2 Omega_m H0^2 (1+z)        [Poisson, per leg],
      D(z) = linear growth (here we use the matter-dominated D ~ a normalisation
             via the canoes growth M[0]/(1+z) so the growth matches the LOW/HIGH
             tables exactly),
      B_delta = 2 F2 P P + 2 cyclic         [standard SPT, eq: appendix bdelta tree],
      F2 = 5/7 + (mu/2)(k_a/k_b + k_b/k_a) + (2/7) mu^2.

The spin-0 angular kernel is 2pi J0(v gamma): this is the squeezed limit of the
flat-sky orientation integral (the leg paired with k3=v/chi carries the gamma
separation; the two coincident points contribute no extra phase).  This is the
hand-derived analogue of canoes' alpha-kernel with spins=(0,0,0); we derive and
code it independently here.

FOLD (single source plane z_s=5; same convergence efficiency as ours-side):
  Z_ref(gamma) = INT_0^{lambda_s} dlambda K(lambda, lambda_s)^3 zeta_TTT(gamma, lambda),
  K(lambda, lambda_s) = a^2(chi) chi (chi_s - chi)/chi_s   [cosmology.tex eq: K comoving kernel].
We fold in chi: Z_ref = INT dchi a^8 K_chi^3 zeta_TTT(chi), K_chi = chi(chi_s-chi)/chi_s.

SIGN: B_delta>0 and the Poisson cube A(a)^3<0 (odd number of negative Poisson
factors) makes zeta_TTT<0, so Z_ref<0 too -- this reference carries the SAME
(-1)^3 sign as the kappa = -INT Dsachs_1 fold (both negative).  We report
|magnitude| and the sign of each explicitly.

UNITS: PCAMBz0.txt is h-units (k [h/Mpc], P [(Mpc/h)^3]).  We work entirely in
h-units internally (chi in Mpc/h, k in h/Mpc), then convert the final 3PCF to
physical units by h^4 (the equal-shell zeta carries h^4) -- matching canoes
units="physical".  2026-06-09 dimensional correction: the equal-shell
3-cumulant is natively (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, so units=physical applies
h^4 (was h^6; the surplus h^2 broke h-invariance of the dimensionless 3PCF).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.special import jv
from scipy.integrate import cumulative_trapezoid

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.cosmo.background import C_KM_S  # noqa: E402

Z_SOURCE = 5.0
_N_LAMBDA = 40       # source-shell grid (z in (0, z_s])
_N_ELL = 220         # ell quadrature nodes (log-spaced)
_N_PHI = 96          # orientation-angle nodes
_ELL_MIN = 1e-3
_ELL_MAX = 1000.0    # match canoes ell_high_max (LOW+HIGH combined cap)


def _growth_canoes(z: np.ndarray, cosmo) -> np.ndarray:
    """Linear growth D(z) consistent with the canoes tables (M[0]/(1+z)).

    canoes builds D via the lightcone amplitude; we reproduce the same D(z) by
    calling build_lightcone so the reference growth matches LOW/HIGH exactly.
    """
    from canoes.sachs.kappa3 import build_lightcone, chi_of_lambda  # type: ignore
    # We have z directly; get chi(z) in h-units and build the lightcone.
    chi_h = np.array([chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h) for zi in z]) * cosmo.h
    lightcone, one_plus_z = build_lightcone(cosmo, chi_h)
    D = lightcone.M[0, :] / one_plus_z
    return np.asarray(D, dtype=np.float64)


def _f2_spt(k_a, k_b, mu):
    """Standard symmetric SPT F2 kernel (NOT the effective F2)."""
    safe_a = np.where(k_a > 0, k_a, 1.0)
    safe_b = np.where(k_b > 0, k_b, 1.0)
    ratio = safe_a / safe_b + safe_b / safe_a
    return 5.0 / 7.0 + 0.5 * mu * ratio + (2.0 / 7.0) * mu * mu


def _build_grid(cosmo):
    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, Z_SOURCE, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])  # physical Mpc
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(Z_SOURCE / _N_LAMBDA, Z_SOURCE, _N_LAMBDA)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    return z, chi_phys, lam_phys, float(chif[-1]), float(lam_f[-1])


def _zeta_TTT_squeezed(gamma_rad, chi_phys_h, z, D, pk, cosmo):
    """Hand-rolled spin-0 squeezed zeta_TTT(gamma, chi) in h-units.

    Returns array (n_gamma, n_chi).  Internal k,chi are in h-units (the P(k)
    table is h-units); the h^4 physical-unit conversion is applied by the caller.
    """
    H0_inv_Mpc = 100.0 / C_KM_S  # H0 in h/Mpc (h-units): H0 = 100 h km/s/Mpc -> /c
    # In h-units H0 = 100 * h / c[km/s] * (1/h) = 100/c per (h/Mpc).  canoes uses
    # _H0_H_PER_MPC = 100/c (h-units H0).  Poisson amp A = -1.5 Om H0^2 (1+z).
    # log-spaced ell quadrature (trapezoid in log ell).
    ell = np.geomspace(_ELL_MIN, _ELL_MAX, _N_ELL)
    phi = np.linspace(0.0, 2.0 * np.pi, _N_PHI, endpoint=False)
    dphi = 2.0 * np.pi / _N_PHI
    u = ell[:, None, None]          # pairs with k2
    v = ell[None, :, None]          # pairs with k3 (carries gamma)
    ph = phi[None, None, :]
    cph = np.cos(ph)
    w = np.sqrt(np.maximum(u * u + v * v + 2.0 * u * v * cph, 0.0))  # third side k1
    # cos angles between the Fourier legs (matter triangle k1+k2+k3=0):
    safe_w = np.maximum(w, _ELL_MIN)
    cos_12 = -(u + v * cph) / safe_w   # angle(k1,k2)
    cos_13 = -(v + u * cph) / safe_w   # angle(k1,k3)
    cos_23 = cph                       # angle(k2,k3)
    # squeezed spin-0 angular kernel: 2 pi J0(v * gamma) per (this is the
    # orientation integral for spins=(0,0,0) with two coincident points).
    one_plus_z = 1.0 + z
    A_poisson = -1.5 * cosmo.Omega_m * (H0_inv_Mpc ** 2) * one_plus_z  # (n_z,)
    prefactor = 1.0 / ((2.0 * np.pi) ** 4)

    n_g = gamma_rad.size
    n_z = z.size
    out = np.zeros((n_g, n_z), dtype=np.float64)
    # ell-grid measures: d^2 ell_u d^2 ell_v / (2pi)^4 -> u du dphi_u * v dv dphi_v;
    # the absolute orientation (phi_u) integral is done analytically into the
    # 2 pi J0 kernel, leaving u v du dv dphi (relative angle phi).
    # log-trapezoid weights: int f dell = int f ell dln(ell).
    dlnell = np.gradient(np.log(ell))
    w_u = (ell * dlnell)[:, None, None]
    w_v = (ell * dlnell)[None, :, None]
    base_meas = w_u * w_v * dphi * u * v   # u*v from the radial Jacobian d^2ell

    for iz in range(n_z):
        chi = float(chi_phys_h[iz])
        k1 = safe_w / chi
        k2 = np.maximum(u, _ELL_MIN) / chi
        k3 = np.maximum(v, _ELL_MIN) / chi
        p1 = np.asarray(pk(k1), dtype=np.float64)
        p2 = np.asarray(pk(k2), dtype=np.float64)
        p3 = np.asarray(pk(k3), dtype=np.float64)
        f12 = _f2_spt(k1, k2, cos_12)
        f13 = _f2_spt(k1, k3, cos_13)
        f23 = _f2_spt(k2, k3, cos_23)
        b_delta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
        # per-leg field response: Phi00 = A(a) (1+z)^4 delta (spin-0).
        response = (A_poisson[iz] * one_plus_z[iz] ** 4) ** 3
        radial = (prefactor * (D[iz] ** 4) / chi ** 4) * response
        # angular kernel depends on gamma via 2 pi J0(v gamma):
        for ig in range(n_g):
            ang = 2.0 * np.pi * jv(0, v * gamma_rad[ig])
            integrand = base_meas * b_delta * ang
            out[ig, iz] = radial * float(np.sum(integrand))
    return out


def main() -> int:
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    h = cosmo.h

    z, chi_phys, lam_phys, chi_s, lam_s = _build_grid(cosmo)
    chi_phys_h = chi_phys * h  # h-units for the internal k = ell/chi
    chi_s_h = chi_s * h
    D = _growth_canoes(z, cosmo)
    print(f"[ref] z_s={Z_SOURCE}, chi_s={chi_s:.2f} Mpc ({chi_s_h:.1f} Mpc/h), "
          f"lambda_s={lam_s:.2f} Mpc")
    print(f"[ref] growth D(z): D(0.125)={D[0]:.4f} D(z_s)={D[-1]:.4f}")

    gamma_arcmin = np.unique(np.concatenate(
        [np.array([1.0]), np.geomspace(0.5, 300.0, 9)]))
    gamma_rad = np.radians(gamma_arcmin / 60.0)
    print(f"[ref] gamma grid (arcmin): {np.array2string(gamma_arcmin, precision=3)}")
    print(f"[ref] ell quad: n_ell={_N_ELL} [{_ELL_MIN:g},{_ELL_MAX:g}], n_phi={_N_PHI}")

    # zeta in h-units, convert to physical (units='physical' applies h^4).
    # 2026-06-09 dimensional correction: equal-shell zeta is natively (h/Mpc)^4
    # = A(a)^3 (P.P)/chi^4, so the physical-unit power is h^4 (was h^6).
    zeta_h = _zeta_TTT_squeezed(gamma_rad, chi_phys_h, z, D, pk, cosmo)
    zeta_TTT = zeta_h * (h ** 4)

    # --- fold (chi route, identical to ours-side chi route) --------------------
    a = 1.0 / (1.0 + z)
    K_chi = chi_phys * (chi_s - chi_phys) / chi_s  # physical Mpc
    integrand = (K_chi ** 3 * a ** 8)[None, :] * zeta_TTT
    Z_ref = np.trapezoid(integrand, chi_phys, axis=1)

    out = _HERE / "outputs" / "stage0_spt_reference.npz"
    np.savez(
        out, gamma_arcmin=gamma_arcmin, gamma_rad=gamma_rad,
        z=z, chi_phys=chi_phys, lam_phys=lam_phys, D=D,
        zeta_TTT=zeta_TTT, Z_ref=Z_ref, chi_s=chi_s, lam_s=lam_s,
        Omega_m=cosmo.Omega_m, h=cosmo.h, n_s=cosmo.n_s,
    )
    print(f"[ref] saved -> {out}")

    i1 = int(np.argmin(np.abs(gamma_arcmin - 1.0)))
    print()
    print(f"[ref] === SINGLE-GAMMA (gamma={gamma_arcmin[i1]:.3f}') ===")
    print(f"[ref]   Z_ref = {Z_ref[i1]:+.6e}   sign={'+' if Z_ref[i1] > 0 else '-'}")
    print()
    print("[ref] === SWEEP ===")
    for ig in range(gamma_arcmin.size):
        print(f"[ref]   gamma={gamma_arcmin[ig]:8.3f}'  Z_ref={Z_ref[ig]:+.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
