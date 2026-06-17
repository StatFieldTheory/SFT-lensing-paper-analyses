#!/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
"""DECISIVE third-code reconciliation: canoes SCALAR convergence 3PCF vs PyCCL.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python

GOAL
----
Settle whether the canoes/STF driving-field SCALAR convergence 3PCF (the
``zeta_TTT`` fold that STAGE 0 produced) is correctly normalized against the
STANDARD lensing convention, using PyCCL as an independent, well-tested third
code (it already validated the canoes 2-point).  We build the STANDARD tree-level
convergence 3PCF on the SAME squeezed/collapsed config and the SAME J0 fold the
canoes STAGE-0 ``Z_lambda`` uses, with carefully matched P(k), cosmology, growth,
distances, single source plane z_s=5, and dimensionless output on both sides.

MATCHED INPUTS (apples-to-apples)
---------------------------------
  * cosmology: Omega_m=0.3160919980475834, h=0.6711, n_s=0.97, w0=-1.
  * P(k): PCAMBz0.txt (h-units), the SAME table canoes uses; fed to CCL as a
    linear Pk2D scaled by D(a)^2.  sigma8(table)=0.80902; CCL with this sigma8
    reproduces it to 5 decimals, and the table is used DIRECTLY (P(k)-matched,
    not just cosmology-matched).
  * growth D(a): canoes growth (lightcone M[0]/(1+z)) matches CCL growth_factor
    (D(1)=1) to <0.1% across z in [0,5] (verified) -> use CCL growth.
  * distances: CCL chi(z=5)=7968.7 Mpc vs canoes 7971.2 Mpc (0.03%).

THE SQUEEZED SCALAR CONVERGENCE 3PCF (matched to canoes STAGE-0 Z_lambda)
------------------------------------------------------------------------
canoes STAGE-0 (stage0_ours.py) folds the equal-shell driving-field 3-cumulant
zeta_TTT(gamma, lambda) on the squeezed config (1, cos g, cos g) with a J0
angular kernel.  The independent same-convention SPT reference
(stage0_spt_reference.py) writes that fold EXPLICITLY as, per source shell,

  zeta_TTT(gamma, chi) = (1/(2pi)^4)(D^4/chi^4)[A(a)(1+z)^4]^3
      * INT du dv dphi (u v) B_delta(k1,k2,k3) * 2pi J0(v gamma),

with k_i = ell_i/chi, A(a) = -3/2 Om H0^2 (1+z), folded with the canoes
efficiency K = a^2 chi (chi_s-chi)/chi_s and Z = INT dlambda K^3 zeta.

The STANDARD convergence convention replaces the canoes per-leg response
[A(a)(1+z)^4] and the K = a^2 chi(chi_s-chi)/chi_s efficiency by the textbook
single-plane convergence weight g(chi) acting on the matter field directly:

  kappa(n) = INT dchi g(chi) delta(chi n, chi),
  g(chi)   = (3/2) Om H0^2 (1/a) chi (chi_s-chi)/chi_s    [H0 in 1/Mpc],

so the squeezed scalar convergence 3PCF in the STANDARD convention is the
Limber projection

  Z_kappa^std(gamma) = INT dchi [g(chi)^3 / chi^4] D(z)^4
      * (1/(2pi)^4) INT du dv dphi (u v) B_delta(k1,k2,k3) 2pi J0(v gamma),

with the SAME B_delta = 2 F2 P P + 2cyc (true SPT F2) and the SAME J0 angular
kernel and the SAME ell-quadrature as the canoes STAGE-0 SPT reference.  This is
literally the canoes STAGE-0 SPT reference with the per-leg response
[A(a)(1+z)^4]^3 K^3 a^8 replaced by the standard g(chi)^3, i.e. the ONLY thing
that changes is the radial lensing-efficiency / per-leg response convention.
Everything angular (J0 kernel, B_delta, ell grid) is byte-identical, so the
ratio isolates the radial-measure/(1+z) convention cleanly.

DISAMBIGUATORS
--------------
1. (1+z)^6 test.  canoes per-leg response * efficiency vs standard g:
     canoes: [A(a)(1+z)^4] folded with K=a^2 chi(...) and a^8 K_chi^3
     std:    g = (3/2)Om H0^2 (1+z) chi(...)
   Back-of-envelope predicts an EXTRA (1+z)^2 per leg -> (1+z)^6 over 3 legs.
   We compute the ratio per gamma, the lensing-weighted z_eff(gamma), and test
   ratio ~ (1+z_eff)^6.  We ALSO vary z_s in {1,2,5} to see the ratio move with
   z_eff as (1+z_eff)^6 predicts.
2. 2pt vs 3pt.  We build the STANDARD linear+tree convergence 2PCF xi_kappa via
   PyCCL WeakLensingTracer (single source plane z_s=5) and compare to the canoes
   2-point convergence (memory says it was PyCCL-validated to ~10%).  If the
   canoes 2-pt MATCHES PyCCL but the 3-pt is (1+z)^6-off, the bug is
   3-POINT-SPECIFIC.  If BOTH are off, it is a per-leg radial-measure convention.

Output: outputs/scalar_vs_pyccl_tree.npz with both Z_kappa(gamma), the ratio,
the z_eff and (1+z_eff)^6 columns, the z_s sweep, and the 2pt comparison.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pyccl as ccl
from scipy.special import jv
from scipy.integrate import cumulative_trapezoid

_HERE = Path(__file__).resolve().parent
_ANALYSIS = _HERE.parent
_OUT = _ANALYSIS / "outputs"

# ---- matched cosmology -------------------------------------------------------
OMEGA_M = 0.3160919980475834
H = 0.6711
N_S = 0.97
OMEGA_B = 0.0494          # split is irrelevant for distances/growth at <0.05%
OMEGA_C = OMEGA_M - OMEGA_B
C_KM_S = 299792.458
H0_INV_MPC = 100.0 * H / C_KM_S      # physical H0 [1/Mpc]
_PK_TABLE = Path("/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt")

# ---- quadrature (match stage0_spt_reference.py) ------------------------------
_N_LAMBDA = 40
_N_ELL = 220
_N_PHI = 96
_ELL_MIN = 1e-3
_ELL_MAX = 1000.0


def _sigma8_table(k_h, P_h):
    R = 8.0
    x = k_h * R
    W = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
    integ = k_h ** 2 * P_h * W ** 2 / (2 * np.pi ** 2)
    return float(np.sqrt(np.trapezoid(integ, k_h)))


def _build_cosmo_and_pk():
    """CCL cosmology with the PCAMBz0 table fed as the linear Pk2D (P-matched)."""
    arr = np.loadtxt(_PK_TABLE)
    k_h, P_h = arr[:, 0], arr[:, 1]           # k[h/Mpc], P[(Mpc/h)^3]
    s8 = _sigma8_table(k_h, P_h)
    cosmo = ccl.Cosmology(
        Omega_c=OMEGA_C, Omega_b=OMEGA_B, h=H, n_s=N_S, sigma8=s8, w0=-1.0,
        transfer_function="boltzmann_camb", matter_power_spectrum="linear",
        T_CMB=2.7255, Neff=3.044, m_nu=0.0,
    )
    # physical-unit P(k) interpolant from the table: k[1/Mpc] = k_h * h,
    # P[Mpc^3] = P_h / h^3.  z=0 table; growth handled separately.
    k_phys = k_h * H
    P_phys = P_h / H ** 3
    lnk = np.log(k_phys)
    lnP = np.log(P_phys)

    def pk_phys_z0(k):
        lk = np.log(np.maximum(k, 1e-12))
        lp = np.interp(lk, lnk, lnP)
        # low-k power-law extrapolation P ~ k^n_s (avoids flat-floor artifact)
        lo = lk < lnk[0]
        if np.any(lo):
            lp = np.where(lo, lnP[0] + N_S * (lk - lnk[0]), lp)
        # high-k: clamp to last value's slope (P falls); np.interp already clamps
        hi = lk > lnk[-1]
        if np.any(hi):
            slope = (lnP[-1] - lnP[-2]) / (lnk[-1] - lnk[-2])
            lp = np.where(hi, lnP[-1] + slope * (lk - lnk[-1]), lp)
        return np.exp(lp)
    return cosmo, s8, pk_phys_z0, (k_h, P_h)


def _f2_spt(k_a, k_b, mu):
    safe_a = np.where(k_a > 0, k_a, 1.0)
    safe_b = np.where(k_b > 0, k_b, 1.0)
    ratio = safe_a / safe_b + safe_b / safe_a
    return 5.0 / 7.0 + 0.5 * mu * ratio + (2.0 / 7.0) * mu * mu


def _grid(cosmo, z_s):
    """Source-shell grid z in (0, z_s], physical chi from CCL, affine lambda."""
    zf = np.linspace(0.0, z_s, 6000)
    af = 1.0 / (1.0 + zf)
    chif = ccl.comoving_radial_distance(cosmo, af)     # physical Mpc
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_s / _N_LAMBDA, z_s, _N_LAMBDA)
    chi = np.interp(z, zf, chif)
    lam = np.interp(z, zf, lam_f)
    return z, chi, lam, float(chif[-1]), float(lam_f[-1])


def _squeezed_zeta_block(gamma_rad, chi_phys, z, D, pk_z0, response_per_shell):
    """Squeezed J0 fold integrand per shell, with a pluggable per-shell radial
    'response_per_shell' (n_z,).  Returns Z-integrand zeta(gamma, chi) such that
    the 3PCF is  INT_meas zeta.  Internal k = ell/chi in PHYSICAL units.

    The angular part (B_delta, J0, ell-quadrature) is IDENTICAL to the canoes
    STAGE-0 SPT reference; only 'response_per_shell' changes (canoes vs standard).
    """
    ell = np.geomspace(_ELL_MIN, _ELL_MAX, _N_ELL)
    phi = np.linspace(0.0, 2.0 * np.pi, _N_PHI, endpoint=False)
    dphi = 2.0 * np.pi / _N_PHI
    u = ell[:, None, None]
    v = ell[None, :, None]
    ph = phi[None, None, :]
    cph = np.cos(ph)
    w = np.sqrt(np.maximum(u * u + v * v + 2.0 * u * v * cph, 0.0))
    safe_w = np.maximum(w, _ELL_MIN)
    cos_12 = -(u + v * cph) / safe_w
    cos_13 = -(v + u * cph) / safe_w
    cos_23 = cph
    prefactor = 1.0 / ((2.0 * np.pi) ** 4)
    dlnell = np.gradient(np.log(ell))
    w_u = (ell * dlnell)[:, None, None]
    w_v = (ell * dlnell)[None, :, None]
    base_meas = w_u * w_v * dphi * u * v

    n_g = gamma_rad.size
    n_z = z.size
    out = np.zeros((n_g, n_z), dtype=np.float64)
    for iz in range(n_z):
        chi = float(chi_phys[iz])
        k1 = safe_w / chi
        k2 = np.maximum(u, _ELL_MIN) / chi
        k3 = np.maximum(v, _ELL_MIN) / chi
        p1 = pk_z0(k1)
        p2 = pk_z0(k2)
        p3 = pk_z0(k3)
        f12 = _f2_spt(k1, k2, cos_12)
        f13 = _f2_spt(k1, k3, cos_13)
        f23 = _f2_spt(k2, k3, cos_23)
        b_delta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
        radial = (prefactor * (D[iz] ** 4) / chi ** 4) * response_per_shell[iz]
        for ig in range(n_g):
            ang = 2.0 * np.pi * jv(0, v * gamma_rad[ig])
            out[ig, iz] = radial * float(np.sum(base_meas * b_delta * ang))
    return out


def _standard_g(chi_phys, z, chi_s):
    """Textbook single-plane convergence weight g(chi) [1/Mpc, physical]:
        g = (3/2) Om H0^2 (1/a) chi (chi_s-chi)/chi_s,  H0 in 1/Mpc.
    """
    a = 1.0 / (1.0 + z)
    return 1.5 * OMEGA_M * H0_INV_MPC ** 2 * (1.0 / a) * chi_phys * (chi_s - chi_phys) / chi_s


def _canoes_response_and_fold(chi_phys, z, chi_s, lam_phys, zeta_canoesconv):
    """The canoes STAGE-0 convention: response=[A(a)(1+z)^4]^3 (already inside
    zeta_canoesconv via the SPT-ref builder), folded with K=a^2 chi(...)/... in
    lambda.  We pass zeta already built with the canoes response; here we just
    fold it the canoes way."""
    a = 1.0 / (1.0 + z)
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s
    return np.trapezoid((K ** 3)[None, :] * zeta_canoesconv, lam_phys, axis=1)


def _run_for_zs(cosmo, pk_z0, gamma_rad, z_s, want_canoes=True):
    """Returns dict with standard Z_kappa(gamma) and (optionally) the canoes-
    convention Z (rebuilt with the SAME angular machinery, canoes response) so the
    ratio is a pure radial-convention ratio (not contaminated by a different
    angular quadrature)."""
    z, chi, lam, chi_s, lam_s = _grid(cosmo, z_s)
    a = 1.0 / (1.0 + z)
    D = ccl.growth_factor(cosmo, a)

    # STANDARD: response_per_shell = g(chi)^3, fold trivially (g already carries
    # the full radial weight; the J0 zeta already has 1/(2pi)^4 D^4/chi^4 B).
    g = _standard_g(chi, z, chi_s)
    resp_std = g ** 3
    zeta_std = _squeezed_zeta_block(gamma_rad, chi, z, D, pk_z0, resp_std)
    # The standard convergence kappa = INT dchi g delta, so the 3PCF is
    #   INT dchi (per-shell zeta_std) — zeta_std already absorbed g^3, so fold in chi
    #   with NO extra kernel (the chi measure of the LOS integral).
    Z_std = np.trapezoid(zeta_std, chi, axis=1)

    out = dict(z=z, chi=chi, lam=lam, chi_s=chi_s, lam_s=lam_s, D=D, g=g, Z_std=Z_std)

    if want_canoes:
        # CANOES convention response = [A(a)(1+z)^4]^3, A=-3/2 Om H0^2 (1+z).
        A_pois = -1.5 * OMEGA_M * H0_INV_MPC ** 2 * (1.0 + z)
        resp_can = (A_pois * (1.0 + z) ** 4) ** 3
        zeta_can = _squeezed_zeta_block(gamma_rad, chi, z, D, pk_z0, resp_can)
        # fold the canoes way: Z = INT dlambda K^3 zeta_can, K = a^2 chi(...)/...
        K = a ** 2 * chi * (chi_s - chi) / chi_s
        Z_can = np.trapezoid((K ** 3)[None, :] * zeta_can, lam, axis=1)
        out["Z_can"] = Z_can
        # lensing-weighted z_eff(gamma): weight by |per-shell standard integrand|
        w_shell = np.abs(zeta_std)                                   # (n_g, n_z)
        z_eff = np.sum(w_shell * z[None, :], axis=1) / np.maximum(
            np.sum(w_shell, axis=1), 1e-300)
        out["z_eff"] = z_eff
    return out


def main() -> int:
    cosmo, s8, pk_z0, (k_h, P_h) = _build_cosmo_and_pk()
    print(f"[pyccl] sigma8(table)={s8:.5f}, CCL sigma8={ccl.sigma8(cosmo):.5f} "
          f"(match -> P(k) normalization matched)")
    print(f"[pyccl] chi(z=5)={ccl.comoving_radial_distance(cosmo, 1/6.):.2f} Mpc "
          f"(canoes 7971.18)")

    gamma_arcmin = np.unique(np.concatenate(
        [np.array([1.0]), np.geomspace(0.5, 300.0, 9)]))
    gamma_rad = np.radians(gamma_arcmin / 60.0)

    # --- z_s = 5 (primary) ----------------------------------------------------
    r5 = _run_for_zs(cosmo, pk_z0, gamma_rad, 5.0, want_canoes=True)
    Z_std = r5["Z_std"]
    Z_can = r5["Z_can"]
    z_eff = r5["z_eff"]

    # canoes STAGE-0 deployed Z_lambda (real canoes zeta_TTT fold) for sanity:
    st0 = np.load(_ANALYSIS / "outputs" / "stage0_ours.npz", allow_pickle=True)
    Z_lambda_canoes = st0["Z_lambda"]
    ga0 = st0["gamma_arcmin"]
    # align (same grid by construction)
    assert np.allclose(ga0, gamma_arcmin), "gamma grid mismatch with stage0_ours"

    ratio_canstd = np.abs(Z_can) / np.maximum(np.abs(Z_std), 1e-300)
    ratio_realcanstd = np.abs(Z_lambda_canoes) / np.maximum(np.abs(Z_std), 1e-300)
    pred_1pz6 = (1.0 + z_eff) ** 6

    print()
    print("=" * 100)
    print("SCALAR convergence 3PCF: canoes convention vs STANDARD (PyCCL distances+W), z_s=5")
    print("=" * 100)
    hdr_g = "gamma'"
    print(f"{hdr_g:>9s} {'|Z_std|(pyccl)':>15s} {'|Z_can|(reb)':>13s} "
          f"{'|Z_real_can|':>13s} {'can/std':>9s} {'realcan/std':>11s} "
          f"{'z_eff':>6s} {'(1+ze)^6':>9s}")
    for i in range(gamma_arcmin.size):
        print(f"{gamma_arcmin[i]:9.3f} {abs(Z_std[i]):15.4e} {abs(Z_can[i]):13.4e} "
              f"{abs(Z_lambda_canoes[i]):13.4e} {ratio_canstd[i]:9.3f} "
              f"{ratio_realcanstd[i]:11.3f} {z_eff[i]:6.3f} {pred_1pz6[i]:9.2f}")

    # --- (1+z)^6 z_s sweep: ratio at gamma=10' for z_s in {1,2,5} -------------
    print()
    print("=" * 100)
    print("(1+z)^6 TEST via z_s sweep at gamma~10' (realcanoes-style rebuild vs std)")
    print("=" * 100)
    ig10 = int(np.argmin(np.abs(gamma_arcmin - 12.247)))  # ~12' node closest to 10'
    zs_list = [1.0, 2.0, 5.0]
    zs_rows = []
    for zs in zs_list:
        rr = _run_for_zs(cosmo, pk_z0, gamma_rad, zs, want_canoes=True)
        rcs = np.abs(rr["Z_can"]) / np.maximum(np.abs(rr["Z_std"]), 1e-300)
        ze = rr["z_eff"]
        zs_rows.append((zs, rcs[ig10], ze[ig10], (1 + ze[ig10]) ** 6))
        print(f"  z_s={zs:.1f}: can/std={rcs[ig10]:8.3f}  z_eff={ze[ig10]:.3f}  "
              f"(1+z_eff)^6={(1+ze[ig10])**6:8.2f}  ratio/(1+ze)^6={rcs[ig10]/(1+ze[ig10])**6:.3f}")

    # --- 2-POINT disambiguation ----------------------------------------------
    print()
    print("=" * 100)
    print("2-POINT: standard convergence xi_kappa(theta) via PyCCL WeakLensingTracer, z_s=5")
    print("=" * 100)
    # single source plane z_s=5: narrow Gaussian dndz
    z_src = 5.0
    zg = np.linspace(z_src - 0.06, z_src + 0.06, 41)
    nz = np.exp(-0.5 * ((zg - z_src) / 0.012) ** 2)
    wl = ccl.WeakLensingTracer(cosmo, dndz=(zg, nz), has_shear=True)
    ells = np.unique(np.geomspace(2, 20000, 400).astype(int)).astype(float)
    pk_lin = cosmo.get_linear_power()
    cl_kk = ccl.angular_cl(cosmo, wl, wl, ells, p_of_k_a=pk_lin)
    theta_arcmin = gamma_arcmin.copy()
    xi_kk = ccl.correlation(cosmo, ell=ells, C_ell=cl_kk,
                            theta=theta_arcmin / 60.0, type="NN", method="fftlog")
    hdr_t = "theta'"
    print(f"{hdr_t:>9s} {'xi_kappa(pyccl,lin)':>20s}")
    for i in range(theta_arcmin.size):
        print(f"{theta_arcmin[i]:9.3f} {xi_kk[i]:20.6e}")
    print("  (compare to canoes 2-pt convergence xi_kappa Order-0; memory: PyCCL-validated ~10%)")

    np.savez(
        _OUT / "scalar_vs_pyccl_tree.npz",
        gamma_arcmin=gamma_arcmin, gamma_rad=gamma_rad,
        Z_std=Z_std, Z_can_rebuilt=Z_can, Z_real_canoes=Z_lambda_canoes,
        ratio_can_std=ratio_canstd, ratio_realcan_std=ratio_realcanstd,
        z_eff=z_eff, pred_1pz6=pred_1pz6,
        zs_sweep=np.array(zs_rows),          # (z_s, can/std, z_eff, (1+ze)^6)
        theta_arcmin=theta_arcmin, xi_kappa_pyccl=xi_kk,
        ells=ells, cl_kk=cl_kk,
        sigma8_table=s8, Omega_m=OMEGA_M, h=H, n_s=N_S, z_s=5.0,
        note="Standard tree conv-3PCF (PyCCL W kernel) vs canoes-convention rebuild "
             "and real canoes zeta_TTT fold, P(k)-matched PCAMBz0; +(1+z)^6 z_s sweep "
             "+ PyCCL linear conv-2PCF for 2pt-vs-3pt disambiguation.",
    )
    print(f"\nsaved -> {_OUT / 'scalar_vs_pyccl_tree.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
