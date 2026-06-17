#!/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
"""CHECK A: is the hand-rolled PyCCL tree-3pt reference TRUSTWORTHY?

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python

STRATEGY
--------
The "(1+z)^4 over-normalization" conclusion rests entirely on the hand-rolled
PyCCL tree-3pt reference (probe_scalar_vs_pyccl_tree.py) being correctly
normalized.  The project lead's standing warning: comparing a lambda-integral
to a chi-integral REQUIRES matching the Jacobian/measure, so the *reference
itself* may carry a measure error.

We test the reference's OWN radial machinery by reducing it to the 2-point.
We re-use the SAME ingredients as probe_scalar_vs_pyccl_tree.py:
  * SAME cosmology (Omega_m, h, n_s, w0).
  * SAME P(k) table PCAMBz0.txt, fed identically (k_phys=k_h*h, P_phys=P_h/h^3,
    z=0 table, growth D(a)^2 applied separately) -- the SAME pk_phys_z0 closure.
  * SAME distances: CCL comoving_radial_distance on the SAME z-grid.
  * SAME growth: CCL growth_factor (D(1)=1).
  * SAME single source plane z_s=5, SAME chi_s.
  * SAME textbook lensing weight g(chi) = (3/2)Om H0^2 (1/a) chi (chi_s-chi)/chi_s
    [H0 in 1/Mpc] -- byte-identical _standard_g.

The reference's 2-point analogue is the Limber projection of g:
    C_l^kk = INT dchi g(chi)^2 / chi^2 * D(z)^2 * P_delta(l/chi, z=0) ,
where P_delta(k,z) = D(z)^2 P_delta(k,z=0) and the g already carries the full
radial weight (this is EXACTLY the 2-pt analogue of the reference's 3-pt
Z_std = INT dchi g^3/chi^4 D^4 (bispectrum fold): drop one g, drop one 1/chi,
drop two D's, replace the bispectrum J0 fold by the power spectrum).

CCL-NATIVE ANCHOR
-----------------
pyccl.WeakLensingTracer (single narrow source plane z_s=5) + angular_cl with the
SAME P(k) (linear, from the SAME cosmology) gives the rock-solid 2-pt C_l^kk.
We also build xi_kappa(theta) on both sides for a real-space cross-check.

VERDICT
-------
ratio = C_l^kk(reference machinery) / C_l^kk(CCL-native).
  ~1-2%  -> reference W kernel / chi-measure / Limber are CCL-correct ->
            the 3-pt reference is TRUSTWORTHY -> the (1+z)^4 lives in canoes
            (zeta or the fold), NOT in the reference.
  (1+z)^2-ish off -> reference carries a measure bug -> (1+z)^4 (partly) spurious.

Output: outputs/checkA_reference_2pt_vs_ccl.npz
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyccl as ccl
from scipy.integrate import cumulative_trapezoid

_HERE = Path(__file__).resolve().parent
_ANALYSIS = _HERE.parent
_OUT = _ANALYSIS / "outputs"

# ---- matched cosmology (IDENTICAL to probe_scalar_vs_pyccl_tree.py) ----------
OMEGA_M = 0.3160919980475834
H = 0.6711
N_S = 0.97
OMEGA_B = 0.0494
OMEGA_C = OMEGA_M - OMEGA_B
C_KM_S = 299792.458
H0_INV_MPC = 100.0 * H / C_KM_S
_PK_TABLE = Path("/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt")

# source plane (IDENTICAL)
Z_S = 5.0
_N_LAMBDA = 40  # reference's shell count (we densify for Limber convergence too)


def _sigma8_table(k_h, P_h):
    R = 8.0
    x = k_h * R
    W = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
    integ = k_h ** 2 * P_h * W ** 2 / (2 * np.pi ** 2)
    return float(np.sqrt(np.trapezoid(integ, k_h)))


def _build_cosmo_and_pk():
    """IDENTICAL to probe_scalar_vs_pyccl_tree._build_cosmo_and_pk."""
    arr = np.loadtxt(_PK_TABLE)
    k_h, P_h = arr[:, 0], arr[:, 1]
    s8 = _sigma8_table(k_h, P_h)
    cosmo = ccl.Cosmology(
        Omega_c=OMEGA_C, Omega_b=OMEGA_B, h=H, n_s=N_S, sigma8=s8, w0=-1.0,
        transfer_function="boltzmann_camb", matter_power_spectrum="linear",
        T_CMB=2.7255, Neff=3.044, m_nu=0.0,
    )
    k_phys = k_h * H
    P_phys = P_h / H ** 3
    lnk = np.log(k_phys)
    lnP = np.log(P_phys)

    def pk_phys_z0(k):
        lk = np.log(np.maximum(k, 1e-12))
        lp = np.interp(lk, lnk, lnP)
        lo = lk < lnk[0]
        if np.any(lo):
            lp = np.where(lo, lnP[0] + N_S * (lk - lnk[0]), lp)
        hi = lk > lnk[-1]
        if np.any(hi):
            slope = (lnP[-1] - lnP[-2]) / (lnk[-1] - lnk[-2])
            lp = np.where(hi, lnP[-1] + slope * (lk - lnk[-1]), lp)
        return np.exp(lp)
    return cosmo, s8, pk_phys_z0, (k_h, P_h)


def _grid(cosmo, z_s, n_lambda):
    """IDENTICAL grid construction to probe._grid (z linear in (0,z_s])."""
    zf = np.linspace(0.0, z_s, 6000)
    af = 1.0 / (1.0 + zf)
    chif = ccl.comoving_radial_distance(cosmo, af)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_s / n_lambda, z_s, n_lambda)
    chi = np.interp(z, zf, chif)
    lam = np.interp(z, zf, lam_f)
    return z, chi, lam, float(chif[-1]), float(lam_f[-1])


def _standard_g(chi_phys, z, chi_s):
    """IDENTICAL to probe._standard_g: textbook single-plane conv weight g(chi)."""
    a = 1.0 / (1.0 + z)
    return (1.5 * OMEGA_M * H0_INV_MPC ** 2 * (1.0 / a)
            * chi_phys * (chi_s - chi_phys) / chi_s)


def _reference_clkk(cosmo, pk_z0, ells, z_s, n_lambda):
    """The reference's OWN radial machinery reduced to the 2-point.

    C_l^kk = INT dchi [g(chi)^2 / chi^2] D(z)^2 P_delta(l/chi, z=0),
    chi-integral with the SAME g, SAME distances, SAME D, SAME pk_z0 closure.
    Returns C_l and also the (z, chi, g, D) grid for provenance.
    """
    z, chi, lam, chi_s, lam_s = _grid(cosmo, z_s, n_lambda)
    a = 1.0 / (1.0 + z)
    D = ccl.growth_factor(cosmo, a)
    g = _standard_g(chi, z, chi_s)
    cl = np.zeros_like(ells, dtype=np.float64)
    for il, ell in enumerate(ells):
        k = (ell + 0.5) / chi  # Limber: k = (l+1/2)/chi
        pk = pk_z0(k) * D ** 2  # P_delta(k, z) = D^2 P_delta(k, z=0)
        integ = g ** 2 / chi ** 2 * pk
        cl[il] = np.trapezoid(integ, chi)
    return cl, dict(z=z, chi=chi, lam=lam, chi_s=chi_s, g=g, D=D)


def main() -> int:
    cosmo, s8, pk_z0, (k_h, P_h) = _build_cosmo_and_pk()
    print(f"[checkA] sigma8(table)={s8:.5f}, CCL sigma8={ccl.sigma8(cosmo):.5f}")
    print(f"[checkA] chi(z=5)={ccl.comoving_radial_distance(cosmo, 1/6.):.2f} Mpc")

    # ---- ell grid -----------------------------------------------------------
    ells = np.unique(np.geomspace(10, 5000, 60).astype(int)).astype(float)

    # ---- reference 2-pt: SAME machinery, with grid-density convergence check -
    cl_ref_40, _ = _reference_clkk(cosmo, pk_z0, ells, Z_S, 40)
    cl_ref_200, ginfo = _reference_clkk(cosmo, pk_z0, ells, Z_S, 200)
    cl_ref_600, _ = _reference_clkk(cosmo, pk_z0, ells, Z_S, 600)
    # use the densest as the reference value; report the n_lambda=40 drift
    cl_ref = cl_ref_600

    # ---- CCL-NATIVE 2-pt: WeakLensingTracer single source plane z_s=5 -------
    z_src = Z_S
    zg = np.linspace(z_src - 0.06, z_src + 0.06, 41)
    nz = np.exp(-0.5 * ((zg - z_src) / 0.012) ** 2)
    wl = ccl.WeakLensingTracer(cosmo, dndz=(zg, nz), has_shear=True)
    pk_lin = cosmo.get_linear_power()  # SAME cosmology's linear P(k)
    cl_ccl = ccl.angular_cl(cosmo, wl, wl, ells, p_of_k_a=pk_lin)

    ratio = cl_ref / np.maximum(cl_ccl, 1e-300)

    print()
    print("=" * 88)
    print("2-POINT C_l^kk: reference machinery (g, chi-measure, Limber) vs CCL-native, z_s=5")
    print("=" * 88)
    print(f"{'ell':>8s} {'Cl_ref(n=600)':>15s} {'Cl_ref(n=40)':>15s} "
          f"{'Cl_CCL':>15s} {'ref/CCL':>9s}")
    sel = np.unique(np.linspace(0, ells.size - 1, 20).astype(int))
    for i in sel:
        print(f"{ells[i]:8.0f} {cl_ref[i]:15.5e} {cl_ref_40[i]:15.5e} "
              f"{cl_ccl[i]:15.5e} {ratio[i]:9.4f}")

    # plateau region (l ~ 100..2000) median ratio = the clean verdict number
    band = (ells >= 100) & (ells <= 2000)
    ratio_band = ratio[band]
    print()
    print(f"[verdict band l in 100..2000] median ref/CCL = {np.median(ratio_band):.4f}, "
          f"mean = {np.mean(ratio_band):.4f}, "
          f"min..max = {ratio_band.min():.4f}..{ratio_band.max():.4f}")
    # n_lambda=40 vs n_lambda=600 self-consistency of the reference machinery
    drift = np.abs(cl_ref_40 - cl_ref_600) / np.maximum(np.abs(cl_ref_600), 1e-300)
    print(f"[reference grid drift n40 vs n600] median = {np.median(drift[band]):.4f} "
          f"(quantifies the reference's own n_lambda=40 radial discretization error)")

    # ---- real-space cross-check: xi_kappa(theta) ----------------------------
    theta_arcmin = np.unique(np.concatenate(
        [np.array([1.0]), np.geomspace(0.5, 300.0, 9)]))
    ells_xi = np.unique(np.geomspace(2, 20000, 400).astype(int)).astype(float)
    cl_ref_xi, _ = _reference_clkk(cosmo, pk_z0, ells_xi, Z_S, 600)
    cl_ccl_xi = ccl.angular_cl(cosmo, wl, wl, ells_xi, p_of_k_a=pk_lin)
    xi_ref = ccl.correlation(cosmo, ell=ells_xi, C_ell=cl_ref_xi,
                             theta=theta_arcmin / 60.0, type="NN", method="fftlog")
    xi_ccl = ccl.correlation(cosmo, ell=ells_xi, C_ell=cl_ccl_xi,
                             theta=theta_arcmin / 60.0, type="NN", method="fftlog")
    xi_ratio = xi_ref / np.where(np.abs(xi_ccl) > 0, xi_ccl, 1e-300)
    print()
    print("=" * 88)
    print("REAL-SPACE xi_kappa(theta): reference machinery vs CCL-native, z_s=5")
    print("=" * 88)
    print(f"{'theta':>8s} {'xi_ref':>15s} {'xi_CCL':>15s} {'ref/CCL':>9s}")
    for i in range(theta_arcmin.size):
        print(f"{theta_arcmin[i]:8.3f} {xi_ref[i]:15.5e} {xi_ccl[i]:15.5e} "
              f"{xi_ratio[i]:9.4f}")
    band_xi = (theta_arcmin >= 1.0) & (theta_arcmin <= 100.0)
    print()
    print(f"[xi verdict band 1'..100'] median ref/CCL = {np.median(xi_ratio[band_xi]):.4f}")

    np.savez(
        _OUT / "checkA_reference_2pt_vs_ccl.npz",
        ells=ells, cl_ref=cl_ref, cl_ref_n40=cl_ref_40, cl_ref_n200=cl_ref_200,
        cl_ccl=cl_ccl, ratio=ratio,
        ratio_band_median=float(np.median(ratio_band)),
        ratio_band_mean=float(np.mean(ratio_band)),
        ref_grid_drift_band_median=float(np.median(drift[band])),
        theta_arcmin=theta_arcmin, xi_ref=xi_ref, xi_ccl=xi_ccl, xi_ratio=xi_ratio,
        xi_ratio_band_median=float(np.median(xi_ratio[band_xi])),
        sigma8_table=s8, Omega_m=OMEGA_M, h=H, n_s=N_S, z_s=Z_S,
        note="Reference's OWN radial machinery (same g, same distances, same growth, "
             "same PCAMBz0 P(k), same single source plane z_s=5) reduced to the 2-point "
             "Limber C_l^kk, vs CCL-native WeakLensingTracer angular_cl. ratio~1 => "
             "reference W/chi-measure/Limber are CCL-correct => 3pt reference trustworthy.",
    )
    print(f"\nsaved -> {_OUT / 'checkA_reference_2pt_vs_ccl.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
