#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Pin the canoes-B vs reference-B_kappa ratio as a function of (l1,l2,l3).

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

We have established the 4-17x lives in canoes' HIGH B/radial integrand, not the
orientation kernel.  This probe extracts the APPLES-TO-APPLES per-triangle ratio

    R(l1,l2,l3) = [canoes convergence-3pt weight] / [reference B_kappa]

where BOTH are the scalar coefficient that multiplies the SAME orientation kernel
in the radial+angular fold:

  reference:  B_kappa(l1,l2,l3) (dimensionless convergence bispectrum; the
              reference folds  INT (u du)(v dv) dphi/(2pi)^4 B_kappa * kernel).
  canoes:     W_can(l1,l2,l3) = INT dlambda K(lambda)^3 * [growth^4 * B_delta(k_i)
              * resp_PPP / chi^4] * h^6 * radial_factor,  k_i = l_i / chi.
              (the deployed fold is  sum_{u,v,phi} measure/(2pi)^4 * W_can * kernel;
              the (2pi)^-4 + measure are SHARED with the reference, so R = W_can/B_kappa
              is exactly the per-triangle B/normalisation discrepancy.)

If R is ~constant in (l1,l2,l3) -> a pure normalisation offset (same for all spins,
cannot explain spin-2-only 4-17x; would also break scalar STAGE-0, which it does
not).  If R grows as a POWER of l (e.g. ~l^4 = (l^2)^2 for the two spin-2-equivalent
legs) -> a per-leg eth^2 / sqrt(l^2(l^2-2)) ell-weight is the lever, and it is
spin-2-specific because the scalar/PPP responses already carry the right ell-power.

NOTE: W_can uses the PPP (psi0^3) response, the SAME response the deployed D-channel
uses.  So R is the scalar-coefficient discrepancy the J_2 kernel then amplifies.
Loads only our own trusted reference npz inputs (no untrusted pickle).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import BkappaLimber  # noqa: E402
sys.path.insert(0, str(Path("/Users/zzhang/projects/canoes")))
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    chi_of_lambda, build_lightcone, _f2_kernel_numpy,
    _kappa3_limber_response_product_h, _lambda_jacobian,
    _kappa3_radial_density_conversion, _H0_H_PER_MPC,
)
from scipy.integrate import cumulative_trapezoid


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    bk = BkappaLimber(n_chi=200)

    # LOS lambda grid + convergence kernel, identical to stage1_ours_v2
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1 / (1 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5 / 40, 5, 40)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * cosmo.h
    chi_s = float(chif[-1])
    a = 1 / (1 + z)
    Kg = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s
    chi_arr = np.maximum(chi_of_lambda(cosmo, lam_h, convention="project"), 1e-6)
    lc, opz = build_lightcone(cosmo, chi_arr)
    poisson_amp = -1.5 * cosmo.Omega_m * (_H0_H_PER_MPC ** 2) * opz
    djdc = _lambda_jacobian(opz, convention="project")
    rf, _, _ = _kappa3_radial_density_conversion(djdc, "lambda")

    def W_can(l1, l2, l3):
        """canoes convergence-3pt weight; matches the deployed fold's scalar coeff.
        The deployed sum is over (u,v,phi) of measure * W_lambda * kernel, then
        trapz over lambda with K^3; here we collapse to the per-triangle scalar
        weight INT dlambda K^3 W_lambda (trapezoid), to compare to B_kappa which
        the reference multiplies by the SAME (u,v,phi) measure and kernel."""
        c12 = (l3**2 - l1**2 - l2**2) / (2*l1*l2)
        c13 = (l2**2 - l1**2 - l3**2) / (2*l1*l3)
        c23 = (l1**2 - l2**2 - l3**2) / (2*l2*l3)
        integ = np.zeros_like(lam_phys)
        for il, cv in enumerate(chi_arr):
            cf = float(cv); k1, k2, k3 = l1/cf, l2/cf, l3/cf
            p1, p2, p3 = float(pk(k1)), float(pk(k2)), float(pk(k3))
            f12 = float(_f2_kernel_numpy(np.array([k1]), np.array([k2]), np.array([c12]))[0])
            f13 = float(_f2_kernel_numpy(np.array([k1]), np.array([k3]), np.array([c13]))[0])
            f23 = float(_f2_kernel_numpy(np.array([k2]), np.array([k3]), np.array([c23]))[0])
            bdelta = 2*(f12*p1*p2 + f13*p1*p3 + f23*p2*p3)
            growth = float(lc.M[0, il] / opz[il])
            resp = _kappa3_limber_response_product_h("PPP", float(poisson_amp[il]), float(opz[il]))
            integ[il] = (Kg[il]**3) * (growth**4) * bdelta / (cf**4) * resp * rf[il]
        return float(np.trapezoid(integ, lam_phys)) * (cosmo.h ** 6)

    print("=" * 96)
    print("R(l1,l2,l3) = W_can / B_kappa   (apples-to-apples convergence-3pt weight)")
    print("  Constant R -> normalisation (would break scalar; it does not).")
    print("  R ~ l^4 -> per-spin-2-leg eth^2/(l^2) ell-weight is the lever.")
    print("=" * 96)
    print(f"{'(l1,l2,l3)':>20s} {'W_can':>13s} {'B_kappa':>13s} {'R':>12s} "
          f"{'R/R0':>10s} {'(l/100)^4':>10s}")
    # equilateral sweep to expose ell-scaling cleanly; l1=l2=l3=L
    Ls = [80., 100., 150., 200., 300., 500., 800.]
    R0 = None
    rows = []
    for L in Ls:
        l1 = l2 = l3 = L
        wc = W_can(l1, l2, l3); bkv = float(bk(l1, l2, l3))
        R = wc / bkv
        if R0 is None:
            R0 = R
        print(f"{f'({L:.0f},{L:.0f},{L:.0f})':>20s} {wc:13.4e} {bkv:13.4e} {R:12.4e} "
              f"{R/R0:10.4f} {(L/100.)**4:10.4f}")
        rows.append((L, wc, bkv, R))
    # log-log slope of R vs L
    Larr = np.array([r[0] for r in rows]); Rarr = np.array([r[3] for r in rows])
    slope = np.polyfit(np.log(Larr), np.log(np.abs(Rarr)), 1)[0]
    print(f"\n  log-log slope d ln|R| / d ln L = {slope:.3f}  "
          f"(0 = pure normalisation; +4 = (l^2)^2 two-spin-2-leg eth^2)")
    np.savez(OUT / "HIGH_B_ratio_ellscaling.npz",
             L=Larr, W_can=np.array([r[1] for r in rows]),
             B_kappa=np.array([r[2] for r in rows]), R=Rarr, slope=slope,
             note="Per-triangle canoes-W/ref-B_kappa ratio vs equilateral L; "
                  "slope diagnoses normalisation (0) vs per-spin-2-leg ell-weight.")
    print(f"\nsaved -> {OUT / 'HIGH_B_ratio_ellscaling.npz'}")


if __name__ == "__main__":
    main()
