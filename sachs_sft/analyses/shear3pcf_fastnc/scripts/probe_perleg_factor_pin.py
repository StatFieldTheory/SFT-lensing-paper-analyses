#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Pin the EXACT per-shell radial-factor discrepancy canoes-W_can vs B_kappa.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

The decomposition probe showed the per-shell canoes/ref radial ratio is
L-INDEPENDENT (=52.19 at chi idx 20).  Here we decompose that single number into
the ratio of EACH factor between the two convergence-3pt weight assemblies, to
identify which factor(s) make canoes' W_can != B_kappa.  Both are scalar (spin-0)
convergence weights; whatever the discrepancy is, it is the SAME for TTT/PPP/D
(spin-independent) and is therefore NOT the spin-2 4-17x by itself -- but pinning
it tells us whether the B_kappa reference and the canoes/STAGE-0 reference are
two DIFFERENT normalisations (in which case the spin-2 verdict is confounded) or
the same.

canoes per-shell folded weight (after b_delta cancels), h-units:
    S_can(chi) = K_g^3 * a^2 * growth^4 / chi^4 * [A(a)(1+z)^4]^3 * radial_factor * h^6
    K_g = a^2 chi (chi_s - chi)/chi_s,   A(a) = -1.5 Om H0^2 (1+z),
    (the a^2 is the dlambda = a^2 dchi Jacobian since the fold is INT dlambda K_g^3 zeta)

B_kappa per-shell weight, h-units:
    S_ref(chi) = g(chi)^3 * (1/chi) * (1+z)^3 * growth^4
    g(chi) = 1.5 H0^2 Om (1 - chi/chi_s)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import H0_HUNIT, OMEGA_M  # noqa: E402

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    chi_of_lambda, build_lightcone, _lambda_jacobian,
    _kappa3_radial_density_conversion, _H0_H_PER_MPC,
)
from scipy.integrate import cumulative_trapezoid  # noqa: E402


def main():
    cosmo = FiducialCosmology()
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5.0 / 40, 5.0, 40)
    chi_phys = np.interp(z, zf, chif)
    lam_h = np.interp(z, zf, lam_f) * cosmo.h
    chi_s = float(chif[-1])
    chi_s_h = chi_s * cosmo.h
    a = 1.0 / (1.0 + z)
    opz = 1.0 + z

    chi_arr = np.maximum(chi_of_lambda(cosmo, lam_h, convention="project"), 1e-6)  # Mpc/h
    lc, opz_c = build_lightcone(cosmo, chi_arr)
    growth = np.asarray(lc.M[0] / opz_c, dtype=np.float64)
    djdc = _lambda_jacobian(opz_c, convention="project")
    rf, _, _ = _kappa3_radial_density_conversion(djdc, "lambda")

    i = 20
    print(f"shell {i}: z={z[i]:.4f}, chi_h={chi_arr[i]:.2f} Mpc/h, "
          f"1+z={opz[i]:.4f}, a={a[i]:.4f}")
    print(f"  H0_HUNIT(ref)={H0_HUNIT:.6e}  _H0_H_PER_MPC(canoes)={_H0_H_PER_MPC:.6e}")
    print(f"  ratio H0^2: {(_H0_H_PER_MPC/H0_HUNIT)**2:.6f}")

    # canoes per-leg lensing efficiency embedded in the fold:
    #   per-leg  =  K_g * A(a)*(1+z)^4    (the cube is the 3-leg product)
    #   plus per-shell scalar carriers a^2 (lambda Jacobian), 1/chi^4, growth^4,
    #   radial_factor, h^6.
    Kg = (a ** 2) * chi_arr * (chi_s_h - chi_arr) / chi_s_h     # h-units (chi in Mpc/h)
    A_can = -1.5 * cosmo.Omega_m * (_H0_H_PER_MPC ** 2) * opz   # A(a)
    perleg_can = Kg * A_can * (opz ** 4)                        # canoes per-leg efficiency
    # reference per-leg lensing efficiency:
    g_ref = 1.5 * (H0_HUNIT ** 2) * OMEGA_M * (1.0 - chi_arr / chi_s_h)  # g(chi)
    perleg_ref = g_ref * (opz ** 1)                            # B_kappa per-leg: g * (1+z)^?? see below

    # full per-shell scalar weights (b_delta and the (u,v,phi) kernel cancel):
    h6 = cosmo.h ** 6
    S_can = (perleg_can ** 3) * (a[i] ** 2 if False else 1.0)  # placeholder; build elementwise below
    S_can = (Kg ** 3) * (a ** 2) * (growth ** 4) / (chi_arr ** 4) \
        * ((A_can * opz ** 4) ** 3) * rf * h6
    S_ref = (g_ref ** 3) * (1.0 / chi_arr) * (opz ** 3) * (growth ** 4)

    print("\nPer-shell scalar weight at shell 20:")
    print(f"  S_can = {S_can[i]:.6e}")
    print(f"  S_ref = {S_ref[i]:.6e}")
    print(f"  S_can/S_ref = {S_can[i]/S_ref[i]:.6e}")

    print("\nFactor-by-factor per-leg ratio (canoes / ref), shell 20:")
    # canoes per-leg efficiency vs ref per-leg efficiency g:
    print(f"  K_g (canoes)               = {Kg[i]:.6e}  [Mpc/h]")
    print(f"  A(a)(1+z)^4 (canoes)       = {(A_can*opz**4)[i]:.6e}")
    print(f"  K_g*A(a)(1+z)^4 (canoes)   = {perleg_can[i]:.6e}")
    print(f"  g(chi) (ref)               = {g_ref[i]:.6e}")
    print(f"  per-leg canoes/ref(g)      = {perleg_can[i]/g_ref[i]:.6e}")
    print(f"  (per-leg)^3                = {(perleg_can[i]/g_ref[i])**3:.6e}")
    # remaining scalar carriers:
    rem_can = (a[i] ** 2) * (growth[i] ** 4) / (chi_arr[i] ** 4) * rf[i] * h6
    rem_ref = (1.0 / chi_arr[i]) * (opz[i] ** 3) * (growth[i] ** 4)
    print(f"\n  canoes remaining carriers  a^2/chi^4*rf*h^6 = {rem_can:.6e}")
    print(f"  ref remaining carriers     (1/chi)(1+z)^3    = {rem_ref:.6e}")
    print(f"  remaining canoes/ref       = {rem_can/rem_ref:.6e}")
    print(f"\n  TOTAL (per-leg)^3 * remaining ratio = "
          f"{((perleg_can[i]/g_ref[i])**3)*(rem_can/rem_ref):.6e}")
    print(f"  (should equal S_can/S_ref = {S_can[i]/S_ref[i]:.6e})")

    print("\nradial_factor rf[20] =", rf[i], " (should be 1.0 for lambda measure)")
    print("lambda Jacobian a^2[20] =", a[i] ** 2)


if __name__ == "__main__":
    main()
