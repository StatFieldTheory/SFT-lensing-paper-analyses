"""Build the PyCCL **FKEM non-Limber** reference for kappa-kappa, xi_+/-,
kappa-gamma at z_s = 5.

This is the non-Limber counterpart of the archived ``build_pyccl_shear_reference.py``
(which used ``ccl.angular_cl`` with the default ``l_limber=-1``, i.e. pure
Limber for every ell). Here C_L^{kappa-kappa} is computed with the
Fang-Krause-Eifler-Mishra-Sharma (FKEM) non-Limber integrator for ell below a
threshold, falling back to Limber at high ell where Limber is exact. The
result is the genuinely non-Limber angular spectrum, the apples-to-apples
match for the SFT corr_op (full non-Limber) C propagator.

The four real-space observables use the SAME direct full-sky Legendre /
Wigner-d sums as the Limber build, so any difference between the two
references is the pure Limber-vs-non-Limber effect on C_L, not a transform
artifact:

    xi_kappa(gamma)        = sum_L (2L+1)/(4 pi) C_L P_L(cos gamma)
    xi_plus(gamma)         = sum_L (2L+1)/(4 pi) C_L d^L_{2, 2}(gamma)
    xi_minus(gamma)        = sum_L (2L+1)/(4 pi) C_L d^L_{2,-2}(gamma)
    xi_kappa_gamma(gamma)  = sum_L (2L+1)/(4 pi) C_L d^L_{0, 2}(gamma)

Why FKEM only below a threshold? The non-Limber correction is a low-ell
effect (~+12% at ell=2, decaying to <0.1% by ell~30 for this z_s=5 kernel).
Limber is accurate to sub-percent at high ell, so running FKEM up to ell=4500
would be pointless and slow. ``--l-limber`` sets the crossover; the build
prints the FKEM/Limber C_L ratio so the crossover can be checked to be smooth.

Output schema (``pyccl_xi_shear_reference_fkem.npz``) matches the Limber build
plus ``l_limber`` / ``non_limber_method`` provenance keys.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
# spin_rotation lives in the shared sachs_sft scripts/ dir.
_SHARED = _HERE.parent.parent / "scripts"
if str(_SHARED) not in sys.path:
    sys.path.insert(0, str(_SHARED))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z-s", type=float, default=5.0)
    parser.add_argument("--lmax", type=int, default=4500)
    parser.add_argument(
        "--l-limber",
        type=int,
        default=250,
        help="ell below this uses FKEM non-Limber; at/above uses Limber "
        "(Limber is exact at high ell). Default 250 is well past the "
        "<0.1%% non-Limber residual at ell~30 for z_s=5.",
    )
    parser.add_argument("--fkem-nchi", type=int, default=512)
    parser.add_argument("--gamma-min-arcmin", type=float, default=0.5)
    parser.add_argument("--gamma-max-arcmin", type=float, default=5000.0)
    parser.add_argument("--n-gamma", type=int, default=40)
    parser.add_argument(
        "--out",
        type=Path,
        default=_HERE / "pyccl_xi_shear_reference_fkem.npz",
    )
    args = parser.parse_args()

    import pyccl as ccl
    from scipy.special import eval_legendre

    from spin_rotation import _wigner_d_array

    cosmo = ccl.Cosmology(
        Omega_c=0.12029 / 0.6711**2,
        Omega_b=0.02207 / 0.6711**2,
        h=0.6711,
        n_s=0.97,
        sigma8=0.81,
        transfer_function="boltzmann_camb",
        matter_power_spectrum="linear",
    )
    print(
        f"[pyccl-fkem-ref] cosmology: Omega_m={cosmo['Omega_m']:.4f}, "
        f"h={cosmo['h']:.4f}"
    )

    tracer_kappa = ccl.CMBLensingTracer(cosmo, z_source=args.z_s)
    ells = np.arange(2, args.lmax + 1, dtype=np.intp)
    ells_f = ells.astype(float)

    print(
        f"[pyccl-fkem-ref] C_l^kk: FKEM non-Limber for ell<{args.l_limber}, "
        f"Limber for ell>={args.l_limber} (ell={int(ells[0])}..{int(ells[-1])}, "
        f"z_s={args.z_s})"
    )
    # FKEM non-Limber below the crossover, Limber above. p_of_k_a is linear
    # here (matter_power_spectrum='linear'), so p_of_k_a_lin defaults match.
    Cl_kk = np.asarray(
        ccl.angular_cl(
            cosmo,
            tracer_kappa,
            tracer_kappa,
            ells_f,
            l_limber=args.l_limber,
            non_limber_integration_method="FKEM",
            fkem_Nchi=args.fkem_nchi,
        ),
        dtype=np.float64,
    )

    # Provenance / smoothness check: pure-Limber C_L for the FKEM/Limber ratio.
    Cl_limb = np.asarray(
        ccl.angular_cl(cosmo, tracer_kappa, tracer_kappa, ells_f),
        dtype=np.float64,
    )
    print(f"[pyccl-fkem-ref] FKEM/Limber C_L ratio at low ell:")
    for L in (2, 3, 5, 10, 20, 30, 50, 100, args.l_limber, args.l_limber + 5):
        if L <= ells[-1]:
            j = L - 2
            print(f"    ell={L:5d}  FKEM/Limber={Cl_kk[j] / Cl_limb[j]:.4f}")

    gamma_arcmin = np.logspace(
        np.log10(args.gamma_min_arcmin),
        np.log10(args.gamma_max_arcmin),
        args.n_gamma,
    )
    gamma_rad = gamma_arcmin * (np.pi / 180.0 / 60.0)

    print(
        f"[pyccl-fkem-ref] direct full-sky sums for N_gamma={args.n_gamma} "
        f"(P_L + 3 Wigner-d per gamma)"
    )

    n_g = args.n_gamma
    xi_kappa = np.zeros(n_g, dtype=np.float64)
    xi_plus = np.zeros(n_g, dtype=np.float64)
    xi_minus = np.zeros(n_g, dtype=np.float64)
    xi_kappa_gamma = np.zeros(n_g, dtype=np.float64)

    weight = (2.0 * ells_f + 1.0) / (4.0 * np.pi)
    w_Cl = weight * Cl_kk

    for ig, g_rad in enumerate(gamma_rad):
        cos_g = float(np.cos(g_rad))
        P_L = eval_legendre(ells, cos_g)
        d22 = _wigner_d_array(ells, 2, 2, cos_g)
        d2m2 = _wigner_d_array(ells, 2, -2, cos_g)
        d02 = _wigner_d_array(ells, 0, 2, cos_g)
        xi_kappa[ig] = float(np.sum(w_Cl * P_L))
        xi_plus[ig] = float(np.sum(w_Cl * d22))
        xi_minus[ig] = float(np.sum(w_Cl * d2m2))
        xi_kappa_gamma[ig] = float(np.sum(w_Cl * d02))

    print(
        f"[pyccl-fkem-ref] ranges  "
        f"xi_kk:  [{xi_kappa.min():+.3e}, {xi_kappa.max():+.3e}]\n"
        f"                        "
        f"xi_+:   [{xi_plus.min():+.3e}, {xi_plus.max():+.3e}]\n"
        f"                        "
        f"xi_-:   [{xi_minus.min():+.3e}, {xi_minus.max():+.3e}]\n"
        f"                        "
        f"xi_kg+: [{xi_kappa_gamma.min():+.3e}, {xi_kappa_gamma.max():+.3e}]"
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        ells=ells,
        Cl_kappa=Cl_kk,
        gamma_arcmin=gamma_arcmin,
        xi_kappa=xi_kappa,
        xi_plus=xi_plus,
        xi_minus=xi_minus,
        xi_kappa_gamma=xi_kappa_gamma,
        z_s=args.z_s,
        Omega_m=float(cosmo["Omega_m"]),
        h=float(cosmo["h"]),
        sigma8=float(cosmo["sigma8"]),
        n_s=float(cosmo["n_s"]),
        l_limber=args.l_limber,
        non_limber_method="FKEM",
        convention="full_sky_wigner_d_legendre_sum",
    )
    print(
        f"[pyccl-fkem-ref] saved {args.out}  "
        f"({args.out.stat().st_size / 1e6:.2f} MB)"
    )


if __name__ == "__main__":
    main()
