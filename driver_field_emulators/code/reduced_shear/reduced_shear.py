"""The reduced-shear correction, as a yardstick for the FK amplitude.

WHY THIS COMPARISON. The FK term is a squeezed matter bispectrum entering
the two-point function because the observable is a nonlinear functional of
the field. Weak lensing already has an effect with exactly that structure,
and it is not the post-Born correction: post-Born is the perturbed path,
which is the FF term's analogue. The structural analogue of FK is the
REDUCED-SHEAR correction. What a survey measures is

    g = gamma / (1 - kappa) ~ gamma (1 + kappa),

so the observed two-point function picks up a term ~ 2 <gamma gamma kappa>,
a matter bispectrum in a squeezed configuration entering a two-point
statistic through the nonlinearity of the observable. It is a well studied
effect (Dodelson, Shapiro and White 2006; Shapiro 2009; Krause and Hirata
2010) and it is known to be of order a percent at the multipoles where
cosmic shear is measured. If the corrected FK amplitude is comparable, that
is a positive result: the same order of magnitude reached from a different
starting point.

WHAT IS COMPUTED. In the flat sky, writing the shear as
gamma(l) = kappa(l) exp(2 i phi_l), the leading correction to the E-mode
power spectrum is

    dC_l = 2 INT d^2 L / (2 pi)^2  cos(2 (phi_l - phi_L))
                 B_kappa(L, |l - L|, l),

with the convergence bispectrum taken in the Limber form that
`code/normalisation/verify_chain.py` showed the pipeline's own factors
reduce to exactly,

    B_kappa(l1,l2,l3) = INT d chi  W(chi)^3 / chi^4  B_delta(l1/chi, ..., z),
    W = (3/2) Omega_m H0^2 (1+z) chi (chi_s - chi) / chi_s .

CONVENTION WARNING. The prefactor and the cos(2 delta phi) projection are
written as they are standardly quoted, but they have NOT been checked
against the source papers here, which are not in the local library. An
error of a factor of two in the prefactor would not change the conclusion
of an order-of-magnitude comparison, and would change a quoted ratio. Treat
the ratio as good to a factor of order unity until the convention is
verified against Krause and Hirata.

Usage::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:.. $CAN/.venv/bin/python reduced_shear.py --z-s 5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from b_model import FIDUCIAL, PK_TABLE_3PT, make_b_delta_fn  # noqa: E402
from b_model.halofit_backend import _halofit_instance  # noqa: E402
from canoes.cosmo.lightcone import build_lightcone  # noqa: E402


def chi_of_z(cosmo, z_target: float) -> float:
    """Comoving distance in Mpc/h, on the same background the pipeline uses."""
    grid = np.linspace(1.0, 20000.0, 40000)
    _lc, opz = build_lightcone(cosmo, grid, distance_unit="Mpc/h")
    return float(np.interp(z_target, opz - 1.0, grid))

#: H0 in h/Mpc, so that chi in Mpc/h and H0 chi are both h-free.
H0_H_PER_MPC = 1.0 / 2997.92458


def lensing_weight(chi, chi_s, opz, omega_m):
    """W(chi) = (3/2) Omega_m H0^2 (1+z) chi (chi_s - chi) / chi_s."""
    return (1.5 * omega_m * H0_H_PER_MPC**2 * opz
            * chi * np.maximum(chi_s - chi, 0.0) / chi_s)


def radial_grid(cosmo, z_s, n_chi, ell_top, k_max):
    """Comoving grid, cut inside where the Limber map would leave the table."""
    chi_s = chi_of_z(cosmo, z_s)
    chi_min = max(2.0, 1.02 * ell_top / k_max)
    chi = np.exp(np.linspace(np.log(chi_min), np.log(chi_s * (1 - 1e-6)), n_chi))
    _lc, opz = build_lightcone(cosmo, chi, distance_unit="Mpc/h")
    return chi, opz - 1.0, opz, chi_s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--z-s", type=float, default=5.0)
    ap.add_argument("--ells", nargs="+", type=float,
                    default=[100, 200, 400, 800, 1600, 3000])
    ap.add_argument("--n-chi", type=int, default=160)
    ap.add_argument("--n-L", type=int, default=48)
    ap.add_argument("--n-phi", type=int, default=32)
    ap.add_argument("--L-min", type=float, default=2.0)
    ap.add_argument("--L-max", type=float, default=2.0e4)
    ap.add_argument("--k-max", type=float, default=200.0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    cosmo = FIDUCIAL
    b_delta = make_b_delta_fn("bihalofit")
    hf = _halofit_instance(cosmo, PK_TABLE_3PT)

    ell_top = 2.0 * args.L_max
    chi, z, opz, chi_s = radial_grid(cosmo, args.z_s, args.n_chi, ell_top, args.k_max)
    W = lensing_weight(chi, chi_s, opz, cosmo.Omega_m)
    dchi = np.gradient(chi)

    def C_kappa(ell):
        """Limber convergence power spectrum with the nonlinear P(k)."""
        k = ell / chi
        ok = k < args.k_max
        p = np.zeros_like(chi)
        if ok.any():
            p[ok] = np.asarray(
                hf.get_pkhalofit(k[ok], np.array([0.0]))).ravel()[:ok.sum()] * 0.0
            # halofit wants z per k; evaluate shell by shell
            p[ok] = np.array([float(np.ravel(hf.get_pkhalofit(
                np.array([kk]), np.array([zz])))[0]) for kk, zz in zip(k[ok], z[ok])])
        return float(np.sum(dchi * W**2 / chi**2 * p))

    # polar grid for the d^2 L integral
    Lg = np.exp(np.linspace(np.log(args.L_min), np.log(args.L_max), args.n_L))
    phg = np.linspace(0.0, 2 * np.pi, args.n_phi, endpoint=False)
    dlnL = np.gradient(np.log(Lg))
    dph = 2 * np.pi / args.n_phi

    print(f"z_s = {args.z_s},  chi_s = {chi_s:.1f} Mpc/h,  "
          f"{args.n_chi} shells, L in [{args.L_min:g}, {args.L_max:g}]\n")
    print(f"{'ell':>7} {'C_l^kk':>13} {'dC_l reduced':>15} {'dC/C':>10}")
    rows = []
    for ell in args.ells:
        Lx = np.outer(Lg, np.cos(phg))
        Ly = np.outer(Lg, np.sin(phg))
        l3 = np.hypot(ell - Lx, -Ly)                    # |l - L|
        l2 = np.broadcast_to(Lg[:, None], l3.shape)     # |L|
        cos2 = np.cos(2 * np.arctan2(Ly, Lx))           # phi_l = 0 by choice
        # B_kappa for every (L, phi) pair, accumulated over shells
        acc = np.zeros_like(l3)
        for j in range(chi.size):
            k1, k2, k3 = ell / chi[j], l2 / chi[j], l3 / chi[j]
            good = (k2 < args.k_max) & (k3 < args.k_max) & (k1 < args.k_max)
            if not good.any():
                continue
            b = np.zeros_like(l3)
            b[good] = b_delta(np.full(good.sum(), k1), k2[good], k3[good], z[j])
            acc += dchi[j] * W[j] ** 3 / chi[j] ** 4 * b
        integ = np.sum(dlnL[:, None] * Lg[:, None] ** 2 * dph * cos2 * acc)
        dC = 2.0 * integ / (2 * np.pi) ** 2
        C = C_kappa(ell)
        print(f"{ell:7.0f} {C:13.4e} {dC:15.4e} {100*dC/C:9.3f}%")
        rows.append((ell, C, dC, dC / C))

    if args.out:
        np.savez(args.out, ell=np.array([r[0] for r in rows]),
                 C_kappa=np.array([r[1] for r in rows]),
                 dC_reduced=np.array([r[2] for r in rows]),
                 ratio=np.array([r[3] for r in rows]), z_s=args.z_s)
        print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
