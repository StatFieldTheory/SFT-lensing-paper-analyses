#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""DECISIVE (A)-locus test: does the canoes 2-pt C_l^{Phi00Phi00}(chi,chi) diagonal
carry the per-leg (1+z)^4 photon-energy factor, or has it already reduced to the
standard (1+z)^1 lensing window per leg?

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

If the 2-pt per-shell harmonic object C_l^{Phi00}(chi,chi) scales (vs the standard
Limber C_l^{kappa} integrand W^2/chi^2 P_m) as (1+z)^0 per shell, then the 2-pt
REDUCES the per-leg (1+z)^4 internally (in the harmonic projection), and the 3-pt
equal-shell zeta SHOULD do the same -> the surviving (1+z)^4 in the 3-pt is a BUG.
If instead the 2-pt per-shell carries (1+z)^2 per leg (uncompensated) and only the
fold compensates, the picture is different.

We compute, at fixed ell, the canoes 2-pt diagonal Cl_pp(chi,chi) over a z-grid and
compare its z-scaling to the standard single-plane Limber convergence integrand
   I_std(chi) = [g(chi)]^2 / chi^2 * P_m((ell+0.5)/chi) * D(z)^2,
   g = (3/2)Om H0^2 (1+z) chi(chi_s-chi)/chi_s,
both divided by P_m * D^2 (clustering) to isolate the per-shell RESPONSE z-scaling.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_HERE = Path(__file__).resolve().parent
_OUT = _HERE.parent / "outputs"
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.cosmo.pk import _CambTablePk  # noqa: E402
from canoes.sachs import compute_kappa2_grid  # noqa: E402

H0 = 100.0 / 299792.458


def _load_pk_phi(cosmo):
    """P_Phi(k) = (A(a=1)/k^2)^2 P_delta? No -- canoes kappa2 wants P_Phi at z=0.
    Build P_Phi from P_delta via Poisson at a=1: Phi = A(1) delta / k^2,
    A(1) = -3/2 Om H0^2.  P_Phi(k) = (A(1)/k^2)^2 P_delta(k)."""
    arr = np.loadtxt("/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt")
    k_h, p_m = arr[:, 0], arr[:, 1]
    A1 = -1.5 * cosmo.Omega_m * H0 ** 2
    p_phi = (A1 / k_h ** 2) ** 2 * p_m
    return _CambTablePk(k_h.copy(), p_phi.copy(), n_s=float(cosmo.n_s),
                        extrapolation_low_k=True), (k_h, p_m)


def main() -> int:
    cosmo = FiducialCosmology()
    h = cosmo.h
    pk_phi, (k_h, p_m) = _load_pk_phi(cosmo)

    # z-grid of source planes (diagonal: chi_s = chi shell, single-plane response)
    z = np.linspace(0.4, 4.5, 10)
    chi_h = np.array([chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h) for zi in z]) * h
    ell = np.array([200], dtype=np.intp)

    # canoes 2-pt DIAGONAL Cl_pp(chi_i, chi_i) at each shell, single ell.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = compute_kappa2_grid(
            ell, chi_h, pk_phi_today=pk_phi, cosmo=cosmo, units="physical",
            radial_measure="lambda", lambda_convention="project",
        )
    cl_pp = np.asarray(out.Cl_phi00_phi00)  # (n_ell, n_chi, n_chi)
    diag = np.array([cl_pp[0, i, i] for i in range(z.size)])  # Cl_pp(chi_i,chi_i)

    # divide out the clustering P_m(k)*D^2 to isolate the per-shell RESPONSE scaling.
    from canoes.sachs.kappa3 import build_lightcone, chi_of_lambda  # noqa
    lc, opz = build_lightcone(cosmo, chi_h)
    D = np.asarray(lc.M[0] / opz)
    chi_phys = chi_h / h
    kP = (ell[0] + 0.5) / chi_phys
    Pm = np.interp(kP, k_h * h, p_m / h ** 3)  # physical P_m at k=(l+.5)/chi

    # canoes per-shell response (Cl_pp diag / [P_m D^2 / chi^2]) -- the (1+z) scaling:
    resp_canoes = diag / (Pm * D ** 2 / chi_phys ** 2)
    # standard single-plane convergence integrand response per shell (chi_s=chi -> 0;
    # instead use the per-leg window response g/(efficiency) = (1+z) scaling):
    # standard per-leg response ~ (1+z)^1 (the W window's photon-energy power).
    opz_arr = 1.0 + z
    # normalise both to shell 0 and read the z-power:
    rc = resp_canoes / resp_canoes[0]
    p_canoes = np.polyfit(np.log(opz_arr), np.log(np.abs(rc) + 1e-300), 1)[0]

    print("=" * 78)
    print("2-pt per-shell RESPONSE z-scaling (canoes Cl_pp diag / clustering)")
    print("=" * 78)
    print(f"{'z':>7s} {'1+z':>7s} {'Cl_pp_diag':>13s} {'resp/clustering':>15s} "
          f"{'(1+z)^4':>9s} {'(1+z)^2':>9s}")
    for i in range(z.size):
        print(f"{z[i]:7.3f} {opz_arr[i]:7.3f} {diag[i]:13.4e} {resp_canoes[i]:15.4e} "
              f"{opz_arr[i]**4:9.3f} {opz_arr[i]**2:9.3f}")
    print(f"\n  canoes 2-pt per-shell response z-power (fit) = {p_canoes:.3f}")
    print(f"  [(1+z)^4 = pure Sachs response; (1+z)^? = standard window photon-energy]")
    print(f"  If ~4: the 2-pt per-shell Cl_pp KEEPS the full (1+z)^4 (reduction is in")
    print(f"         the FOLD, so the 3-pt should also be compensated by ITS fold).")
    print(f"  If ~1: the 2-pt Cl_pp already reduced to the standard window -> 3-pt zeta")
    print(f"         should too -> the 3-pt's (1+z)^4 surplus is a builder bug.")
    np.savez(_OUT / "twopt_perleg_reduction.npz", z=z, diag=diag,
             resp_canoes=resp_canoes, p_canoes=p_canoes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
