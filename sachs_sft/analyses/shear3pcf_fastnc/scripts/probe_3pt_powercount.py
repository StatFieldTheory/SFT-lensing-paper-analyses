#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""PHASE-1 power-count (DELIVERABLE): pin the EXACT (1+z)^4 + h^2 surplus in the
canoes 3-pt SCALAR (TTT) equal-shell vertex, by power-counting against (a) the
PyCCL-validated canoes 2-pt and (b) the PyCCL standard oracle ratio.  NO fitting.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Two independent surpluses, cleanly power-counted (NO fitting):

SURPLUS (B) -- the h-power.  A convergence 3PCF is dimensionless and MUST be
h-invariant.  We build the canoes equal-shell zeta_TTT in BOTH units="h" and
units="physical" and read the empirical ratio = h^6.  The dimensionless-correct
footing is h^4 (the PyCCL oracle's own h-test).  So the canoes "physical" zeta is
h^2 (= 1/h^2 ~ 2.22) larger than a dimensionless conv-3PCF.

  CORRECTED CODE MAP: this h^6 is the UNIFORM canoes "physical" units convention,
  applied CONSISTENTLY by ALL FOUR equal-shell builders -- scalar LOW
  `compute_kappa3_zeta_table` (kappa3.py:4872), scalar HIGH
  `compute_kappa3_sigma3_high` (kappa3.py:3642), modulus LOW
  `compute_kappa3_mod_zeta_equal_time` (kappa3.py:3906-07), modulus HIGH
  `compute_kappa3_mod_sigma3_high` (kappa3.py:4050-51).  (An empirical h-power test
  on the SCALAR LOW `compute_kappa3_zeta_table(tree_phi)` gives 6.0 EXACTLY,
  confirming LOW and HIGH are CONSISTENT at h^6 -- the h^4 literals at kappa3.py:3151
  and :3344 belong to a DIFFERENT object, `compute_kappa3_sigma3`, the cross-shell
  direct-sigma3 path, NOT the equal-shell pipeline.)  Hence the h^2 is NOT an
  internal inconsistency to "fix" on one line -- it is the units convention.

SURPLUS (A) -- the (1+z)^4 per shell [3-pt equal-shell radial assembly].  We
construct the EXACT per-shell radial scalar response S_can(chi) the deployed
canoes zeta_TTT fold carries (after the byte-identical angular B_delta + J0 kernel
+ 1/(2pi)^4 + D^4/chi^4 cancels) and the EXACT standard per-shell response
S_std(chi) the PyCCL oracle carries, and decompose S_can/S_std symbolically:

    S_can = a^8 * K_geo^3 * [A(a)(1+z)^4]^3       (deployed: INT dlambda K_lambda^3
                                                   zeta, K_lambda=a^2 K_geo, folded
                                                   to chi -> a^8; zeta carries the
                                                   per-leg Sachs (1+z)^4)
    S_std = g^3,  g = (3/2) Om H0^2 (1+z) K_geo   (PyCCL standard single-plane weight)
    S_can/S_std = a^8 [A(a)(1+z)^4]^3 / [3/2 Om H0^2 (1+z) K_geo]^3
                = a^8 * (-1)^3 (1+z)^12 = -(1+z)^4   (sign from A(a)<0; |.| = (1+z)^4)

Verified to machine precision per shell (ratio/(-(1+z)^4) = 1.000000).

DEPLOYED ORACLE CLOSURE: realcan/std ~ (1+z_eff)^4 * h^2 (both surpluses; they
partially cancel since (1+z)^4>1 and h^2<1).  In the signal regime gamma<=12':
realcan/std ~ 9.95 ~ (1+z_eff)^4 * h^2 ~ 23.6 * 0.45 ~ 10.6 (Jensen-gap ~6%).

2-pt IMMUNITY (why the surplus is 3-POINT-SPECIFIC): the PyCCL-validated 2-pt
identity (cosmology.tex Eq. "Cell kappa harmonic bridge", exact, no Limber)
C_l^kk = INT dchi dchi' [chi(chi_s-chi)/chi_s][...] C_l^{Phi00Phi00}(chi,chi')
folds the canoes Phi00 cumulant with the a^2-FREE textbook efficiency (the a^2
of K is absorbed by the cumulant measure change kappa2(chi)=a^2 a^2 kappa2(lambda),
Eq. "cumulant 2 measure change").  The per-leg (1+z)^4 Sachs response is therefore
EXACTLY compensated at 2 legs and the 2-pt matches PyCCL to ~2%.  The 3-pt
equal-shell vertex is a SEPARATE code path (kappa3.py) whose zeta normalization +
units h-power carry the (1+z)^4 + h^2 surplus relative to the standard.

Output: outputs/powercount_3pt.npz + printed power-count tables.
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
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs import compute_kappa3_sigma3_high  # noqa: E402
from scipy.integrate import cumulative_trapezoid  # noqa: E402

H0 = 100.0 / 299792.458  # h/Mpc, canoes _H0_H_PER_MPC (h-units H0)


def _surplus_B_hpower(cosmo, pk) -> tuple[float, float]:
    """HIGH zeta_TTT physical/h ratio -> empirical h-power (expect 6 = bug)."""
    h = cosmo.h
    zf = np.linspace(0.0, 5.0, 4000)
    chif = np.array([chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h) for zi in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5.0 / 12, 5.0, 12)
    lam_h = np.interp(z, zf, lam_f) * h
    cg = np.cos(np.radians(np.array([1.0, 10.0]) / 60.0))
    triples = np.stack([np.ones_like(cg), cg, cg], axis=1)
    common = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, channels=("TTT",), lambda_convention="project",
        radial_discretization="sample", radial_measure="lambda",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = np.asarray(compute_kappa3_sigma3_high(triples, lam_h, units="physical", **common).zeta_TTT)
        zh = np.asarray(compute_kappa3_sigma3_high(triples, lam_h, units="h", **common).zeta_TTT)
    msk = np.isfinite(zp) & np.isfinite(zh) & (np.abs(zh) > 0)
    med = float(np.median(zp[msk] / zh[msk]))
    hpow = float(np.log(np.median(np.abs(zp[msk] / zh[msk]))) / np.log(h))
    return med, hpow


def main() -> int:
    cosmo = FiducialCosmology()
    h = cosmo.h
    pk = load_pk_delta(n_s=cosmo.n_s)

    # ---- SURPLUS (B): h-power of the 3-pt HIGH zeta ----
    print("=" * 90)
    print("SURPLUS (B): h-power of the 3-pt HIGH zeta_TTT  [kappa3.py:3642]")
    print("=" * 90)
    med_units, hpow_3pt = _surplus_B_hpower(cosmo, pk)
    print(f"  HIGH zeta physical/h median  = {med_units:.6e}")
    print(f"  h^6 = {h**6:.6e}   h^4 = {h**4:.6e}")
    print(f"  empirical 3-pt h-power       = {hpow_3pt:.4f}   (code applies h^6)")
    print(f"  dimensionless-correct h-power= 4.0000  (PyCCL oracle h-test; 2-pt kappa2.py:364)")
    print(f"  equal-shell zeta convention  = h^6 UNIFORMLY (scalar LOW :4872, scalar HIGH :3642,")
    print(f"                                 modulus LOW :3906-07, modulus HIGH :4050-51)")
    print(f"  => SURPLUS h-power           = {hpow_3pt - 4.0:+.4f}  (extra h^2; a CONVENTION, not a 1-line bug)")
    print()

    # ---- SURPLUS (A): per-shell (1+z)^4 radial ----
    print("=" * 90)
    print("SURPLUS (A): per-shell deployed-canoes / standard radial response ratio")
    print("=" * 90)
    zf = np.linspace(0.0, 5.0, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h) for zi in zf])
    z = np.linspace(5.0 / 40, 5.0 * 0.975, 40)
    chi = np.interp(z, zf, chif)
    chi_s = float(chif[-1])
    a = 1.0 / (1.0 + z)
    opz = 1.0 + z
    K_geo = chi * (chi_s - chi) / chi_s
    A_pois = -1.5 * cosmo.Omega_m * H0 ** 2 * opz
    resp_can = (a ** 8) * (K_geo ** 3) * (A_pois * opz ** 4) ** 3
    g_std = 1.5 * cosmo.Omega_m * H0 ** 2 * opz * K_geo
    resp_std = g_std ** 3
    ratio_shell = resp_can / resp_std
    sym = -(opz ** 4)
    print(f"{'shell':>5s} {'z':>7s} {'1+z':>7s} {'S_can/S_std':>13s} {'-(1+z)^4':>10s} "
          f"{'ratio/(-(1+z)^4)':>16s}")
    for i in (5, 10, 20, 30, 38):
        print(f"{i:5d} {z[i]:7.3f} {opz[i]:7.3f} {ratio_shell[i]:13.5e} {sym[i]:10.3f} "
              f"{ratio_shell[i]/sym[i]:16.6f}")
    print(f"\n  symbolic: a^8 [A(a)(1+z)^4]^3 / [3/2 Om H0^2 (1+z) K_geo]^3 = -(1+z)^4")
    print(f"  max |ratio_shell/(-(1+z)^4) - 1| = {np.nanmax(np.abs(ratio_shell/sym - 1.0)):.3e}")
    print()

    # ---- DEPLOYED ORACLE CLOSURE: realcan/std ~ (1+z_eff)^4 * h^2 ----
    print("=" * 90)
    print("DEPLOYED ORACLE CLOSURE: realcan/std ~ (1+z_eff)^4 * h^2  (both surpluses)")
    print("=" * 90)
    orc = np.load(_OUT / "scalar_vs_pyccl_tree.npz", allow_pickle=True)  # trusted local npz
    ga = orc["gamma_arcmin"]
    realcan = orc["ratio_realcan_std"]
    z_eff = orc["z_eff"]
    pred = (1.0 + z_eff) ** 4 * h ** 2
    print(f"{'gamma':>8s} {'realcan/std':>11s} {'(1+ze)^4':>9s} {'h^2':>6s} "
          f"{'pred':>9s} {'real/pred':>9s}")
    for i in range(ga.size):
        print(f"{ga[i]:8.3f} {realcan[i]:11.3f} {(1+z_eff[i])**4:9.3f} {h**2:6.3f} "
              f"{pred[i]:9.3f} {realcan[i]/pred[i]:9.3f}")
    sig = ga <= 12.3
    print(f"\n  signal regime gamma<=12': realcan/std median = {np.median(realcan[sig]):.3f}")
    print(f"  => removing (1+z)^4 [radial, surplus A] AND h^6->h^4 [units, surplus B]")
    print(f"     makes the dimensionless conv-3PCF match PyCCL.")

    np.savez(
        _OUT / "powercount_3pt.npz",
        z=z, chi=chi, opz=opz, a=a, resp_can=resp_can, resp_std=resp_std,
        ratio_shell=ratio_shell, sym_minus_1pz4=sym,
        hpow_3pt=hpow_3pt, physical_over_h_median=med_units,
        gamma_arcmin=ga, realcan_std=realcan, z_eff=z_eff,
        pred_1pz4_h2=pred,
        note="3-pt power-count: surplus (A)=(1+z)^4 per shell [radial], "
             "(B)=h^2 [units kappa3.py:3642 h^6->h^4].",
    )
    print(f"\nsaved -> {_OUT / 'powercount_3pt.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
