"""STAGE 0 (ours side): scalar convergence 3PCF from the driving-field 3-cumulant.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

GOAL
----
Build a GENUINE scalar convergence 3PCF on the COLLAPSED (squeezed, degenerate-
triangle) configuration ``(1, cos gamma, cos gamma)`` at a single source plane
z_s = 5, directly from the canoes equal-shell driving-field 3-cumulant
``zeta_TTT(gamma, lambda)`` (appendix `eq: appendix equal shell approx`, D7).

The fold is the 3-point analogue of the convergence 2PCF fold
(cosmology.tex `eq: xi kappa folded` + appendix `eq: appendix kappa3 cross shell`
+ `eq: appendix equal shell approx`):

    Z_kappa(gamma) = int_0^{lambda_s} dlambda  K(lambda, lambda_s)^3
                       (1+z)^-4  zeta_TTT(gamma, lambda)

with the single-source-plane convergence efficiency (cosmology.tex
`eq: K comoving kernel`, lines 1064-1068):

    K(lambda, lambda_s) = Dbar^2(lambda) int_lambda^{lambda_s} dlambda1 / Dbar^2(lambda1)
                        = a^2(chi) * chi * (chi_s - chi) / chi_s,   Dbar = a*chi.

THE (1+z)^-4 BORN PER-LEG REDUCTION (2026-06-09; fix/spin2-zetaD-2_2_-2-norm).
The single-plane K^3 fold against canoes' per-leg (1+z)^4 Sachs response
over-counts the textbook convergence by a per-shell (1+z)^(2n-2) = (1+z)^4 for
n=3 legs (derived from first principles in ``_born_leg_reduction`` and
anchored decisively against a CCL-native 2-pt control in
``b2_kfold_2pt_vs_ccl.py``: the 2-leg analogue over-counts by exactly (1+z)^2,
and dividing pointwise by (1+z)^2 restores Cl(hand)/Cl(CCL) = 0.99).  We
therefore multiply the fold integrand pointwise by (1+z(lambda))^-4 -- the Born
E=E0/a per-leg reduction that the validated cross-shell production
(corr_op O0 / equal_time_limber FK; O0/CCL = 0.855, no excess) already
performs.  canoes (zeta) and the kernel K are unchanged; this is an
analysis-side bookkeeping reduction only.

CONVENTIONS (asserted in code):
  * zeta is built fresh per the canonical recipe in
    callables/kappa3_vertex/equal_time_limber/build_equal_time_limber_table.py:
      LOW  : compute_kappa3_zeta_table(..., ell_max=60, Nmax=128, spt_kind="tree_phi",
             radial_measure="lambda", units="physical", lambda_convention="project",
             channels=("TTT",))
      HIGH : compute_kappa3_sigma3_high(..., ell_cut=60, ell_high_max=1000,
             n_ell=96, n_phi=64, radial_measure="lambda", units="physical",
             channels=("TTT",))
      combine: kappa3_combine_low_high(low, high).
  * radial_measure="lambda": zeta is a lambda-density 3-leg object.  Do NOT apply
    (dchi/dlambda)^3 (per-leg (1+z)^4 Sachs prefactor already bakes in the
    Jacobian -- radial-measure memory 2026-05-17).
  * units="physical": zeta carries h^4, chi/lambda carry h^-1; the kernel K is
    built in physical Mpc, so the fold is fully physical-unit consistent.
    (2026-06-09 dimensional correction: the equal-shell 3-cumulant is natively
    (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, so units=physical applies h^4, not h^6.)

We also run a CHI-route cross-check (integrate in chi instead of lambda; each leg
picks up a^2(chi), kernel uses chi(chi_s-chi)/chi_s without the a^2 prefactor) and
verify the two routes agree.

Output: outputs/stage0_ours.npz with the gamma grid, lambda grid, zeta_TTT(gamma,
lambda), K(lambda), and the folded Z_kappa(gamma) (lambda route + chi route).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.cosmo.background import chi as chi_of_z  # noqa: E402

Z_SOURCE = 5.0
_DEPLOYED_NPZ = _ETL / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"

# Build zeta on a DENSE affine grid that integrates well to lambda_s(z=5).  We
# reuse the deployed table's z-shell spacing but truncate at z_s=5 and add a
# few intermediate shells for a smooth quadrature.
_N_LAMBDA = 40  # affine-grid resolution (z in (0, 5])


def _squeezed_triples(gamma_rad: np.ndarray) -> np.ndarray:
    """(1, cos gamma, cos gamma): two coincident directions, third at gamma."""
    cg = np.cos(gamma_rad)
    return np.stack([np.ones_like(cg), cg, cg], axis=1)  # (n_gamma, 3)


def _build_lambda_grid(cosmo, n_lambda: int):
    """Affine grid lambda(z) for z in (0, z_s] with lambda = int_0^z a^2 dchi.

    Returns (z, chi_phys, lambda_phys, lambda_h) all length n_lambda, plus the
    source-plane (z_s, chi_s, lambda_s).
    """
    from scipy.integrate import cumulative_trapezoid

    Om, h = cosmo.Omega_m, cosmo.h
    # Fine z-grid 0..z_s for the affine integral.
    zf = np.linspace(0.0, Z_SOURCE, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])  # physical Mpc
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)  # physical Mpc
    # Target shell redshifts: linear in z over (0, z_s], avoid z=0 (chi=0).
    z = np.linspace(Z_SOURCE / n_lambda, Z_SOURCE, n_lambda)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * h
    chi_s = float(chif[-1])
    lam_s = float(lam_f[-1])
    return z, chi_phys, lam_phys, lam_h, Z_SOURCE, chi_s, lam_s


def _convergence_kernel(chi_phys: np.ndarray, z: np.ndarray, chi_s: float):
    """K(lambda, lambda_s) = a^2(chi) * chi * (chi_s - chi)/chi_s  (physical Mpc)."""
    a = 1.0 / (1.0 + z)
    return a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s


def _born_leg_reduction(z: np.ndarray, n_legs: int) -> np.ndarray:
    """Per-shell Born per-leg photon-energy reduction factor (1+z)^-(2 n_legs - 2).

    THE (1+z)^(2n-2) FIX (2026-06-09; branch fix/spin2-zetaD-2_2_-2-norm).
    ------------------------------------------------------------------------
    The hand-rolled single-plane fold ``int dlambda K(lambda)^n zeta(lambda)``
    folds the *collapsed* single-source-plane convergence efficiency
    ``K = a^2 chi (chi_s-chi)/chi_s`` against the per-leg Sachs response that
    canoes bakes into zeta (the (1+z)^4 photon-energy / Poisson factor per leg).
    A first-principles match of this fold's per-shell implied convergence Cl
    integrand to the textbook CCL Limber integrand
    ``C_l^kk = int dchi W(chi)^2/chi^2 P_delta``, with the textbook lensing
    efficiency ``W(chi) = (3/2) Om H0^2 (1+z) g(chi)`` (ONE (1+z) per leg, the
    Born E=E0/a reduction), gives the EXACT pointwise over-count

        integrand_hand / integrand_textbook = a^(2n+2) (1+z)^(4n)
                                             = (1+z)^(2n-2)   per shell,

    independent of ell.  (Derivation: per leg the hand fold carries (1+z)^4 from
    canoes' Sachs prefactor and (1+z) from A(a)=-3/2 Om H0^2 (1+z); K=a^2 g and
    the affine measure dlambda=a^2 dchi together refund a^(2n+2); the textbook
    carries only (1+z)^n for n legs.)

    DECISIVE CCL-NATIVE ANCHOR (b2_kfold_2pt_vs_ccl.py, PyCCL env, ran live):
    the 2-leg analogue of this very fold gives Z2_hand/CCL = 6.0 and a grid-
    independent implied-Cl(hand)/Cl(CCL) = 3.6..7.1 across ell=30..3000, i.e.
    eff (1+z) = sqrt = 1.9..2.66 = per-shell (1+z) at the kernel-relevant z
    (z~0.9..1.7, NOT z_s=5).  Dividing the 2-pt integrand pointwise by (1+z)^2
    drives the implied-Cl ratio to 0.99 (flat across all ell).  The validated
    CROSS-SHELL corr_op O0/CCL(1') = 0.855 (~1, NO excess) is the control: the
    production sft-wick cross-shell R/C radial fold distributes (1+z) correctly
    across independent source positions and needs NO such factor; only the
    hand-rolled single-plane K^n collapse over-counts.  The same-family
    stage0_spt_reference (which shares this K^n fold and the same per-leg
    (1+z)^4) carries the SAME excess, so the old ours/ref cross-check could not
    see it -- only the CCL-native anchor + corr_op control expose it.

    Therefore the fold integrand is multiplied pointwise (per shell) by
    (1+z(lambda))^-(2 n_legs - 2): the Born per-leg reduction that the validated
    cross-shell production already performs.  This is the ONLY change; the
    kernel K, the affine measure dlambda=a^2 dchi, and canoes' zeta are all
    left untouched (canoes is read-only and correct).
    """
    return (1.0 + z) ** (-(2 * n_legs - 2))


def main() -> int:
    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)

    # --- assert the deployed table is units==physical, measure==lambda ----------
    import json
    dep = np.load(_DEPLOYED_NPZ, allow_pickle=True)
    meta = json.loads(dep["cosmo_meta"].tobytes().decode("utf-8", errors="replace"))
    assert meta.get("units") == "physical", f"deployed units != physical: {meta.get('units')}"
    assert meta.get("radial_measure") == "lambda", f"deployed measure != lambda"
    print(f"[ours] deployed-table assert OK: units={meta['units']}, "
          f"measure={meta['radial_measure']}, has_mod={meta.get('has_modulus_channels')}")

    # --- affine grid + source plane --------------------------------------------
    z, chi_phys, lam_phys, lam_h, z_s, chi_s, lam_s = _build_lambda_grid(cosmo, _N_LAMBDA)
    print(f"[ours] z_s={z_s}, chi_s={chi_s:.2f} Mpc, lambda_s={lam_s:.2f} Mpc")
    print(f"[ours] affine grid: n={_N_LAMBDA}, z in [{z[0]:.3f}, {z[-1]:.3f}], "
          f"chi in [{chi_phys[0]:.1f}, {chi_phys[-1]:.1f}] Mpc")

    # --- gamma grid: single 1' first, then a log-spaced sweep ------------------
    gamma_single = np.array([1.0])  # arcmin
    gamma_sweep = np.geomspace(0.5, 300.0, 9)  # arcmin (0.5' .. few hundred ')
    gamma_arcmin = np.unique(np.concatenate([gamma_single, gamma_sweep]))
    gamma_rad = np.radians(gamma_arcmin / 60.0)
    triples = _squeezed_triples(gamma_rad)
    print(f"[ours] gamma grid (arcmin): {np.array2string(gamma_arcmin, precision=3)}")

    common_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda", channels=("TTT",), parallelism="sequential",
        shell_block_size=4, triple_chunk_size=4096, n_workers=1,
        fuse_tree_phi_channels=True,
    )
    common_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, channels=("TTT",), units="physical",
        lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda",
    )

    print("[ours] computing LOW (exact Wigner-3j, ell<=60) ...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low = compute_kappa3_zeta_table(triples, lam_h, **common_low)
        print("[ours] computing HIGH (flat-sky Limber, (60,1000]) ...")
        high = compute_kappa3_sigma3_high(triples, lam_h, **common_high)
    combined = kappa3_combine_low_high(low, high)
    zeta_TTT = np.asarray(combined.zeta_TTT, dtype=np.float64)  # (n_gamma, n_lambda)
    print(f"[ours] zeta_TTT shape={zeta_TTT.shape}, "
          f"finite={int(np.sum(np.isfinite(zeta_TTT)))}/{zeta_TTT.size}")

    # --- fold in lambda: Z = int dlambda K^3 (1+z)^-4 zeta_TTT -------------------
    # The (1+z)^-(2*3-2) = (1+z)^-4 factor is the Born per-leg photon-energy
    # reduction (see _born_leg_reduction); it removes the per-shell (1+z)^4
    # over-count of the single-plane K^3 fold so this scalar 3PCF reduces to the
    # textbook convergence response that the validated cross-shell production
    # (corr_op O0 / equal_time_limber FK) carries.
    K = _convergence_kernel(chi_phys, z, chi_s)  # physical Mpc, (n_lambda,)
    born3 = _born_leg_reduction(z, n_legs=3)  # (1+z)^-4, (n_lambda,)
    integrand_lambda = (K ** 3 * born3)[None, :] * zeta_TTT  # (n_gamma, n_lambda)
    Z_lambda = np.trapezoid(integrand_lambda, lam_phys, axis=1)  # physical-unit 3PCF

    # --- fold in chi (cross-check): pure change of integration variable -------
    #   Z = int dlambda K_lambda^3 zeta.  Substitute dlambda = a^2(chi) dchi and
    #   K_lambda = a^2(chi) K_chi with K_chi = chi(chi_s-chi)/chi_s (textbook
    #   efficiency, no a^2 prefactor).  Then
    #       Z = int dchi a^2 (a^2 K_chi)^3 zeta = int dchi a^8 K_chi^3 zeta,
    #   where zeta is the SAME lambda-density value (it is NOT re-scaled; the
    #   measure change supplies the one extra a^2 from dlambda=a^2 dchi).  This
    #   is the SAME integral as the lambda route, evaluated with an independent
    #   quadrature -> a quadrature-consistency check, not a new physics input.
    a = 1.0 / (1.0 + z)
    K_chi = chi_phys * (chi_s - chi_phys) / chi_s  # textbook efficiency (no a^2)
    integrand_chi = (K_chi ** 3 * a ** 8 * born3)[None, :] * zeta_TTT
    Z_chi = np.trapezoid(integrand_chi, chi_phys, axis=1)

    rel_route = np.abs(Z_lambda - Z_chi) / np.maximum(np.abs(Z_lambda), 1e-300)
    print(f"[ours] lambda-route vs chi-route max rel diff: {np.max(rel_route):.2e}")

    out = _HERE / "outputs" / "stage0_ours.npz"
    np.savez(
        out,
        gamma_arcmin=gamma_arcmin, gamma_rad=gamma_rad,
        z=z, chi_phys=chi_phys, lam_phys=lam_phys, K=K,
        zeta_TTT=zeta_TTT, Z_lambda=Z_lambda, Z_chi=Z_chi,
        z_s=z_s, chi_s=chi_s, lam_s=lam_s,
        Omega_m=cosmo.Omega_m, h=cosmo.h, n_s=cosmo.n_s,
    )
    print(f"[ours] saved -> {out}")

    # --- report the single 1' value FIRST --------------------------------------
    i1 = int(np.argmin(np.abs(gamma_arcmin - 1.0)))
    print()
    print(f"[ours] === SINGLE-GAMMA (gamma={gamma_arcmin[i1]:.3f}') ===")
    print(f"[ours]   Z_kappa (lambda route) = {Z_lambda[i1]:+.6e}")
    print(f"[ours]   Z_kappa (chi route)    = {Z_chi[i1]:+.6e}")
    print(f"[ours]   sign = {'+' if Z_lambda[i1] > 0 else '-'}")
    print()
    print("[ours] === SWEEP ===")
    for ig in range(gamma_arcmin.size):
        print(f"[ours]   gamma={gamma_arcmin[ig]:8.3f}'  Z_lambda={Z_lambda[ig]:+.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
