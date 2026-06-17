"""STAGE 1 v2 (ours side) -- JACFIX: (1+z)^-4 comes from the canoes vertex, NOT a manual patch.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

This is a byte-for-byte copy of stage1_ours_v2_FIXED.py EXCEPT for the radial
fold weight: the manual per-shell Born per-leg reduction (1+z)^-4 (the `born3`
factor) is REMOVED (set to ones).

WHY (2026-06-09 radial-Jacobian fix landed in canoes):
  The canoes equal-shell kappa3 vertex was corrected so that
  `radial_measure="lambda"` now returns the (1+z)^-4 = (d lambda/d chi)^2
  Jacobian (the two collapsed comoving matter deltas) instead of the previous
  `ones` -- see canoes/src/canoes/sachs/kappa3.py::_kappa3_radial_density_conversion.
  Every compute_kappa3_* call in build_zeta_channels() below uses
  radial_measure="lambda", so the FRESH canoes queries now ALREADY carry
  (1+z)^-4 internally.  Applying the manual `born3 = (1+z)^-4` in the fold on top
  of that would DOUBLE-COUNT the factor.  We therefore set born3 = ones here: the
  (1+z)^-4 enters ONCE, from the vertex.

The deployed equal_time_limber kappa3 table (the FIXED one, Jun 9) likewise
carries this (1+z)^-4; a consistency probe established the deployed callable ==
a fresh targeted canoes query bit-for-bit, so these fresh queries reflect the
deployed object.

The net result should be ~unchanged from the prior validated stage1_ours_v2_FIXED
run (the factor just moved from the manual patch into the vertex): the
frame-invariant canoes/fastnc (and canoes/ref) ratio should still land near
~0.70x (median), NOT 7.77x.

Output: outputs/stage1_ours_v2_JACFIX.npz (same schema as stage1_ours_v2_FIXED.npz).
"""
from __future__ import annotations

import json
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

_GRID_JSON = _HERE / "stage1_grid.json"
_N_LAMBDA = 40


# ---------------------------------------------------------------------------
# SLOT-LABEL + CENTROID-PHASE FIX (2026-06-10, J_2 angular-residual closeout)
#
# Probe C (scripts/derive_natural_component_leg_map.wl) established, FIRST-
# PRINCIPLES, the Schneider-Lombardi / fastnc convention for the SAS isoceles
# family (t1=t2=gamma, t3=2 gamma sin(phi/2)):
#   * the APEX vertex is self-mirror, so its conjugated natural component must
#     be REAL (fastnc: Gamma^1 real to 1e-16), and
#   * the two BASE vertices are a mirror pair, so their conjugated components
#     are complex conjugates (fastnc: Gamma^2 = conj(Gamma^3)).
# The xAct derivation gives the per-vertex centroid reference angle EXACTLY:
#     alpha_apex   = 0                              (=> apex phase exp(-2i*0)=1, REAL)
#     alpha_baseB  = -pi + arctan(3 tan(phi/2))
#     alpha_baseC  = +pi - arctan(3 tan(phi/2))     (= -alpha_baseB mod 2pi).
#
# In the SAS cosine triple (cos t3, cos gamma, cos gamma), legs 1,2 are the BASE
# endpoints (separation t3) and leg 3 is the APEX.  The previous `_amb` great-
# circle->centroid angles were NOT these centroid angles, so even the correctly-
# bound apex channel was not real.  We therefore use the centroid angles directly.
#
# Per-vertex centroid angle, indexed by physical leg (1=baseB, 2=baseC, 3=apex):
def _centroid_alpha(leg: int, phi: np.ndarray) -> np.ndarray:
    """Centroid projection reference angle alpha_leg(phi) [rad], from the xAct
    derivation derive_natural_component_leg_map.wl section 0/5."""
    arc = np.arctan2(3.0 * np.tan(np.asarray(phi, dtype=np.float64) / 2.0), 1.0)
    if leg == 3:      # APEX (self-mirror, on the symmetry axis)
        return np.zeros_like(arc)
    if leg == 1:      # BASE vertex B
        return -np.pi + arc
    if leg == 2:      # BASE vertex C (= -alpha_baseB)
        return np.pi - arc
    raise ValueError(leg)


def centroid_phase(conj_leg: int, phi: np.ndarray) -> np.ndarray:
    """Spin-2 centroid projection phase for the natural component that conjugates
    physical leg `conj_leg` (s_j = -1 on the conjugated leg, +1 elsewhere):
        prod_j exp(-2 i s_j alpha_j).
    conj_leg=0 means the all-unconjugated PPP component (s_j = +1 for all)."""
    a1 = _centroid_alpha(1, phi)
    a2 = _centroid_alpha(2, phi)
    a3 = _centroid_alpha(3, phi)
    s = {0: (1, 1, 1), 1: (-1, 1, 1), 2: (1, -1, 1), 3: (1, 1, -1)}[conj_leg]
    return np.exp(-2j * (s[0] * a1 + s[1] * a2 + s[2] * a3))


def _build_lambda_grid(cosmo, z_source: float, n_lambda: int):
    from scipy.integrate import cumulative_trapezoid

    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, z_source, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_source / n_lambda, z_source, n_lambda)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * h
    chi_s = float(chif[-1])
    lam_s = float(lam_f[-1])
    return z, chi_phys, lam_phys, lam_h, chi_s, lam_s


def _convergence_kernel(chi_phys, z, chi_s):
    a = 1.0 / (1.0 + z)
    return a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s


def _sas_cosine_triples(gamma_rad: float, phi_rad: np.ndarray) -> np.ndarray:
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad)
    ct3 = np.cos(t3)
    return np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)


# Spin-weight triples keyed by the PHYSICAL leg that is conjugated (the -2 slot):
#   leg 1 = BASE vertex B, leg 2 = BASE vertex C, leg 3 = APEX.
# This key is the SAME index passed to centroid_phase(), so the spin re-binding
# and the centroid-projection phase stay consistent under the slot relabel below.
_DMOD_SPINS = {
    1: (-2, 2, 2),   # conjugate BASE leg 1
    2: (2, -2, 2),   # conjugate BASE leg 2
    3: (2, 2, -2),   # conjugate APEX leg 3
}


def main() -> int:
    grid = json.loads(_GRID_JSON.read_text())
    z_source = float(grid["z_s"])
    gamma_arcmin = np.array(grid["gamma_arcmin"], dtype=np.float64)
    phi_deg = np.array(grid["phi_deg"], dtype=np.float64)
    phi_rad = np.radians(phi_deg)

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)

    z, chi_phys, lam_phys, lam_h, chi_s, lam_s = _build_lambda_grid(
        cosmo, z_source, _N_LAMBDA)
    K = _convergence_kernel(chi_phys, z, chi_s)
    print(f"[JACFIX] z_s={z_source}, chi_s={chi_s:.2f} Mpc, lam_s={lam_s:.2f} Mpc")

    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    from canoes.sachs.kappa3 import (  # noqa: E402
        compute_kappa3_mod_zeta_equal_time,
        compute_kappa3_mod_sigma3_high,
    )

    common_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda", channels=("PPP",), parallelism="sequential",
        shell_block_size=4, triple_chunk_size=4096, n_workers=1,
        fuse_tree_phi_channels=True,
    )
    common_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, channels=("PPP",), units="physical",
        lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda",
    )
    mod_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_measure="lambda",
    )
    mod_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, units="physical", lambda_convention="project",
        radial_measure="lambda",
    )

    def build_zeta_channels(triples):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            low = compute_kappa3_zeta_table(triples, lam_h, **common_low)
            high = compute_kappa3_sigma3_high(triples, lam_h, **common_high)
            comb = kappa3_combine_low_high(low, high)
            zPPP = np.asarray(comb.zeta_PPP, dtype=np.float64)
            zD = {}
            for mu, sw in _DMOD_SPINS.items():
                dlo = compute_kappa3_mod_zeta_equal_time(
                    triples, lam_h, _spin_weights_dmod=sw, **mod_low)[1]
                dhi = compute_kappa3_mod_sigma3_high(
                    triples, lam_h, _spin_weights_dmod=sw, **mod_high)[1]
                zD[mu] = np.asarray(dlo + dhi, dtype=np.float64)
        return zPPP, zD[1], zD[2], zD[3]

    # JACFIX: NO manual (1+z)^-4 here.  The canoes vertex (radial_measure="lambda")
    # now carries the (1+z)^-4 = (d lambda/d chi)^2 Jacobian internally, so applying
    # it again in the fold would double-count.  born3 = ones => factor enters ONCE,
    # from the vertex.
    born3 = np.ones_like(z)

    def fold(zeta):
        return np.trapezoid((K ** 3 * born3)[None, :] * zeta, lam_phys, axis=1)

    def natural_components(gamma_arcmin_val, phi_rad_arr):
        gamma_rad = np.radians(gamma_arcmin_val / 60.0)
        triples = _sas_cosine_triples(gamma_rad, phi_rad_arr)
        zPPP, zD1, zD2, zD3 = build_zeta_channels(triples)
        # fD1,fD2 = BASE-conjugated real folds; fD3 = APEX-conjugated real fold.
        fPPP, fD1, fD2, fD3 = fold(zPPP), fold(zD1), fold(zD2), fold(zD3)
        sign = -1.0
        # SLOT RELABEL (Schneider-Lombardi / fastnc SAS convention):
        #   Gamma^0 <- all-unconjugated PPP (real on isoceles).
        #   Gamma^1 <- APEX-conjugated channel (leg 3): alpha_apex=0 => REAL.
        #   Gamma^2 <- BASE leg-1 conjugated ; Gamma^3 <- BASE leg-2 conjugated,
        #   so Gamma^2 = conj(Gamma^3) (the two base vertices are a mirror pair).
        # Both the spin re-binding (_DMOD_SPINS key) and the centroid phase
        # (centroid_phase conj_leg) use the SAME physical-leg index, so this is a
        # consistent simultaneous re-binding of spins AND phase, not a rename.
        G0 = sign * fPPP * centroid_phase(0, phi_rad_arr)
        G1 = sign * fD3 * centroid_phase(3, phi_rad_arr)   # APEX -> slot 1 (REAL)
        G2 = sign * fD1 * centroid_phase(1, phi_rad_arr)   # BASE leg 1 -> slot 2
        G3 = sign * fD2 * centroid_phase(2, phi_rad_arr)   # BASE leg 2 -> slot 3
        diag = dict(fPPP=fPPP, fD1=fD1, fD2=fD2, fD3=fD3)
        return (np.asarray(G0, np.complex128), np.asarray(G1, np.complex128),
                np.asarray(G2, np.complex128), np.asarray(G3, np.complex128), diag)

    n_g, n_p = gamma_arcmin.size, phi_deg.size
    Gamma0 = np.zeros((n_g, n_p), np.complex128)
    Gamma1 = np.zeros((n_g, n_p), np.complex128)
    Gamma2 = np.zeros((n_g, n_p), np.complex128)
    Gamma3 = np.zeros((n_g, n_p), np.complex128)
    diag_fPPP = np.zeros((n_g, n_p))
    diag_fD = np.zeros((n_g, n_p))
    for ig, gam in enumerate(gamma_arcmin):
        G0, G1, G2, G3, diag = natural_components(gam, phi_rad)
        Gamma0[ig], Gamma1[ig], Gamma2[ig], Gamma3[ig] = G0, G1, G2, G3
        diag_fPPP[ig] = diag["fPPP"]
        diag_fD[ig] = diag["fD1"]
        print(f"[JACFIX]  gamma={gam:7.2f}' done  |G1_apex(phi120)|={abs(G1[-1]):.4e}"
              f"  ImG1/|G1|(phi60)={abs(G1[4].imag)/max(abs(G1[4]),1e-300):.1e}")

    out = _HERE / "outputs" / "stage1_ours_v2_JACFIX.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    convention_note = (
        "JACFIX rebuild of stage1_ours_v2_FIXED with the manual (1+z)^-4 Born "
        "per-leg reduction REMOVED (born3=ones). The (1+z)^-4 = (d lambda/d chi)^2 "
        "Jacobian now comes ONLY from the canoes equal-shell kappa3 vertex "
        "(radial_measure='lambda' returns (1+z)^-4 after the 2026-06-09 "
        "_kappa3_radial_density_conversion fix), avoiding the double-count that "
        "applying born3 on top would cause. SLOT-LABEL + CENTROID-PHASE FIX "
        "(2026-06-10): the natural-component slots follow the Schneider-Lombardi/"
        "fastnc SAS convention -- the APEX-conjugated channel (leg 3, spin "
        "(2,2,-2)) feeds Gamma^1 (real on isoceles since alpha_apex=0), and the "
        "two BASE channels feed {Gamma^2,Gamma^3} so Gamma^2=conj(Gamma^3). The "
        "spin-2 projection phase uses the xAct-derived per-vertex CENTROID angles "
        "(derive_natural_component_leg_map.wl) directly, replacing the previous "
        "_amb great-circle->centroid angles. "
        "zeta_D via compute_kappa3_mod_zeta_equal_time/_sigma3_high; PPP via "
        "compute_kappa3_zeta_table/_sigma3_high. units=physical, "
        "radial_measure=lambda, fold int dlambda K^3 zeta (no manual born3), "
        "K=a^2 chi(chi_s-chi)/chi_s, z_s=5. Net result ~unchanged vs the prior "
        "stage1_ours_v2_FIXED run (factor moved from manual patch into vertex): "
        "frame-invariant canoes/fastnc ~0.70x median expected, NOT 7.77x."
    )
    np.savez(
        out,
        gamma_arcmin=gamma_arcmin, phi_deg=phi_deg,
        Gamma0=Gamma0, Gamma1=Gamma1, Gamma2=Gamma2, Gamma3=Gamma3,
        projection="cent", convention_note=convention_note,
        z_s=z_source, chi_s=chi_s, lam_s=lam_s,
        Omega_m=cosmo.Omega_m, h=cosmo.h, n_s=cosmo.n_s,
        diag_fPPP=diag_fPPP, diag_fD=diag_fD,
    )
    print(f"\n[JACFIX] saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
