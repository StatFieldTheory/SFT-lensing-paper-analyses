"""STAGE 1 v2 (ours side): cosmic-SHEAR 3PCF natural components Gamma^0..3
from the driving-field spin-2 cumulants, using the DERIVED great-circle ->
centroid projection map (scripts/derive_shear3pcf_natural_components.wl).

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

WHAT CHANGED vs stage1_ours.py (v1)
----------------------------------
v1 multiplied the canonical-frame cumulants by fastnc's OWN x2cent(mu,t1,t2,phi)
phase.  That is WRONG: x2cent converts fastnc's FFT-native x-projection to the
centroid projection; it is NOT the bridge from canoes' canonical GREAT-CIRCLE
frame to the centroid projection.  v2 replaces x2cent by the per-leg spin-2
rotation DERIVED from first principles:

    Gamma^mu = Z_channel(mu) * exp(-2 i Sum_j s_j^mu (a_j - b_j))

  - a_j = centroid projection reference angle at vertex j (vertex->centroid),
  - b_j = canonical great-circle reference angle (geodesic to leg 1),
  - s_j^mu = +1 for an unconjugated leg, -1 for the conjugated leg of Gamma^mu,
  - Z_channel: Gamma^0 <- zeta_PPP (all +2), Gamma^{1,2,3} <- zeta_D (one leg
    conjugated; the conjugated leg is leg mu).

The channel assignment (which cumulant feeds which Gamma^mu) is fixed by the
helicity bookkeeping confirmed against fastnc's own Bessel-order map
(fastnc.py:365): Gamma^0 = all-same-helicity (total spin 6), Gamma^{1,2,3} =
one helicity flipped (total spin 2).

The DERIVED closed-form phases (from the .wl, pure functions of phi):
    a1 - b1 = -arccot(3 tan(phi/2))
    a2 - b2 = -2pi + arccot(3 tan(phi/2))
    a3 - b3 = -pi/2 + arctan(cot(phi/2))
with leg-1,2 = base vertices, leg-3 = apex (SAS triple (cos t3, cos g, cos g)).

NOTE ON THE SHIPPED NPZ. `outputs/stage1_ours_v2.npz` was produced by reusing
v1's deployed-faithful cumulant folds (recovered as Gamma_v1 / x2cent, which is
exactly the real radial fold since |x2cent|=1) and re-applying the DERIVED phase
above -- the prompt explicitly permits reusing the v1 values ("they are
deployed-faithful per the consistency check").  This is bit-identical in
magnitude to a fresh recompute (the cumulants are unchanged; only the phase
differs) and avoids the ~hours-long 6-channel spin-aware re-build.  Running this
script top-to-bottom recomputes the channels from canoes for an independent
check (slow).

WARNING (2026-06-09; fix/spin2-zetaD-2_2_-2-norm). The `fold()` here NOW applies
the (1+z)^-4 Born per-leg reduction (see _born_leg_reduction).  The previously
shipped `stage1_ours_v2.npz` was built BEFORE this fix and therefore carries the
per-shell (1+z)^4 over-count; it is STALE.  Re-run this script (or, preferably,
stage1_ours_v2_FIXED.py which does a genuine fresh recompute) to regenerate the
corrected npz.  Reusing v1 folds is NO LONGER bit-identical because v1's folds
also predate the reduction.

HONEST CAVEAT (do not paper over): the per-component MAGNITUDES |Gamma^mu| are
set entirely by the channel cumulants and are phase-invariant; the projection
phase derived here cannot change them.  Whether the comparison CLOSES is decided
by stage1_compare_v2.py.  If the magnitudes do not close, the residual is a
genuine projection/normalization difference between canoes' great-circle
driving-field cumulants and fastnc's fixed-frame spin-2 projection of the scalar
convergence bispectrum, reported in stage1_compare_v2.py + the docs.
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
# DERIVED great-circle -> centroid spin-2 projection phase
# (scripts/derive_shear3pcf_natural_components.wl, section 8).  a_j - b_j are
# pure functions of phi; the per-leg helicity sign s_j^mu selects the channel.
# ---------------------------------------------------------------------------
def _amb(phi: np.ndarray):
    """Return (a1-b1, a2-b2, a3-b3) in radians, pure functions of phi.

    From the .wl:  arccot(x) = arctan(1/x);  cot(phi/2) = 1/tan(phi/2).
    """
    half = phi / 2.0
    d1 = -np.arctan2(1.0, 3.0 * np.tan(half))                 # -arccot(3 tan)
    d2 = -2.0 * np.pi + np.arctan2(1.0, 3.0 * np.tan(half))   # -2pi + arccot(3 tan)
    d3 = -0.5 * np.pi + np.arctan2(1.0, np.tan(half))         # -pi/2 + arctan(cot)
    return d1, d2, d3


def derived_phase(mu: int, phi: np.ndarray) -> np.ndarray:
    """Per-component great-circle->centroid spin-2 phase exp(-2i sum s_j (a_j-b_j))."""
    d1, d2, d3 = _amb(np.asarray(phi, dtype=np.float64))
    s = {0: (1, 1, 1), 1: (-1, 1, 1), 2: (1, -1, 1), 3: (1, 1, -1)}[mu]
    return np.exp(-2j * (s[0] * d1 + s[1] * d2 + s[2] * d3))


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


def _born_leg_reduction(z, n_legs: int = 3):
    """Per-shell Born per-leg reduction (1+z)^-(2 n_legs - 2); see stage0_ours.py.

    THE (1+z)^(2n-2) FIX (2026-06-09; fix/spin2-zetaD-2_2_-2-norm).  The single-
    plane K^n fold against canoes' per-leg (1+z)^4 Sachs response over-counts the
    textbook lensing response by a per-shell (1+z)^(2 n_legs - 2) (anchored
    against a CCL-native 2-pt control in b2_kfold_2pt_vs_ccl.py: the 2-leg fold
    over-counts by exactly (1+z)^2; dividing pointwise by (1+z)^2 restores
    Cl(hand)/Cl(CCL)=0.99).  Shear 3PCF = 3-leg fold -> (1+z)^-4.  Born E=E0/a
    per-leg reduction; canoes' zeta and kernel K unchanged.  Phase-independent:
    scales every Gamma^mu identically (overall normalization only).
    """
    return (1.0 + z) ** (-(2 * n_legs - 2))


def _sas_cosine_triples(gamma_rad: float, phi_rad: np.ndarray) -> np.ndarray:
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad)
    ct3 = np.cos(t3)
    return np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)


# leg conjugated by each modulus component: Gamma^1 conjugates leg... we keep the
# SAME helicity pattern as the phase sign s_j^mu: Gamma^mu conjugates leg mu.
_DMOD_SPINS = {
    1: (-2, 2, 2),   # Gamma^1: conjugate leg 1
    2: (2, -2, 2),   # Gamma^2: conjugate leg 2
    3: (2, 2, -2),   # Gamma^3: conjugate leg 3
}


def main() -> int:
    grid = json.loads(_GRID_JSON.read_text())
    z_source = float(grid["z_s"])
    gamma_arcmin = np.array(grid["gamma_arcmin"], dtype=np.float64)
    phi_deg = np.array(grid["phi_deg"], dtype=np.float64)
    phi_rad = np.radians(phi_deg)
    prim = grid["primary_point"]
    g_prim = float(prim["gamma_arcmin"])
    phi_prim = float(prim["phi_deg"])

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    assert abs(cosmo.Omega_m - grid["cosmology"]["Omega_m"]) < 1e-9
    assert abs(cosmo.h - grid["cosmology"]["h"]) < 1e-9

    z, chi_phys, lam_phys, lam_h, chi_s, lam_s = _build_lambda_grid(
        cosmo, z_source, _N_LAMBDA)
    K = _convergence_kernel(chi_phys, z, chi_s)
    print(f"[s1-ours-v2] z_s={z_source}, chi_s={chi_s:.2f} Mpc, lam_s={lam_s:.2f} Mpc")

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

    born3 = _born_leg_reduction(z, n_legs=3)  # (1+z)^-4 Born per-leg reduction

    def fold(zeta):
        return np.trapezoid((K ** 3 * born3)[None, :] * zeta, lam_phys, axis=1)

    def natural_components(gamma_arcmin_val, phi_rad_arr):
        gamma_rad = np.radians(gamma_arcmin_val / 60.0)
        triples = _sas_cosine_triples(gamma_rad, phi_rad_arr)
        zPPP, zD1, zD2, zD3 = build_zeta_channels(triples)
        fPPP, fD1, fD2, fD3 = fold(zPPP), fold(zD1), fold(zD2), fold(zD3)
        # overall sign (-1)^3 from gamma = -int Dsachs (matches STAGE-0 kappa sign)
        sign = -1.0
        G0 = sign * fPPP * derived_phase(0, phi_rad_arr)
        G1 = sign * fD1 * derived_phase(1, phi_rad_arr)
        G2 = sign * fD2 * derived_phase(2, phi_rad_arr)
        G3 = sign * fD3 * derived_phase(3, phi_rad_arr)
        diag = dict(fPPP=fPPP, fD1=fD1, fD2=fD2, fD3=fD3)
        return (np.asarray(G0, np.complex128), np.asarray(G1, np.complex128),
                np.asarray(G2, np.complex128), np.asarray(G3, np.complex128), diag)

    # ---- PRIMARY point first --------------------------------------------
    print(f"\n[s1-ours-v2] === PRIMARY (gamma={g_prim}', phi={phi_prim}deg) ===")
    G0p, G1p, G2p, G3p, dp = natural_components(g_prim, np.array([np.radians(phi_prim)]))
    for mu, G in enumerate([G0p, G1p, G2p, G3p]):
        print(f"[s1-ours-v2]   Gamma^{mu} = {G[0]:+.6e}  |.|={abs(G[0]):.4e}")
    inv = np.sqrt(sum(abs(G[0]) ** 2 for G in [G0p, G1p, G2p, G3p]))
    print(f"[s1-ours-v2]   frame-invariant = {inv:.4e}")

    # ---- full grid ------------------------------------------------------
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
        print(f"[s1-ours-v2]  gamma={gam:7.2f}' done")

    out = _HERE / "outputs" / "stage1_ours_v2.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    convention_note = (
        "OUR-SIDE v2 cosmic-shear 3PCF natural components Gamma^0..3, "
        "Schneider-Lombardi, CENTROID projection. SAME driving-field channel "
        "assignment as v1 (Gamma^0 <- zeta_PPP all-+2; Gamma^{1,2,3} <- zeta_D "
        "with leg mu conjugated), BUT the projection phase is the DERIVED "
        "great-circle->centroid spin-2 rotation exp(-2i sum_j s_j (a_j-b_j)) "
        "(scripts/derive_shear3pcf_natural_components.wl), REPLACING v1's "
        "incorrect use of fastnc's x2cent (which bridges fastnc's x-projection, "
        "not canoes' canonical great-circle frame). a_j=centroid ref angle, "
        "b_j=geodesic-to-leg-1 ref angle; a_j-b_j pure functions of phi: "
        "a1-b1=-arccot(3 tan(phi/2)), a2-b2=-2pi+arccot(3 tan(phi/2)), "
        "a3-b3=-pi/2+arctan(cot(phi/2)). Channels from canoes "
        "compute_kappa3_mod_zeta_equal_time/_sigma3_high (zeta_D) and "
        "compute_kappa3_zeta_table/_sigma3_high (zeta_PPP), units=physical, "
        "radial_measure=lambda, fold int dlambda K^3 (1+z)^-4 zeta, "
        "K=a^2 chi(chi_s-chi)/chi_s, z_s=5. The (1+z)^-4 is the Born per-leg "
        "reduction removing the single-plane K^3 fold's per-shell (1+z)^4 "
        "over-count (fix/spin2-zetaD-2_2_-2-norm; CCL-native-anchored). "
        "Overall sign (-1)^3 from gamma=-int Dsachs. "
        "MAGNITUDE CAVEAT: |Gamma^mu| are phase-invariant (set by the channel "
        "cumulants); see stage1_compare_v2.npz for the honest closure assessment."
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
    print(f"\n[s1-ours-v2] saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
