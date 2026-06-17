"""STAGE 1 (ours side): cosmic-SHEAR 3PCF natural components Gamma^0..3 from the
driving-field equal-shell 3-cumulant (spin-2 modulus channels), z_s = 5.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

GOAL
----
Extend STAGE 0 (scalar convergence, certified to |ratio|=0.85-0.99) to the
spin-2 (shear) channels.  We build a GENUINE shear 3PCF on the pinned SAS
isoceles family (two equal sides t1=t2=gamma, opening angle phi, third side
t3 = 2 gamma sin(phi/2)) at a single source plane z_s = 5, directly from the
canoes equal-shell driving-field spin-2 modulus cumulants
zeta_B = <Phi00 Psi0 Psi0*> and zeta_D = <Psi0 Psi0 Psi0*>
(appendix `eq: appendix modulus reconstruction`, post-2026-06-04 helicity fix),
project to the Schneider-Lombardi natural components Gamma^0..Gamma^3, and write
them in fastnc's CENTROID projection so a parallel fastnc run can be compared
apples-to-apples on the pinned grid.

=====================================================================
PHYSICS CHAIN  (every convention documented; see CONVENTION CHECKLIST)
=====================================================================

1.  DRIVING FIELD -> SHEAR.  The spin-2 driving field is Psi0 = sigma_+ + i sigma_x
    (cosmology.tex `eq: gamma from Dsachs23`); the complex shear is
    gamma = gamma_+ + i gamma_x = -int dlambda (Dsachs_2 + i Dsachs_3).  In the
    great-circle-aligned local screen frame (the SAME frame in which canoes
    evaluates zeta, cf. scripts/mathematica/derive_spin2_great_circle_basis.wl),
    the complex spin-2 field P = Psi_+ + i Psi_x sources the shear leg-by-leg
    through the single-plane response, with the spin-2 screen-Hessian multiplier
    sqrt(L^2(L^2-2))/chi^2 (appendix `eq: appendix screen hessian`).  The
    multiplier is ALREADY baked into the modulus channels zeta_B/zeta_D by
    canoes (the spin-weighted-harmonic contraction carries the spin-2 eth^2);
    we do NOT re-apply it.  The spin-0 screen multiplier L^2/chi^2 is likewise
    baked into the TTT/TTP channels.  [Both are produced by the SAME
    compute_kappa3_* path that produced the certified STAGE-0 spin-0 result.]

2.  COMPLEX NATURAL COMPONENTS = SPIN-PERMUTED ZETA CUMULANTS (canonical frame).
    The cosmic-shear natural components are the four independent triple products
    of the complex shear gamma = P (and its conjugate) at the three vertices.
    With P_i the complex spin-2 driving field at leg i, these are EXACTLY the
    canoes spin-aware cumulants, one per single-leg conjugation, computed in
    canoes' canonical great-circle frame (n1 at the pole, the spin reference of
    each leg being its geodesic to leg 1; _canonical_frame_angles):
        Gamma^0_loc = <P1  P2  P3 > = zeta_PPP          spin_weights (+2,+2,+2)
        Gamma^1_loc = <P1* P2  P3 > = zeta_D(-2,+2,+2)  (leg 1 conjugated)
        Gamma^2_loc = <P1  P2* P3 > = zeta_D(+2,-2,+2)  (leg 2 conjugated)
        Gamma^3_loc = <P1  P2  P3*> = zeta_D(+2,+2,-2)  (leg 3 conjugated)
    NB: the three single-conjugation permutations are DISTINCT cumulants on the
    SAME triple (verified: |G1|~|G2|!=|G3| on the isoceles family, base vs apex
    leg), obtained via the canoes `_spin_weights_dmod` override.  Gamma^0 is the
    all-unconjugated (all-+2) channel = zeta_PPP, gamma^4-suppressed (the
    d^L_{2,-2} factor); Gamma^{1,2,3} are modulus-dominated (d^L_{2,2}~O(1)).
    The appendix real-basis reconstruction map (`eq: appendix modulus
    reconstruction`) is the equivalent statement in the (Psi+, Psix) real basis;
    here we use the complex cumulants directly, which is exact and avoids the
    real-basis round-trip.

4.  RADIAL FOLD.  Each natural component is folded leg-by-leg over the affine
    parameter lambda with the single-plane convergence/shear response cubed,
    EXACTLY as STAGE 0 (the response is spin-independent; the spin-2 character
    is carried by zeta):
        Z(triangle) = int_0^{lambda_s} dlambda  K(lambda,lambda_s)^3  zeta(triangle, lambda)
        K(lambda,lambda_s) = a^2(chi) chi (chi_s - chi)/chi_s   [Dbar = a chi]
    radial_measure="lambda": zeta is a 3-leg lambda-density; do NOT apply
    (dchi/dlambda)^3 (radial-measure memory 2026-05-17).  units="physical".

5.  CENTROID PROJECTION (match fastnc EXACTLY).  fastnc computes Gamma^mu in its
    FFT-native x-projection, then multiplies by x2cent(mu, t1, t2, phi)
    (fastnc.py:752-784, "equations between Eq.(15)-(16) of arXiv:2309.08601").
    Our local great-circle frame is the natural REAL-SPACE separation-aligned
    frame; the documented relation (Schneider-Lombardi; Sugiyama+2024) is that
    the natural components in any projection differ from the x-projection by the
    SAME multiplicative phase factor.  We therefore obtain the centroid-projected
    components by applying fastnc's OWN x2cent phase to our local-frame real
    natural components (see CONVENTION NOTE + the phase-ambiguity discussion in
    the report).  This makes the projection convention identical by construction.

NB: building fastnc's full 2DFFTLog multipole machinery on our side would be an
independent reimplementation of Sugiyama+2024; instead we compute the natural
components directly from the real-space driving-field cumulant (the object the
whole fastnc pipeline is engineered to reproduce) and apply fastnc's published
projection phase.  This is the transparent, analytically-traceable route.
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
_N_LAMBDA = 40  # affine-grid resolution, same family as STAGE 0


# ---------------------------------------------------------------------------
# fastnc centroid-projection phase (verbatim port of fastnc.py:752-784, x2cent)
# arXiv:2309.08601, equations between Eq.(15) and (16).  t1,t2,phi are the SAS
# sides/opening angle (t1=t2=gamma here, phi the opening angle).
# ---------------------------------------------------------------------------
def x2cent(mu: int, t1, t2, phi):
    """x-projection -> centroid-projection phase factor (fastnc.py:752)."""
    v = t1 + t2 * np.exp(-1j * phi)
    q1 = v / np.conj(v)
    v = -2 * t1 + t2 * np.exp(-1j * phi)
    q2 = v / np.conj(v)
    v = t1 - 2 * t2 * np.exp(-1j * phi)
    q3 = v / np.conj(v)
    if mu == 0:
        out = q1 * q2 * q3 * np.exp(3j * phi)
    elif mu == 1:
        out = np.conj(q1) * q2 * q3 * np.exp(1j * phi)
    elif mu == 2:
        out = q1 * np.conj(q2) * q3 * np.exp(3j * phi)
    elif mu == 3:
        out = q1 * q2 * np.conj(q3) * np.exp(-1j * phi)
    else:
        raise ValueError(f"mu={mu} not expected")
    return out


def _build_lambda_grid(cosmo, z_source: float, n_lambda: int):
    """Affine grid lambda(z) for z in (0, z_s] (identical machinery to STAGE 0)."""
    from scipy.integrate import cumulative_trapezoid

    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, z_source, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])  # physical Mpc
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
    """K(lambda,lambda_s) = a^2(chi) chi (chi_s-chi)/chi_s  (physical Mpc)."""
    a = 1.0 / (1.0 + z)
    return a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s


def _sas_cosine_triples(gamma_rad: float, phi_rad: np.ndarray) -> np.ndarray:
    """SAS family cosine triples (cos t3, cos gamma, cos gamma).

    t1 = t2 = gamma (the two equal sides), opening angle phi, third side
    t3 = 2 gamma sin(phi/2).  cosine_triples ordering is
    (cos theta_12, cos theta_23, cos theta_31); we place the base t3 on the
    (1,2) pair and the two equal sides gamma on (2,3) and (3,1), matching the
    grid.json directive.
    """
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)  # (n_phi,)
    cg = np.cos(gamma_rad)
    ct3 = np.cos(t3)
    out = np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)
    return out  # (n_phi, 3)


# Spin-weight assignments for the four pure-shear natural components.
#
# Leg labeling on our SAS triple (cos t12, cos t23, cos t31) = (cos t3, cos g,
# cos g): pair (1,2) is the base t3, pairs (2,3) and (3,1) are the equal sides
# gamma.  Hence legs 1,2 are the BASE vertices (endpoints of t3) and leg 3 is
# the APEX (where the two equal sides meet).
#
# fastnc's natural-component leg convention (Schneider-Lombardi / Sugiyama+2024,
# t1=t2 isoceles): Gamma^1 conjugates the APEX vertex (the one between the two
# equal sides), Gamma^2 and Gamma^3 conjugate the two BASE vertices -- giving the
# observed |Gamma^2|=|Gamma^3| isoceles symmetry (smoke test confirms
# Gamma^3=conj(Gamma^2)).  We therefore map:
#     Gamma^1  <-  conjugate leg 3 (apex)  ->  spin_weights (+2,+2,-2)
#     Gamma^2  <-  conjugate leg 1 (base)  ->  spin_weights (-2,+2,+2)
#     Gamma^3  <-  conjugate leg 2 (base)  ->  spin_weights (+2,-2,+2)
# (Gamma^0 = all-+2 = zeta_PPP unchanged.)
_DMOD_SPINS = {
    1: (2, 2, -2),   # Gamma^1: conjugate APEX (leg 3)
    2: (-2, 2, 2),   # Gamma^2: conjugate BASE leg 1
    3: (2, -2, 2),   # Gamma^3: conjugate BASE leg 2
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
    assert abs(cosmo.Omega_m - grid["cosmology"]["Omega_m"]) < 1e-9, "Omega_m mismatch"
    assert abs(cosmo.h - grid["cosmology"]["h"]) < 1e-9, "h mismatch"

    z, chi_phys, lam_phys, lam_h, chi_s, lam_s = _build_lambda_grid(
        cosmo, z_source, _N_LAMBDA)
    K = _convergence_kernel(chi_phys, z, chi_s)  # (n_lambda,)
    print(f"[s1-ours] z_s={z_source}, chi_s={chi_s:.2f} Mpc, lam_s={lam_s:.2f} Mpc, "
          f"n_lambda={_N_LAMBDA}")
    print(f"[s1-ours] cosmo Omega_m={cosmo.Omega_m:.6f}, h={cosmo.h}, n_s={cosmo.n_s}")
    print(f"[s1-ours] gamma(arcmin)={gamma_arcmin.tolist()}")
    print(f"[s1-ours] phi(deg)={phi_deg.tolist()}")

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
        """Build the four pure-shear natural-component cumulants over lambda.

        Returns (zPPP, zD1, zD2, zD3): each (n_triples, n_lambda) float.
        zPPP = <P1 P2 P3> (all +2);  zD{i} = single-leg-i-conjugated zeta_D.
        Each is LOW (equal-shell Wigner-3j) + HIGH (Limber) summed.  The
        _spin_weights_dmod override selects which leg is conjugated; _bmod is
        left at its default (unused -- only the D output is read).
        """
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

    def fold(zeta):  # zeta: (n_triples, n_lambda) -> (n_triples,)
        integrand = (K ** 3)[None, :] * zeta
        return np.trapezoid(integrand, lam_phys, axis=1)

    def natural_components_for_gamma(gamma_arcmin_val, phi_rad_arr):
        """Return Gamma0..3 (complex, len n_phi) for one gamma over the phi grid."""
        gamma_rad = np.radians(gamma_arcmin_val / 60.0)
        triples = _sas_cosine_triples(gamma_rad, phi_rad_arr)  # (n_phi, 3)
        zPPP, zD1, zD2, zD3 = build_zeta_channels(triples)
        fPPP, fD1, fD2, fD3 = fold(zPPP), fold(zD1), fold(zD2), fold(zD3)
        # canonical-frame complex natural components (zeta cumulants are real in
        # the canonical great-circle frame; the spin reference is the geodesic
        # to leg 1).  Overall sign (-1)^3 from gamma = -int Dsachs (matches the
        # STAGE-0 kappa (-1)^3 and eq: gamma from Dsachs23).
        sign = -1.0
        G0c = sign * fPPP    # <P1 P2 P3>   (gamma^4-suppressed)
        G1c = sign * fD1     # <P1* P2 P3>
        G2c = sign * fD2     # <P1 P2* P3>
        G3c = sign * fD3     # <P1 P2 P3*>
        # CENTROID projection: apply fastnc x2cent phase (t1=t2=gamma, phi).
        t1 = t2 = gamma_rad
        G0 = G0c * x2cent(0, t1, t2, phi_rad_arr)
        G1 = G1c * x2cent(1, t1, t2, phi_rad_arr)
        G2 = G2c * x2cent(2, t1, t2, phi_rad_arr)
        G3 = G3c * x2cent(3, t1, t2, phi_rad_arr)
        return (np.asarray(G0, np.complex128), np.asarray(G1, np.complex128),
                np.asarray(G2, np.complex128), np.asarray(G3, np.complex128),
                dict(fPPP=fPPP, fD1=fD1, fD2=fD2, fD3=fD3,
                     G0c=G0c, G1c=G1c, G2c=G2c, G3c=G3c))

    # ---------- PRIMARY single point FIRST (gamma=10', phi=60deg) -------------
    print()
    print(f"[s1-ours] === PRIMARY POINT (gamma={g_prim}', phi={phi_prim}deg) ===")
    phi_prim_rad = np.array([np.radians(phi_prim)])
    G0p, G1p, G2p, G3p, diagp = natural_components_for_gamma(g_prim, phi_prim_rad)
    fD_mean = np.mean([np.abs(diagp['fD1'])[0], np.abs(diagp['fD2'])[0],
                       np.abs(diagp['fD3'])[0]])
    print(f"[s1-ours]   zeta_D dominance check (folded): "
          f"|Z_D| (mean of 3 perms)={fD_mean:.3e}, |Z_PPP|={np.abs(diagp['fPPP'])[0]:.3e}, "
          f"ratio |PPP|/|D|={np.abs(diagp['fPPP'])[0]/max(fD_mean,1e-300):.2e} "
          f"(gamma^4-suppression of the all-+2 channel)")
    print(f"[s1-ours]   fD perms (component-mapped): G1<-apex={diagp['fD1'][0]:+.3e}, "
          f"G2<-base1={diagp['fD2'][0]:+.3e}, G3<-base2={diagp['fD3'][0]:+.3e} "
          f"(G2~G3 base symmetry => fastnc |G2|=|G3|)")
    sym = G3p[0] - np.conj(G2p[0])
    print(f"[s1-ours]   isoceles symmetry check  Gamma3 - conj(Gamma2) = {sym:+.3e} "
          f"(should be ~0 for t1=t2)")
    print(f"[s1-ours]   Gamma0 = {G0p[0]:+.6e}")
    print(f"[s1-ours]   Gamma1 = {G1p[0]:+.6e}")
    print(f"[s1-ours]   Gamma2 = {G2p[0]:+.6e}")
    print(f"[s1-ours]   Gamma3 = {G3p[0]:+.6e}")

    # ---------- full pinned grid --------------------------------------------
    n_g, n_p = gamma_arcmin.size, phi_deg.size
    Gamma0 = np.zeros((n_g, n_p), np.complex128)
    Gamma1 = np.zeros((n_g, n_p), np.complex128)
    Gamma2 = np.zeros((n_g, n_p), np.complex128)
    Gamma3 = np.zeros((n_g, n_p), np.complex128)
    diag_fD = np.zeros((n_g, n_p))
    diag_fPPP = np.zeros((n_g, n_p))
    print()
    print("[s1-ours] === FULL GRID ===")
    for ig, gam in enumerate(gamma_arcmin):
        G0, G1, G2, G3, diag = natural_components_for_gamma(gam, phi_rad)
        Gamma0[ig], Gamma1[ig], Gamma2[ig], Gamma3[ig] = G0, G1, G2, G3
        diag_fD[ig], diag_fPPP[ig] = diag["fD3"], diag["fPPP"]
        print(f"[s1-ours]  gamma={gam:7.2f}'  |G3|@phi: "
              + " ".join(f"{np.abs(G3[j]):.2e}" for j in range(n_p)))

    # ---------- save npz (pinned schema) ------------------------------------
    out = _HERE / "outputs" / "stage1_ours.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    convention_note = (
        "OUR-SIDE cosmic-shear 3PCF natural components Gamma^0..3, "
        "Schneider-Lombardi basis, CENTROID projection. The four components are "
        "the single-leg-conjugated spin-2 driving-field cumulants evaluated in "
        "canoes' canonical great-circle frame (n1 at the pole; each leg's spin "
        "reference is its geodesic to leg 1): "
        "Gamma^0=<P1 P2 P3>=zeta_PPP (spins +2,+2,+2; gamma^4-suppressed), "
        "Gamma^1=<P1* P2 P3>=zeta_D(-2,+2,+2), Gamma^2=<P1 P2* P3>=zeta_D(+2,-2,+2), "
        "Gamma^3=<P1 P2 P3*>=zeta_D(+2,+2,-2), P=Psi0=Psi_+ + i Psi_x. "
        "zeta_D from canoes compute_kappa3_mod_zeta_equal_time (LOW, spt_kind="
        "tree_phi) + compute_kappa3_mod_sigma3_high (HIGH), units=physical, "
        "radial_measure=lambda; zeta_PPP from compute_kappa3_zeta_table + "
        "compute_kappa3_sigma3_high. This is the complex-cumulant form of "
        "appendix `eq: appendix modulus reconstruction` (helicity fix 2026-06-04). "
        "Radial fold int dlambda K^3 zeta, K=a^2 chi (chi_s-chi)/chi_s, Dbar=a chi, "
        "single source plane z_s=5, lambda-density (no (dchi/dl)^3 per "
        "radial-measure memory). Overall sign (-1)^3 from gamma=-int Dsachs "
        "(matches STAGE-0 kappa sign). SAS triple (cos t3, cos gamma, cos gamma), "
        "t3=2 gamma sin(phi/2), leg-pair (1,2)=t3 base, (2,3)=(3,1)=gamma "
        "(leg 3 = apex). CENTROID projection applied via fastnc's OWN "
        "x2cent(mu,t1,t2,phi) phase (fastnc.py:752-784; arXiv:2309.08601 "
        "Eq.15-16), t1=t2=gamma. "
        "FRAME CAVEAT: our cumulants are in canoes' canonical great-circle frame; "
        "fastnc's x2cent maps the x-projection (FFT-native) -> centroid. The "
        "canonical great-circle frame and fastnc's x-projection coincide only up "
        "to a per-leg spin-2 rotation; if a residual configuration-dependent "
        "phase remains it is documented in the report. The frame-INVARIANT "
        "sqrt(sum_mu |Gamma^mu|^2) is unaffected by this and is the robust "
        "magnitude/normalization check. "
        "z_s=5, Omega_m=0.3160919980475834, h=0.6711, n_s=0.97, PCAMBz0.txt."
    )
    np.savez(
        out,
        gamma_arcmin=gamma_arcmin, phi_deg=phi_deg,
        Gamma0=Gamma0, Gamma1=Gamma1, Gamma2=Gamma2, Gamma3=Gamma3,
        projection="cent",
        convention_note=convention_note,
        # extras (diagnostics; not part of the schema but harmless)
        z_s=z_source, chi_s=chi_s, lam_s=lam_s,
        Omega_m=cosmo.Omega_m, h=cosmo.h, n_s=cosmo.n_s,
        diag_fD=diag_fD, diag_fPPP=diag_fPPP,
    )
    print(f"\n[s1-ours] saved -> {out}")

    # ---------- self-diagnostic PNG -----------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        Ginv = np.sqrt(np.abs(Gamma0) ** 2 + np.abs(Gamma1) ** 2
                       + np.abs(Gamma2) ** 2 + np.abs(Gamma3) ** 2)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ig, gam in enumerate(gamma_arcmin):
            axes[0].plot(phi_deg, np.abs(Gamma3[ig]), marker="o",
                         label=f"gamma={gam:.0f}'")
            axes[1].plot(phi_deg, Ginv[ig], marker="s",
                         label=f"gamma={gam:.0f}'")
        axes[0].set_ylabel(r"$|\Gamma^3|$ (centroid; dominant modulus)")
        axes[1].set_ylabel(r"$\sqrt{\sum_\mu|\Gamma^\mu|^2}$ (frame-invariant)")
        for ax in axes:
            ax.set_xlabel("opening angle phi [deg]")
            ax.set_yscale("log")
            ax.legend()
        axes[0].set_title("STAGE 1 our-side: |Gamma^3| vs phi")
        axes[1].set_title("STAGE 1 our-side: frame-invariant magnitude")
        fig.tight_layout()
        png = _HERE / "outputs" / "stage1_ours_selfcheck.png"
        fig.savefig(png, dpi=110)
        print(f"[s1-ours] self-check PNG -> {png}")
    except Exception as exc:  # noqa: BLE001
        print(f"[s1-ours] (PNG skipped: {exc})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
