"""STAGE 1 v2 (ours side) -- FRESH REBUILD with the FIXED canoes kernel.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

This is a byte-for-byte copy of stage1_ours_v2.py's channel build + natural
component map (the DERIVED great-circle -> centroid spin-2 projection phase),
EXCEPT it does a genuine FRESH recompute of every channel (no v1 reuse) so the
output reflects the group-zeroing bug fix in
canoes/src/canoes/nuell/correlation/_spin_aware_three_pt.py
(branch fix/spin2-zetaD-2_2_-2-norm).

The bug (pinned by scripts/probe_zetaD_kernel_oracle.py): the (l2,l3)-keyed
group-batched m-sum reduction in SpinAwareZetaEvaluator._build_kernel_numpy read
count_g / the W3j reshape from the group-FIRST triple. For s1=+/-2 channels the
first triple (smallest l1) has l1 < |s1|, so m_count=0 zeroed the ENTIRE group,
dropping ~half the valid (l1>=|s1|) cells -> zeta_D came out ~10^4x too small,
which (via |Gamma|=|fold(zeta_D)|) made the shear natural components 4-17x too
SMALL vs the reference -- equivalently the reference/canoes ratio was 4-17x.

Output: outputs/stage1_ours_v2_FIXED.npz (same schema as stage1_ours_v2.npz).
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


def _amb(phi: np.ndarray):
    half = phi / 2.0
    d1 = -np.arctan2(1.0, 3.0 * np.tan(half))
    d2 = -2.0 * np.pi + np.arctan2(1.0, 3.0 * np.tan(half))
    d3 = -0.5 * np.pi + np.arctan2(1.0, np.tan(half))
    return d1, d2, d3


def derived_phase(mu: int, phi: np.ndarray) -> np.ndarray:
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

    (1+z)^(2n-2) FIX (2026-06-09; fix/spin2-zetaD-2_2_-2-norm).  Removes the
    per-shell (1+z)^(2 n_legs - 2) over-count of the single-plane K^n fold
    against canoes' per-leg (1+z)^4 Sachs response (CCL-native-anchored in
    b2_kfold_2pt_vs_ccl.py: 2-leg over-count = (1+z)^2 exactly; pointwise /(1+z)^2
    restores Cl(hand)/Cl(CCL)=0.99).  Shear 3PCF = 3 legs -> (1+z)^-4.  Born
    E=E0/a per-leg reduction; canoes zeta + kernel K unchanged; phase-independent.
    """
    return (1.0 + z) ** (-(2 * n_legs - 2))


def _sas_cosine_triples(gamma_rad: float, phi_rad: np.ndarray) -> np.ndarray:
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad)
    ct3 = np.cos(t3)
    return np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)


_DMOD_SPINS = {
    1: (-2, 2, 2),
    2: (2, -2, 2),
    3: (2, 2, -2),
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
    print(f"[FIXED] z_s={z_source}, chi_s={chi_s:.2f} Mpc, lam_s={lam_s:.2f} Mpc")

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
        sign = -1.0
        G0 = sign * fPPP * derived_phase(0, phi_rad_arr)
        G1 = sign * fD1 * derived_phase(1, phi_rad_arr)
        G2 = sign * fD2 * derived_phase(2, phi_rad_arr)
        G3 = sign * fD3 * derived_phase(3, phi_rad_arr)
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
        print(f"[FIXED]  gamma={gam:7.2f}' done  |G3(phi120)|={abs(G3[-1]):.4e}")

    out = _HERE / "outputs" / "stage1_ours_v2_FIXED.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    convention_note = (
        "FRESH REBUILD of stage1_ours_v2 with the FIXED SpinAwareZetaEvaluator "
        "group-batched reduction (fix/spin2-zetaD-2_2_-2-norm). Same channel "
        "assignment + DERIVED great-circle->centroid spin-2 projection phase as "
        "stage1_ours_v2.py. zeta_D channels rebuilt via "
        "compute_kappa3_mod_zeta_equal_time/_sigma3_high; PPP via "
        "compute_kappa3_zeta_table/_sigma3_high. units=physical, "
        "radial_measure=lambda, fold int dlambda K^3 (1+z)^-4 zeta, "
        "K=a^2 chi(chi_s-chi)/chi_s, z_s=5. The (1+z)^-4 Born per-leg reduction "
        "(fix/spin2-zetaD-2_2_-2-norm; CCL-native-anchored) removes the single-"
        "plane K^3 fold's per-shell (1+z)^4 over-count; phase-independent so the "
        "relative Gamma^mu structure is unchanged."
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
    print(f"\n[FIXED] saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
