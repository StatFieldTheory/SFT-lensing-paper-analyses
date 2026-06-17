#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Measure the deployed LOW vs HIGH contribution to each channel's folded shear
3PCF at the PRIMARY point, to localise the spin-2 4-17x.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

We have excluded (xAct + probes): the per-leg spin-2 eth^2 weight (->1), the
radial normalisation (L-independent), and the W_can/B_kappa scalar shape (same
for all spins; TTT closes).  The remaining question is whether the deployed
zeta_D is HIGH-dominated (prompt claim: by ~8 orders) and whether the LOW vs
HIGH split DIFFERS between TTT (closes 0.85) and the D-channels (4-17x).

For the PRIMARY SAS triangle (gamma=10', phi=60deg) we build, separately:
  - zeta_TTT  LOW  (compute_kappa3_zeta_table, TTT)        and HIGH
  - zeta_PPP  LOW  (compute_kappa3_zeta_table, PPP)        and HIGH
  - zeta_D3   LOW  (compute_kappa3_mod_zeta_equal_time)    and HIGH
fold each (INT dlambda K^3 zeta) and print |LOW|, |HIGH|, |LOW+HIGH|, HIGH/LOW.

This is the missing measurement: if D is HIGH-dominated AND TTT is LOW-dominated,
the 4-17x lives in the HIGH path's spin-2 sampling (already excluded as a per-leg
factor -> must be the J_2 orientation-kernel sampling of W_can).  If instead D's
LOW and HIGH are comparable and the LOW (operator) is the large piece, the bug is
in compute_bispectra_deri_psi0.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs import (  # noqa: E402
    compute_kappa3_sigma3_high, compute_kappa3_zeta_table, kappa3_combine_low_high,
)
from canoes.sachs.kappa3 import (  # noqa: E402
    compute_kappa3_mod_zeta_equal_time, compute_kappa3_mod_sigma3_high,
)
from scipy.integrate import cumulative_trapezoid  # noqa: E402

_N_LAMBDA = 40


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5.0 / _N_LAMBDA, 5.0, _N_LAMBDA)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * cosmo.h
    chi_s = float(chif[-1])
    a = 1.0 / (1.0 + z)
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s

    gamma_rad = np.radians(10.0 / 60.0)
    phi_rad = np.radians(60.0)
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad)
    triples = np.array([[np.cos(t3), cg, cg]], dtype=np.float64)

    def fold(zeta):
        return float(np.trapezoid((K ** 3) * zeta[0], lam_phys))

    common_low = dict(pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128,
                      spt_kind="tree_phi", units="physical",
                      lambda_convention="project", radial_discretization="sample",
                      radial_measure="lambda", parallelism="sequential",
                      shell_block_size=4, triple_chunk_size=4096, n_workers=1,
                      fuse_tree_phi_channels=True)
    common_high = dict(pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000,
                       ell_min=1e-3, n_ell=96, n_phi=64, units="physical",
                       lambda_convention="project", radial_discretization="sample",
                       radial_measure="lambda")
    mod_low = dict(pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128,
                   spt_kind="tree_phi", units="physical",
                   lambda_convention="project", radial_measure="lambda")
    mod_high = dict(pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000,
                    ell_min=1e-3, n_ell=96, n_phi=64, units="physical",
                    lambda_convention="project", radial_measure="lambda")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # TTT
        ttt_low = compute_kappa3_zeta_table(triples, lam_h, channels=("TTT",), **common_low)
        ttt_high = compute_kappa3_sigma3_high(triples, lam_h, channels=("TTT",), **common_high)
        zTTT_low = np.asarray(ttt_low.zeta_TTT, dtype=np.float64)
        zTTT_high = np.asarray(ttt_high.zeta_TTT, dtype=np.float64)
        # PPP
        ppp_low = compute_kappa3_zeta_table(triples, lam_h, channels=("PPP",), **common_low)
        ppp_high = compute_kappa3_sigma3_high(triples, lam_h, channels=("PPP",), **common_high)
        zPPP_low = np.asarray(ppp_low.zeta_PPP, dtype=np.float64)
        zPPP_high = np.asarray(ppp_high.zeta_PPP, dtype=np.float64)
        # D3 (spins (2,2,-2))
        d_low = compute_kappa3_mod_zeta_equal_time(
            triples, lam_h, _spin_weights_dmod=(2, 2, -2), **mod_low)[1]
        d_high = compute_kappa3_mod_sigma3_high(
            triples, lam_h, _spin_weights_dmod=(2, 2, -2), **mod_high)[1]
        zD_low = np.asarray(d_low, dtype=np.float64)
        zD_high = np.asarray(d_high, dtype=np.float64)

    print("=" * 92)
    print("Deployed folded shear weight (INT dlambda K^3 zeta), PRIMARY (gamma=10',phi=60deg)")
    print("=" * 92)
    print(f"{'channel':>8s} {'|fold LOW|':>14s} {'|fold HIGH|':>14s} "
          f"{'|fold LOW+HIGH|':>16s} {'HIGH/LOW':>10s}")
    rows = {}
    for name, zlo, zhi in [("TTT", zTTT_low, zTTT_high),
                           ("PPP", zPPP_low, zPPP_high),
                           ("D3", zD_low, zD_high)]:
        flo = fold(zlo)
        fhi = fold(zhi)
        fcomb = fold(zlo + zhi)
        rows[name] = (flo, fhi, fcomb)
        hl = abs(fhi) / abs(flo) if flo != 0 else float("inf")
        print(f"{name:>8s} {abs(flo):14.4e} {abs(fhi):14.4e} {abs(fcomb):16.4e} {hl:10.3e}")
    print("=" * 92)
    print("Interpretation: if D3 HIGH/LOW >> 1 and TTT HIGH/LOW << 1, the spin-2 4-17x")
    print("lives in the HIGH (Limber/J_S kernel) path; the LOW operator is subdominant for D.")

    np.savez(OUT / "low_high_split.npz",
             channels=np.array(["TTT", "PPP", "D3"]),
             fold_low=np.array([rows[k][0] for k in ("TTT", "PPP", "D3")]),
             fold_high=np.array([rows[k][1] for k in ("TTT", "PPP", "D3")]),
             fold_comb=np.array([rows[k][2] for k in ("TTT", "PPP", "D3")]),
             note="Deployed LOW vs HIGH folded shear weight at primary SAS point.")
    print(f"saved -> {OUT / 'low_high_split.npz'}")


if __name__ == "__main__":
    main()
