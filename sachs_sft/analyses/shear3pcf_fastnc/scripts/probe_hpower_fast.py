#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""FAST h-power check (HIGH-only, no slow LOW exact-Wigner build).

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Confirms SURPLUS (B): the 3-pt HIGH zeta_TTT carries h^6 (kappa3.py:3642), while
the dimensionally-equivalent equal-shell 3-leg lambda-density should carry h^4
(the SAME power the LOW direct-sigma3 paths apply at kappa3.py:2320/3151/3344 and
the 2-pt C_ell applies at kappa2.py:364, both PyCCL/GL-validated).

We build the HIGH zeta in BOTH units and read the empirical physical/h ratio.
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


def main() -> int:
    cosmo = FiducialCosmology()
    h = cosmo.h
    pk = load_pk_delta(n_s=cosmo.n_s)
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
        n_ell=96, n_phi=64, channels=("TTT",),
        lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = np.asarray(compute_kappa3_sigma3_high(triples, lam_h, units="physical", **common).zeta_TTT)
        zh = np.asarray(compute_kappa3_sigma3_high(triples, lam_h, units="h", **common).zeta_TTT)
    msk = np.isfinite(zp) & np.isfinite(zh) & (np.abs(zh) > 0)
    ratio = zp[msk] / zh[msk]
    med = float(np.median(ratio))
    hpow = np.log(np.median(np.abs(ratio))) / np.log(h)
    print(f"HIGH zeta_TTT physical/h median = {med:.6e}")
    print(f"  h^6 = {h**6:.6e},  h^4 = {h**4:.6e}")
    print(f"  empirical h-power = {hpow:.4f}  (code:kappa3.py:3642 applies h^6)")
    print(f"  2-pt C_ell + LOW sigma3 apply h^4 (kappa2.py:364, kappa3.py:2320/3151/3344)")
    print(f"  => SURPLUS h-power = {hpow - 4.0:.4f}  (expect +2.0 -> extra h^2)")
    np.savez(_OUT / "hpower_fast.npz", physical_over_h_median=med, hpow=hpow,
             h=h, note="HIGH zeta h-power; surplus h^2 vs validated h^4.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
