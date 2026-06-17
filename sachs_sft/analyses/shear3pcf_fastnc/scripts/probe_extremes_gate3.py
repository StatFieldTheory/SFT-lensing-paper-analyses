"""Gate-3 extremes probe: gamma=150', phi=120deg, squeezed -- finite, no blow-up.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Directly exercises the FIXED canoes kappa3 emitters (LOW zeta_table + HIGH
sigma3_high + the three spin-2 modulus D-channels) at the EXTREME corners of the
shear-3PCF SAS family, following the boundary-validation methodology (include
extreme parameter values + endpoints + the genuinely squeezed degenerate limit):

  * gamma = 150 arcmin (the LARGE-separation endpoint of the stage1 grid)
  * phi   = 120 deg     (the WIDE-opening endpoint of the stage1 grid)
  * squeezed: phi -> 0  (t3 -> 0, cosine-triple -> (1, cos g, cos g)) -- the
    DEGENERATE-triangle collapse where the spin-2 phases are delicate.
  * an explicit phi sweep across [eps, 120deg] at gamma=150' to confirm no
    blow-up anywhere on the closure curve.

For each config and each of the 4 emitter channels (PPP + the three (2,2,-2)-
family D-channels Dmod1/2/3) we assert:
  - np.all(np.isfinite(zeta))             (no NaN / Inf)
  - max|zeta| < 1e6  AND  no catastrophic blow-up relative to the gamma=10'
    primary scale (rel growth < 1e6x) -- the boundary-validation "no
    catastrophic blowup" guard.

This is a pure emitter-level finiteness/no-blowup probe (NOT a fold or a
reference compare); it certifies the FIXED group-zeroing + h^4 source does not
introduce NaN/Inf or a runaway at the extreme corners.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.cosmo.background import chi as chi_of_z  # noqa: E402

_HERE = Path(__file__).resolve().parent.parent
_N_LAMBDA = 40
_Z_SOURCE = 5.0
_DMOD_SPINS = {1: (-2, 2, 2), 2: (2, -2, 2), 3: (2, 2, -2)}

# Finiteness / no-blowup thresholds (boundary-validation guards).
_ABS_CEILING = 1.0e6          # max|zeta| absolute ceiling
_REL_BLOWUP = 1.0e6           # max growth vs the gamma=10' primary reference scale


def _build_lambda_h(cosmo, z_source, n_lambda):
    from scipy.integrate import cumulative_trapezoid

    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, z_source, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_source / n_lambda, z_source, n_lambda)
    lam_phys = np.interp(z, zf, lam_f)
    return lam_phys * h


def _sas_cosine_triples(gamma_rad, phi_rad):
    phi_rad = np.atleast_1d(np.asarray(phi_rad, dtype=np.float64))
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad)
    ct3 = np.cos(t3)
    return np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)


def _build_channels(triples, lam_h, cosmo, pk):
    from canoes.sachs import (
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    from canoes.sachs.kappa3 import (
        compute_kappa3_mod_sigma3_high,
        compute_kappa3_mod_zeta_equal_time,
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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low = compute_kappa3_zeta_table(triples, lam_h, **common_low)
        high = compute_kappa3_sigma3_high(triples, lam_h, **common_high)
        comb = kappa3_combine_low_high(low, high)
        out = {"PPP": np.asarray(comb.zeta_PPP, dtype=np.float64)}
        for mu, sw in _DMOD_SPINS.items():
            dlo = compute_kappa3_mod_zeta_equal_time(
                triples, lam_h, _spin_weights_dmod=sw, **mod_low)[1]
            dhi = compute_kappa3_mod_sigma3_high(
                triples, lam_h, _spin_weights_dmod=sw, **mod_high)[1]
            out[f"Dmod{mu}"] = np.asarray(dlo + dhi, dtype=np.float64)
    return out


def main() -> int:
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    lam_h = _build_lambda_h(cosmo, _Z_SOURCE, _N_LAMBDA)

    arcmin2rad = np.pi / 10800.0

    # Reference scale: the gamma=10' primary, phi=60 (well-posed) -- used only to
    # define "no catastrophic blow-up" relative growth.
    g10 = 10.0 * arcmin2rad
    ref_triples = _sas_cosine_triples(g10, np.radians([60.0]))
    ref_ch = _build_channels(ref_triples, lam_h, cosmo, pk)
    ref_scale = {k: float(np.max(np.abs(v))) for k, v in ref_ch.items()}
    print("[ref scale gamma=10' phi=60deg] " +
          "  ".join(f"{k}={ref_scale[k]:.3e}" for k in ref_ch))

    # Extreme configs (the gate-3 corners).
    g150 = 150.0 * arcmin2rad
    EPS_PHI_DEG = 0.5  # squeezed surrogate (do NOT evaluate exactly at phi=0)
    configs = [
        ("gamma=150' phi=120deg (WIDE endpoint)", g150, [120.0]),
        ("gamma=150' phi=0.5deg (SQUEEZED -> degenerate)", g150, [EPS_PHI_DEG]),
        ("gamma=150' phi SWEEP [0.5..120]deg", g150,
         [0.5, 5.0, 20.0, 45.0, 60.0, 90.0, 120.0]),
        ("gamma=10' phi=120deg (small-gamma wide)", g10, [120.0]),
    ]

    all_ok = True
    rows = []
    for label, gamma_rad, phi_deg_list in configs:
        phi_rad = np.radians(phi_deg_list)
        triples = _sas_cosine_triples(gamma_rad, phi_rad)
        ch = _build_channels(triples, lam_h, cosmo, pk)
        print("\n" + "=" * 90)
        print(label)
        cos_t3 = triples[:, 0]
        print(f"  phi_deg={phi_deg_list}")
        print(f"  cos(t3) range=[{cos_t3.min():.6f},{cos_t3.max():.6f}]  "
              f"(t3->0 squeezed gives cos->1)")
        for name, arr in ch.items():
            finite = bool(np.all(np.isfinite(arr)))
            amax = float(np.max(np.abs(arr))) if arr.size else 0.0
            rel = amax / max(ref_scale.get(name, 1e-300), 1e-300)
            ok = finite and (amax < _ABS_CEILING) and (rel < _REL_BLOWUP)
            all_ok = all_ok and ok
            flag = "OK " if ok else "BAD"
            print(f"    {flag} {name:7s}: finite={finite}  max|zeta|={amax:.4e}  "
                  f"rel_vs_ref={rel:.3e}")
            rows.append((label, name, finite, amax, rel, ok))

    print("\n" + "=" * 90)
    print(f"EXTREMES VERDICT: {'ALL FINITE, NO BLOW-UP' if all_ok else 'FAILURE DETECTED'}")
    print(f"  ABS ceiling={_ABS_CEILING:.0e}  rel-blowup guard={_REL_BLOWUP:.0e}x")

    out = _HERE / "outputs" / "extremes_gate3_probe.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out,
        labels=np.array([r[0] for r in rows]),
        channels=np.array([r[1] for r in rows]),
        finite=np.array([r[2] for r in rows]),
        max_abs=np.array([r[3] for r in rows]),
        rel_vs_ref=np.array([r[4] for r in rows]),
        ok=np.array([r[5] for r in rows]),
        all_ok=np.array(all_ok),
        ref_scale=np.array([ref_scale[k] for k in ("PPP", "Dmod1", "Dmod2", "Dmod3")]),
        note=("Gate-3 extremes probe: gamma=150', phi=120deg, squeezed phi->0. "
              "Direct FIXED-canoes emitter finiteness + no-catastrophic-blowup "
              "guards (abs<1e6, rel-vs-gamma10-ref<1e6x)."),
    )
    print(f"saved -> {out}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
