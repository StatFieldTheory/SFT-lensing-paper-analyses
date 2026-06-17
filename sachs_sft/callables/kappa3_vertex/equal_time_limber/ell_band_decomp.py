"""ℓ-band decomposition of the equal-shell driving-field cumulant ζ.

The part of ζ that enters the lensing 2PCF is the squeezed / collapsed 3PCF —
the near-collinear geometry (cos12,cos23,cos31)=(1,cosγ,cosγ). This script
decomposes ζ at that geometry into ℓ-band contributions, for several angular
scales γ, to quantify how much each ℓ region contributes (expectation: high-ℓ
small-scale power dominates ζ even at large γ).

Method (calls the canoes branch functions DIRECTLY on a minimal squeezed-triple
set — fast, exact, no callable/interpolation, no grid-matching):
  * LOW  (exact Wigner-3j, ℓ≤ell_max):  compute_kappa3_zeta_table at an ell_max
    ladder → cumulative → DIFFERENCE for the LOW bands (no lower bound exists).
  * HIGH (flat-sky Born-Limber, (lo,hi]): compute_kappa3_sigma3_high with a
    disjoint window per band → DIRECT (no differencing).
  * gate: Σ bands ≈ full (LOW(≤60)+HIGH(60,1000]) at the SAME triples.

Setup mirrors build_equal_time_limber_table.py (same cosmo, pk, λ-grid in
h-units, units=physical, measure=lambda). Run with the sft-wick env:
  LOKY_MAX_CPU_COUNT=12 CANOES_SUPPRESS_METAL_WARNING=1 PYTHONUNBUFFERED=1 \
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python ell_band_decomp.py [--smoke]
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

_TABLE = _HERE / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"

# ℓ-band edges. LOW (exact Wigner, ℓ≤60) via the ell_max ladder; HIGH (flat-sky,
# (60,1000]) via disjoint windows.
_LOW_LADDER = [30, 60]                       # cumulative ell_max -> bands [≤30],(30,60]
_HIGH_WINDOWS = [(60, 125), (125, 250), (250, 500), (500, 1000)]
_CHANNELS = ("TTT", "TTP", "TPP", "PPP")


def _squeezed_triples(gamma_rad: np.ndarray) -> np.ndarray:
    """(1, cosγ, cosγ) — two coincident directions separated by γ from the third."""
    cg = np.cos(gamma_rad)
    return np.stack([np.ones_like(cg), cg, cg], axis=1)  # (n_gamma, 3)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true",
                    help="3 γ + 1 HIGH window — wiring/timing/sign check.")
    ap.add_argument("--n-gamma", type=int, default=28)
    ap.add_argument("--out", type=Path, default=_HERE / "ell_band_decomp_results.npz")
    args = ap.parse_args()

    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
    )
    # Conjugate-helicity modulus channels (not in canoes.sachs.__all__):
    from canoes.sachs.kappa3 import (  # noqa: E402
        compute_kappa3_mod_sigma3_high,
        compute_kappa3_mod_zeta_equal_time,
    )

    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    lam_phys = np.asarray(np.load(_TABLE, allow_pickle=False)["lambda_shells_Mpc"],
                          dtype=np.float64)
    lam_h = lam_phys * cosmo.h

    # γ grid: log-spaced 0.5' .. 5000' (the FK config's range).
    if args.smoke:
        gamma_arcmin = np.array([0.5, 50.0, 1000.0])
        high_windows = [(60, 1000)]
        low_ladder = [60]
    else:
        gamma_arcmin = np.geomspace(0.5, 5000.0, int(args.n_gamma))
        high_windows = _HIGH_WINDOWS
        low_ladder = _LOW_LADDER
    gamma_rad = np.radians(gamma_arcmin / 60.0)
    triples = _squeezed_triples(gamma_rad)
    ng, nl = triples.shape[0], lam_phys.size
    print(f"[band] cosmo Ω_m={cosmo.Omega_m} h={cosmo.h} n_s={cosmo.n_s}")
    print(f"[band] {ng} squeezed γ-triples × {nl} λ-shells "
          f"[{lam_phys[0]:.0f}..{lam_phys[-1]:.0f}] Mpc")

    common_low = dict(pk_matter_today=pk, cosmo=cosmo, Nmax=128, spt_kind="tree_phi",
                      units="physical", lambda_convention="project",
                      radial_discretization="sample", radial_measure="lambda",
                      channels=_CHANNELS, parallelism="sequential", shell_block_size=4,
                      triple_chunk_size=4096, n_workers=1, fuse_tree_phi_channels=True)
    common_high = dict(pk_matter_today=pk, cosmo=cosmo, ell_min=1e-3, n_ell=96,
                       n_phi=64, channels=_CHANNELS, units="physical",
                       lambda_convention="project", radial_discretization="sample",
                       radial_measure="lambda")

    def _chan_dict(out) -> dict:
        return {ch: np.asarray(getattr(out, f"zeta_{ch}"), dtype=np.float64)
                for ch in _CHANNELS}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        # LOW cumulative ladder.
        low_cum = {}
        for emax in low_ladder:
            t0 = time.perf_counter()
            out = compute_kappa3_zeta_table(triples, lam_h, ell_max=int(emax), **common_low)
            low_cum[emax] = _chan_dict(out)
            print(f"[band] LOW ell_max={emax}: {time.perf_counter()-t0:.1f}s")
        # HIGH disjoint windows.
        high_bands = {}
        for (lo, hi) in high_windows:
            t0 = time.perf_counter()
            out = compute_kappa3_sigma3_high(triples, lam_h, ell_cut=int(lo),
                                             ell_high_max=int(hi), **common_high)
            high_bands[(lo, hi)] = _chan_dict(out)
            print(f"[band] HIGH ({lo},{hi}]: {time.perf_counter()-t0:.1f}s")
        # Full reference HIGH (60, max) for the Σ-bands gate.
        hi_max = high_windows[-1][1]
        t0 = time.perf_counter()
        out = compute_kappa3_sigma3_high(triples, lam_h, ell_cut=60,
                                         ell_high_max=int(hi_max), **common_high)
        high_full = _chan_dict(out)
        print(f"[band] HIGH full (60,{hi_max}]: {time.perf_counter()-t0:.1f}s")

        # Conjugate-helicity MODULUS channels: Bmod=<Φ00 Ψ0 conj(Ψ0)> (spin 0,2,-2),
        # Dmod=<Ψ0 Ψ0 conj(Ψ0)> (spin 2,2,-2). Same triples/λ/cosmo/Pδ as above;
        # LOW(≤max_low) + HIGH(60,hi_max] to match the canonical full_{ch}.
        t0 = time.perf_counter()
        bmod_low, dmod_low = compute_kappa3_mod_zeta_equal_time(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_max=int(low_ladder[-1]), Nmax=128, spt_kind="tree_phi",
            units="physical", lambda_convention="project", radial_measure="lambda")
        print(f"[band] MOD LOW ell_max={low_ladder[-1]}: {time.perf_counter()-t0:.1f}s")
        t0 = time.perf_counter()
        bmod_high, dmod_high = compute_kappa3_mod_sigma3_high(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_cut=60, ell_high_max=int(hi_max), ell_min=1e-3, n_ell=96, n_phi=64,
            units="physical", lambda_convention="project", radial_measure="lambda")
        print(f"[band] MOD HIGH (60,{hi_max}]: {time.perf_counter()-t0:.1f}s")
        full_Bmod = np.asarray(bmod_low + bmod_high, dtype=np.float64)
        full_Dmod = np.asarray(dmod_low + dmod_high, dtype=np.float64)

    # Assemble LOW bands by differencing the cumulative ladder.
    low_edges = [0] + list(low_ladder)
    bands = []  # list of (label, lo, hi, {ch: (ng,nl)})
    prev = None
    for i, emax in enumerate(low_ladder):
        lo = low_edges[i]
        d = {ch: low_cum[emax][ch] - (prev[ch] if prev is not None else 0.0)
             for ch in _CHANNELS}
        bands.append((f"LOW ({lo},{emax}]" if lo else f"LOW [2,{emax}]", lo, emax, d))
        prev = low_cum[emax]
    for (lo, hi) in high_windows:
        bands.append((f"HIGH ({lo},{hi}]", lo, hi, high_bands[(lo, hi)]))

    # Gate: Σ bands vs full = LOW(≤max_low) + HIGH(60,hi_max].
    full = {ch: low_cum[low_ladder[-1]][ch] + high_full[ch] for ch in _CHANNELS}
    band_sum = {ch: sum(b[3][ch] for b in bands) for ch in _CHANNELS}
    print("\n[GATE] Σ-bands vs full (max |rel diff| over γ,λ where |full|>0):")
    for ch in _CHANNELS:
        f = full[ch]; s = band_sum[ch]
        m = np.abs(f) > 0
        rel = np.abs(s[m] - f[m]) / np.abs(f[m])
        print(f"   {ch}: max rel diff = {rel.max() if rel.size else 0.0:.2e} "
              f"(median {np.median(rel) if rel.size else 0.0:.2e})")

    # Save everything for the figure.
    save = dict(gamma_arcmin=gamma_arcmin, lambda_phys=lam_phys,
                z_shells=np.asarray(np.load(_TABLE)["z_shells"]),
                band_labels=np.array([b[0] for b in bands]),
                band_lo=np.array([b[1] for b in bands], float),
                band_hi=np.array([b[2] for b in bands], float))
    for ch in _CHANNELS:
        save[f"full_{ch}"] = full[ch]
        for j, b in enumerate(bands):
            save[f"band{j}_{ch}"] = b[3][ch]
    # Modulus channels (Fig 12 shows these in place of the γ^4-suppressed TPP/PPP).
    save["full_Bmod"] = full_Bmod
    save["full_Dmod"] = full_Dmod
    np.savez(args.out, **save)
    print(f"[band] saved -> {args.out}")

    # Quick TTT table at the far source shell (z≈5.7, dominant).
    il = nl - 1
    print(f"\n[TTT] per-band ζ at λ={lam_phys[il]:.0f} Mpc (z≈{save['z_shells'][il]:.2f}):")
    hdr = "  γ[arcmin] " + " ".join(f"{b[0]:>16}" for b in bands) + f" {'FULL':>13}"
    print(hdr)
    for ig in range(ng):
        row = f"  {gamma_arcmin[ig]:9.2f}"
        for b in bands:
            row += f" {b[3]['TTT'][ig, il]:+16.3e}"
        row += f" {full['TTT'][ig, il]:+13.3e}"
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
