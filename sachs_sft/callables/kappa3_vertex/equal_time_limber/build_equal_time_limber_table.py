"""Reproducible build of the equal_time_limber κ³ vertex table.

The SECOND κ³ implementation alongside R_contracted: the NON-R-contracted
equal-time Limber κ³ that is the source of the validated FK κκ baseline
(``+1.219e-4`` at γ=0.5' after the 2026-06-09 h^4 correction; was ``+3.086e-5``
pre-fix — see the ``fk_kk_baseline`` meta note below).  Composes two canoes API
calls (ports the archived
``analysis_3/build_kappa3_equal_time_limber_l2setup.py``):

  * LOW  : ``compute_kappa3_zeta_table``     — exact Wigner-3j, ℓ ≤ ell_cut (60)
  * HIGH : ``compute_kappa3_sigma3_high``    — flat-sky Born-Limber, ℓ ∈ (60, 1000]
  * combine: ``kappa3_combine_low_high``     — channel-by-channel sum

→ ONE combined ``Kappa3Output`` NPZ (zeta_TTT/TTP/TPP/PPP over
(cosine_triples, λ_shell, λ_shell, λ_shell)).

Contrast with R_contracted:
  * ``already_R_contracted = False`` — sft-wick applies the response R itself.
  * ``equal_time = True``, ``limber_high = True``.
  * indexed by COSINE-triples (cKDTree) + λ_shells, not (L,L,ℓ₃) + λ.

λ-grid dependency: reads the L2 λ-grid from an L2 NPZ so this κ³ shares the
identical source grid as R_contracted.  Point ``--l2-npz`` at the R_contracted
table (or the archived L2 NPZ fallback).  Build R_contracted FIRST.

MEASURE: ``--radial-measure lambda`` (the FK kk baseline is the LIMBER grid
equal_time (1+z)^4/leg lambda-density result: +1.219e-4 h^4-corrected,
was +3.08e-5 pre-fix).  The callable reads the density convention from the
NPZ metadata and FAILS LOUD on a mismatch.

Usage:
    python build_equal_time_limber_table.py --smoke
    python build_equal_time_limber_table.py --l2-npz <R_contracted table.npz>
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

_HERE = Path(__file__).resolve().parent
OUTPUT = _HERE / "equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz"

# Fallback L2 grid if --l2-npz is not given (the archived deployed L2 NPZ).
_ARCHIVED_L2 = (
    Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/_archive/"
         "canoes_pipeline/inputs/l2_callables_2026_05_26/"
         "windowed_squeezed_kappa3_lmax1000_z5_covgrid16_omega0316_h06711_physMpc.npz")
)


def _load_l2_lambda_grid_h(l2_npz: Path, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Read the L2 NPZ lambda_grid (physical Mpc) → (λ_phys, λ_h [Mpc/h])."""
    d = np.load(l2_npz, allow_pickle=False)
    if "lambda_grid" not in d:
        raise KeyError(f"{l2_npz} has no 'lambda_grid' (keys={list(d.keys())})")
    lam_phys = np.asarray(d["lambda_grid"], dtype=np.float64)
    if lam_phys.ndim != 1 or not np.all(np.diff(lam_phys) > 0):
        raise ValueError("L2 lambda_grid must be strictly-increasing 1-D")
    return lam_phys, lam_phys * float(h)


def _triangle_filtered_triples(n_cosine: int) -> np.ndarray:
    """Triangle-filtered (cos γ12, cos γ23, cos γ31) grid (mirrors the archive)."""
    c = np.linspace(-1.0, 1.0, n_cosine)
    c12, c23, c31 = np.meshgrid(c, c, c, indexing="ij")
    tr = np.stack([c12.ravel(), c23.ravel(), c31.ravel()], axis=1)
    g = np.arccos(np.clip(tr, -1, 1))
    g12, g23, g31 = g[:, 0], g[:, 1], g[:, 2]
    valid = ((g12 + g23 >= g31) & (g23 + g31 >= g12) & (g31 + g12 >= g23)
             & (g12 + g23 + g31 <= 2 * np.pi))
    return tr[valid]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--l2-npz", type=Path, default=None,
                    help="L2 NPZ providing the λ-grid (default: archived deployed L2).")
    ap.add_argument("--ell-cut", type=int, default=60)
    ap.add_argument("--ell-high-max", type=int, default=1000)
    ap.add_argument("--ell-min", type=float, default=1e-3)
    ap.add_argument("--n-ell", type=int, default=96)
    ap.add_argument("--n-phi", type=int, default=64)
    ap.add_argument("--Nmax", type=int, default=128)
    ap.add_argument("--n-cosine", type=int, default=15)
    ap.add_argument("--radial-measure", choices=("lambda", "chi"), default="lambda")
    ap.add_argument("--channels", nargs="+", default=["TTT", "TTP", "TPP", "PPP"])
    ap.add_argument("--smoke", action="store_true",
                    help="Tiny (ell-cut=20, ell-high-max=80, n-cosine=7) wiring check.")
    ap.add_argument("--out", type=Path, default=OUTPUT)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    from canoes.sachs import (  # noqa: E402
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    from canoes.sachs.kappa3 import (  # noqa: E402
        compute_kappa3_mod_zeta_equal_time,
        compute_kappa3_mod_sigma3_high,
    )

    l2_npz = args.l2_npz or _ARCHIVED_L2
    if not l2_npz.exists():
        print(f"error: L2 NPZ not found: {l2_npz}", file=sys.stderr)
        return 2
    if args.smoke:
        args.ell_cut, args.ell_high_max, args.n_cosine = 20, 80, 7
        out = _HERE / "_smoke_equal_time_limber.npz"
    else:
        out = args.out

    cosmo = FiducialCosmology()
    lam_phys, lam_h = _load_l2_lambda_grid_h(l2_npz, cosmo.h)
    pk_delta = load_pk_delta(n_s=cosmo.n_s)
    triples = _triangle_filtered_triples(int(args.n_cosine))
    channels = tuple(args.channels)
    print(f"[etl] cosmo Ω_m={cosmo.Omega_m}, h={cosmo.h}, n_s={cosmo.n_s}")
    print(f"[etl] λ-grid ({l2_npz.name}): n={lam_phys.size}, "
          f"[{lam_phys[0]:.1f}..{lam_phys[-1]:.1f}] Mpc; triples={triples.shape[0]}; "
          f"channels={channels}; measure={args.radial_measure}; "
          f"ell LOW≤{args.ell_cut}, HIGH({args.ell_cut},{args.ell_high_max}]")

    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        low = compute_kappa3_zeta_table(
            triples, lam_h, pk_matter_today=pk_delta, cosmo=cosmo,
            ell_max=int(args.ell_cut), Nmax=int(args.Nmax), spt_kind="tree_phi",
            units="physical", lambda_convention="project",
            radial_discretization="sample", radial_measure=args.radial_measure,
            channels=channels, parallelism="sequential", shell_block_size=4,
            triple_chunk_size=4096, n_workers=1, fuse_tree_phi_channels=True,
            verbose=bool(args.verbose))
        print(f"[etl] LOW done in {time.perf_counter()-t0:.1f}s")
        t1 = time.perf_counter()
        high = compute_kappa3_sigma3_high(
            triples, lam_h, pk_matter_today=pk_delta, cosmo=cosmo,
            ell_cut=int(args.ell_cut), ell_high_max=int(args.ell_high_max),
            ell_min=float(args.ell_min), n_ell=int(args.n_ell), n_phi=int(args.n_phi),
            channels=channels, units="physical", lambda_convention="project",
            radial_discretization="sample", radial_measure=args.radial_measure)
        print(f"[etl] HIGH done in {time.perf_counter()-t1:.1f}s")

    combined = kappa3_combine_low_high(low, high)

    # ── Conjugate-helicity modulus channels (Step 3 of spin-2 fix) ─────────
    # Bmod = <Phi00 Psi0 conj(Psi0)>  spin_weights (0, 2, -2)
    # Dmod = <Psi0  Psi0 conj(Psi0)>  spin_weights (2, 2, -2)
    # Both share the same SPT bispectrum as TPP/PPP respectively; only the
    # angular projection changes.  Units/measure identical to zeta_TPP/PPP.
    print("[etl] computing Bmod/Dmod LOW (equal-shell Wigner-3j) ...")
    t_mod = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        bmod_low, dmod_low = compute_kappa3_mod_zeta_equal_time(
            triples, lam_h,
            pk_matter_today=pk_delta, cosmo=cosmo,
            ell_max=int(args.ell_cut), Nmax=int(args.Nmax),
            spt_kind="tree_phi", units="physical",
            lambda_convention="project",
            radial_measure=args.radial_measure,
        )
    print(f"[etl] Bmod/Dmod LOW done in {time.perf_counter()-t_mod:.1f}s")

    print("[etl] computing Bmod/Dmod HIGH (Limber) ...")
    t_mod2 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        bmod_high, dmod_high = compute_kappa3_mod_sigma3_high(
            triples, lam_h,
            pk_matter_today=pk_delta, cosmo=cosmo,
            ell_cut=int(args.ell_cut), ell_high_max=int(args.ell_high_max),
            ell_min=float(args.ell_min), n_ell=int(args.n_ell),
            n_phi=int(args.n_phi), units="physical",
            lambda_convention="project",
            radial_measure=args.radial_measure,
        )
    print(f"[etl] Bmod/Dmod HIGH done in {time.perf_counter()-t_mod2:.1f}s")

    zeta_Bmod = bmod_low + bmod_high
    zeta_Dmod = dmod_low + dmod_high
    print(f"[etl] zeta_Bmod: shape={zeta_Bmod.shape}, "
          f"nonzero={int(np.count_nonzero(zeta_Bmod))}/{zeta_Bmod.size}, "
          f"finite={int(np.sum(np.isfinite(zeta_Bmod)))}")
    print(f"[etl] zeta_Dmod: shape={zeta_Dmod.shape}, "
          f"nonzero={int(np.count_nonzero(zeta_Dmod))}/{zeta_Dmod.size}, "
          f"finite={int(np.sum(np.isfinite(zeta_Dmod)))}")

    meta = dict(combined.cosmo_meta)
    meta.update(
        producer="sachs_sft/kappa3_vertex/equal_time_limber/build_equal_time_limber_table.py",
        artifact="equal_time_limber_lambda_NON_R_contracted",
        equal_time=True, limber_high=True,
        already_R_convolved=False, already_R_contracted=False,
        radial_measure=args.radial_measure, radial_density_measure=args.radial_measure,
        radial_kernel_convention=(
            "lambda_three_leg_native_density" if args.radial_measure == "lambda"
            else "chi_density_three_leg"),
        sachs_prefactor_per_leg="(1+z)^4 bare driving field (SACHS_GEOMETRIC_POWER=4)",
        poisson_policy="P_matter input; tree_phi applies A(a) per leg internally",
        query_coordinate="lambda_project_mpc",
        l2_lambda_grid_npz=l2_npz.name,
        ell_cut=int(args.ell_cut), ell_high_max=int(args.ell_high_max),
        n_cosine=int(args.n_cosine),
        fk_kk_baseline=(
            "+1.219e-4 (γ=0.5', LIMBER grid equal_time (1+z)^4/leg λ-density, ng32; "
            "h^4-corrected 2026-06-09). NOTE: the pre-fix value was +3.086e-5; the "
            "h^6->h^4 units fix in canoes kappa3.py (branch fix/spin2-zetaD-2_2_-2-norm, "
            "the 5 equal-shell sites) and the consequent table/FK rebuild shifted FK kk "
            "by x2.2204 to +1.219e-4. Order-0 CCL kk unchanged (0.998)."
        ),
        has_modulus_channels=True,
    )
    combined = replace(combined, cosmo_meta=meta)
    for ch in ("TTT", "TTP", "TPP", "PPP"):
        arr = np.asarray(getattr(combined, f"zeta_{ch}"))
        print(f"[etl] {'*' if ch in channels else ' '} zeta_{ch}: shape={arr.shape}, "
              f"nonzero={int(np.count_nonzero(arr))}/{arr.size}")
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.save_npz(out)

    # ── Post-save merge: inject modulus arrays into the frozen-schema NPZ ───
    # Kappa3Output.save_npz is schema-frozen (4 canonical channels + cosmo_meta).
    # We extend the on-disk file by reading it back and writing a merged copy.
    #
    # allow_pickle=True is required because save_npz stores cosmo_meta as a
    # pickled 0-d object array (a canoes internal convention).  This is safe
    # here: the file was written by this script seconds ago from trusted local
    # computations; it is NOT user-supplied or downloaded data.
    _existing = dict(np.load(out, allow_pickle=True))
    np.savez(
        out,
        **_existing,
        zeta_Bmod=zeta_Bmod,
        zeta_Dmod=zeta_Dmod,
        has_modulus_channels_flag=np.bool_(True),
    )

    print(f"[etl] total {time.perf_counter()-t0:.1f}s -> {out} "
          f"({out.stat().st_size/1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
