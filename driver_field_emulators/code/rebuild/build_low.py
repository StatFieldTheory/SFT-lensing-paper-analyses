"""Build the LOW (exact Wigner-3j) branch of the kappa3 vertex, once.

LOW is independent of the HIGH cutoff and of the bispectrum model: it is
always tree-level, because its FFTlog machinery requires a separable
bispectrum. One LOW build therefore serves every variant built on the same
triples, shells, ell_cut, Nmax and radial measure.

Usage::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1 \
      $CAN/.venv/bin/python -u build_low.py \
        --triples-npz ../../products/triples_full_r4.npz \
        --out ../../products/pieces/low_r4.npz
"""

from __future__ import annotations

import argparse
import time
import warnings
from pathlib import Path

import numpy as np

from _common import (ARCHIVED_L2, CHANNELS, bootstrap_paths, dump_meta,
                     grid_fingerprint, load_cosmology, load_lambda_grid,
                     load_triples)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--triples-npz", type=Path, required=True)
    ap.add_argument("--l2-npz", type=Path, default=ARCHIVED_L2)
    ap.add_argument("--ell-cut", type=int, default=60)
    ap.add_argument("--Nmax", type=int, default=128)
    ap.add_argument("--radial-measure", choices=("lambda", "chi"), default="lambda")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bootstrap_paths()
    from canoes.sachs import compute_kappa3_zeta_table
    from canoes.sachs.kappa3 import compute_kappa3_mod_zeta_equal_time

    cosmo, pk = load_cosmology()
    lam_phys, lam_h = load_lambda_grid(args.l2_npz, cosmo.h)
    triples = load_triples(args.triples_npz)
    print(f"[low] triples={triples.shape[0]}  shells={lam_phys.size}  "
          f"ell_cut={args.ell_cut}  Nmax={args.Nmax}  "
          f"measure={args.radial_measure}")

    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        low = compute_kappa3_zeta_table(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_max=int(args.ell_cut), Nmax=int(args.Nmax),
            spt_kind="tree_phi", units="physical",
            lambda_convention="project", radial_discretization="sample",
            radial_measure=args.radial_measure, channels=CHANNELS,
            parallelism="sequential", shell_block_size=4,
            triple_chunk_size=4096, n_workers=1,
            fuse_tree_phi_channels=True, verbose=bool(args.verbose))
        print(f"[low] four channels in {time.perf_counter() - t0:.1f}s")

        t1 = time.perf_counter()
        bmod, dmod = compute_kappa3_mod_zeta_equal_time(
            triples, lam_h, pk_matter_today=pk, cosmo=cosmo,
            ell_max=int(args.ell_cut), Nmax=int(args.Nmax),
            spt_kind="tree_phi", units="physical",
            lambda_convention="project",
            radial_measure=args.radial_measure)
        print(f"[low] Bmod/Dmod in {time.perf_counter() - t1:.1f}s")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        cosine_triples=low.cosine_triples,
        lambda_shells=low.lambda_shells,
        chi_shells=low.chi_shells,
        z_shells=low.z_shells,
        zeta_TTT=low.zeta_TTT, zeta_TTP=low.zeta_TTP,
        zeta_TPP=low.zeta_TPP, zeta_PPP=low.zeta_PPP,
        bmod=bmod, dmod=dmod,
        ell_max=np.int64(low.ell_max),
        cosmo_meta=dump_meta(dict(low.cosmo_meta)),
        fingerprint=dump_meta(grid_fingerprint(triples, lam_phys)),
        build=dump_meta(dict(
            kind="low", ell_cut=int(args.ell_cut), Nmax=int(args.Nmax),
            radial_measure=args.radial_measure,
            l2_lambda_grid_npz=args.l2_npz.name,
            triples_npz=args.triples_npz.name)),
    )
    print(f"[low] total {time.perf_counter() - t0:.1f}s -> {args.out} "
          f"({args.out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
