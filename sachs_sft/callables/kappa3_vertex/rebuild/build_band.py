"""Build one HIGH (flat-sky Born-Limber) window of the kappa3 vertex.

The HIGH branch integrates over multipoles in a window (lo, hi] that is
disjoint in max-leg from every other window, so a set of windows covering
(ell_cut, ell_max] sums to the same integral as one run over the whole
range. Building windows separately buys two things: they run in parallel,
and partial sums give every intermediate cutoff without rebuilding.

Node density is what makes the windowing worthwhile rather than merely
convenient. A single Gauss-Legendre rule with n_ell nodes mapped linearly
onto [ell_min, ell_max] coarsens as ell_max grows, so a naive cutoff sweep
at fixed n_ell mixes genuine truncation with a quadrature artifact. One
rule per octave keeps the node density per decade fixed.

Usage::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1 \
      $CAN/.venv/bin/python -u build_band.py \
        --triples-npz products/triples_full_r4.npz \
        --lo 60 --hi 120 --b-model tree \
        --out products/pieces/band_tree_r4_00060_00120.npz
"""

from __future__ import annotations

import argparse
import time
import warnings
from pathlib import Path

import numpy as np

from _common import (ARCHIVED_L2, CHANNELS, LEG_ORDERS, PK_TABLE_DEFAULT,
                     PK_TABLES, bootstrap_paths, canoes_record, dump_meta,
                     grid_fingerprint, load_cosmology, load_lambda_grid,
                     load_triples, pk_table_record)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--triples-npz", type=Path, required=True)
    ap.add_argument("--l2-npz", type=Path, default=ARCHIVED_L2)
    ap.add_argument("--lo", type=int, required=True, help="window lower edge (exclusive)")
    ap.add_argument("--hi", type=int, required=True, help="window upper edge (inclusive)")
    ap.add_argument("--n-ell", type=int, default=96)
    ap.add_argument("--n-phi", type=int, default=64)
    ap.add_argument("--ell-min", type=float, default=1e-3)
    ap.add_argument("--radial-measure", choices=("lambda", "chi"), default="lambda")
    ap.add_argument("--b-model", default="internal_tree",
                    choices=("internal_tree", "tree", "bihalofit", "bihalofit_baryon"))
    ap.add_argument("--pk-table", default=PK_TABLE_DEFAULT, choices=PK_TABLES,
                    help="linear P(k) table (canoes examples/data) for the "
                         "bispectrum model, pk_matter_today and the BiHalofit "
                         "sigma8; recorded with its sha256 in the build metadata")
    ap.add_argument("--leg-order", default="auto", choices=LEG_ORDERS,
                    help="canoes leg_order: 'auto' evaluates spin-sum-zero channels "
                         "(TTT, Bmod) of rows whose closest pair is points 2 and 3 "
                         "in the cyclic order; 'given' keeps the given order")
    ap.add_argument("--chunk-rows", type=int, default=256,
                    help="rows per call into canoes. The HIGH branch holds a "
                         "(rows, n_ell, n_ell, n_phi) working array, which at "
                         "1828 rows and the default quadrature is about 8 GB. "
                         "Rows are independent (verified bit-exactly), so "
                         "chunking is exact and caps the footprint.")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    bootstrap_paths()
    from canoes.sachs import compute_kappa3_sigma3_high
    from canoes.sachs.kappa3 import compute_kappa3_mod_sigma3_high

    cosmo, pk = load_cosmology(args.pk_table)
    lam_phys, lam_h = load_lambda_grid(args.l2_npz, cosmo.h)
    triples = load_triples(args.triples_npz)

    b_delta_fn = None
    if args.b_model != "internal_tree":
        from b_model import make_b_delta_fn
        from b_model.cosmology import fiducial_for_table
        b_delta_fn = make_b_delta_fn(args.b_model, cosmo=fiducial_for_table(args.pk_table),
                                     pk_table=args.pk_table)

    canoes_state = canoes_record()
    print(f"[band] ({args.lo}, {args.hi}]  model={args.b_model}  pk={args.pk_table}  "
          f"leg_order={args.leg_order}  n_ell={args.n_ell}  n_phi={args.n_phi}  "
          f"triples={triples.shape[0]}  canoes={canoes_state['canoes_commit'][:7]}"
          f"{'' if canoes_state['canoes_src_clean'] else ' (src modified)'}")

    shared = dict(pk_matter_today=pk, cosmo=cosmo, ell_min=float(args.ell_min),
                  n_ell=int(args.n_ell), n_phi=int(args.n_phi),
                  units="physical", lambda_convention="project",
                  radial_discretization="sample",
                  radial_measure=args.radial_measure, b_delta_fn=b_delta_fn,
                  leg_order=args.leg_order)

    mod_shared = {k: v for k, v in shared.items() if k != "radial_discretization"}
    chunk = max(1, int(args.chunk_rows))
    blocks = [slice(i, min(i + chunk, triples.shape[0]))
              for i in range(0, triples.shape[0], chunk)]

    t0 = time.perf_counter()
    parts, bmods, dmods, template = [], [], [], None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        for n, block in enumerate(blocks, 1):
            rows = triples[block]
            part = compute_kappa3_sigma3_high(
                rows, lam_h, ell_cut=int(args.lo), ell_high_max=int(args.hi),
                channels=CHANNELS, **shared)
            b_part, d_part = compute_kappa3_mod_sigma3_high(
                rows, lam_h, ell_cut=int(args.lo), ell_high_max=int(args.hi),
                **mod_shared)
            template = template or part
            parts.append({c: np.asarray(getattr(part, f"zeta_{c}"), float)
                          for c in CHANNELS})
            bmods.append(np.asarray(b_part, float))
            dmods.append(np.asarray(d_part, float))
            del part, b_part, d_part
            print(f"[band]   rows {block.start}-{block.stop} "
                  f"({n}/{len(blocks)}) at {time.perf_counter() - t0:.0f}s",
                  flush=True)

    channels_out = {f"zeta_{c}": np.concatenate([p[c] for p in parts], axis=0)
                    for c in CHANNELS}
    bmod = np.concatenate(bmods, axis=0)
    dmod = np.concatenate(dmods, axis=0)
    high = template
    print(f"[band] all rows in {time.perf_counter() - t0:.1f}s")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        cosine_triples=triples,
        lambda_shells=high.lambda_shells,
        chi_shells=high.chi_shells,
        z_shells=high.z_shells,
        bmod=bmod, dmod=dmod,
        **channels_out,
        ell_max=np.int64(high.ell_max),
        cosmo_meta=dump_meta(dict(high.cosmo_meta)),
        fingerprint=dump_meta(grid_fingerprint(triples, lam_phys)),
        build=dump_meta(dict(
            kind="band", lo=int(args.lo), hi=int(args.hi),
            n_ell=int(args.n_ell), n_phi=int(args.n_phi),
            chunk_rows=chunk,
            b_model=args.b_model, radial_measure=args.radial_measure,
            **pk_table_record(args.pk_table),
            leg_order=args.leg_order, **canoes_state,
            l2_lambda_grid_npz=args.l2_npz.name,
            triples_npz=args.triples_npz.name)),
    )
    print(f"[band] total {time.perf_counter() - t0:.1f}s -> {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
