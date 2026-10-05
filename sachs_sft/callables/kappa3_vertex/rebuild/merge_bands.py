"""Concatenate two band pieces built on disjoint row sets.

The permutation-closed triple set is the r4 set followed by 314 extra rows,
verified by exact array comparison. Rows are independent bit-exactly, so a
band on the union is the concatenation of a band on each part, with no
approximation. This reuses the octave bands already built on r4 and needs
only the extra rows built, which is 15% of the work.

Usage::

    python merge_bands.py --parts band_tree_r4_X.npz band_tree_extra_X.npz \
        --triples products/triples_permclosed_r4.npz --out merged.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from _common import (CHANNELS, build_leg_order, build_pk_identity, dump_meta,
                     grid_fingerprint, load_meta)

ROW_ARRAYS = tuple(f"zeta_{c}" for c in CHANNELS) + ("bmod", "dmod")
SHARED = ("lambda_shells", "chi_shells", "z_shells")


def merged_fingerprint(parts, triples: np.ndarray) -> dict:
    """The parts' shell fingerprint with the row-set fields recomputed."""
    base = load_meta(parts[0]["fingerprint"])
    rows = grid_fingerprint(triples, np.array([base["lam_first"], base["lam_last"]]))
    base["n_triples"] = rows["n_triples"]
    base["triples_hash"] = rows["triples_hash"]
    return base


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parts", nargs="+", type=Path, required=True)
    ap.add_argument("--triples", type=Path, required=True,
                    help="the union triple set; row order must match the "
                         "concatenation of the parts, and is checked")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    parts = [np.load(p, allow_pickle=False) for p in args.parts]
    builds = [load_meta(p["build"]) for p in parts]
    for key in ("lo", "hi", "n_ell", "n_phi", "b_model", "radial_measure"):
        values = {b[key] for b in builds}
        if len(values) != 1:
            raise SystemExit(f"parts disagree on {key}: {sorted(values)}")
    pk_identities = {build_pk_identity(b) for b in builds}
    if len(pk_identities) != 1:
        raise SystemExit(f"parts were built on different P(k) table identities: {sorted(pk_identities, key=repr)}")
    pk_table, pk_sha = next(iter(pk_identities))
    leg_orders = {build_leg_order(b) for b in builds}
    if len(leg_orders) != 1:
        raise SystemExit(f"parts were built with different leg orders: {sorted(leg_orders)}")
    for name in SHARED:
        if not all(np.array_equal(parts[0][name], p[name]) for p in parts[1:]):
            raise SystemExit(f"parts were built on different {name}")

    triples = np.asarray(np.load(args.triples)["cosine_triples"], float)
    stacked = np.concatenate([p["cosine_triples"] for p in parts], axis=0)
    if not np.array_equal(stacked, triples):
        raise SystemExit(
            "the concatenated rows are not the union triple set in order; "
            "merging would silently mislabel every row")

    merged = {name: np.concatenate([p[name] for p in parts], axis=0)
              for name in ROW_ARRAYS}
    build = dict(builds[0])
    build["pk_table"] = pk_table
    build.pop("pk_table_sha256", None)
    if pk_sha is not None:
        build["pk_table_sha256"] = pk_sha
    build["merged_from"] = [p.name for p in args.parts]
    build["n_rows"] = int(triples.shape[0])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, cosine_triples=triples,
             **{n: parts[0][n] for n in SHARED},
             ell_max=parts[0]["ell_max"],
             cosmo_meta=parts[0]["cosmo_meta"],
             # The fingerprint must be rebuilt from the same lambda array the
             # parts used, which is the PHYSICAL grid, not the h-native one
             # stored in lambda_shells. Take the parts' own fingerprint and
             # update only the row-set fields, so the shell fields stay
             # byte-consistent with the LOW piece assemble.py will check
             # against.
             fingerprint=dump_meta(merged_fingerprint(parts, triples)),
             build=dump_meta(build), **merged)
    print(f"[merge] ({build['lo']}, {build['hi']}]  "
          f"{' + '.join(str(p['cosine_triples'].shape[0]) for p in parts)}"
          f" = {triples.shape[0]} rows -> {args.out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
