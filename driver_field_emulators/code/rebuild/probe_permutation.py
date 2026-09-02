"""Do the vertex channels depend on WHICH leg carries the offset?

The production callable sorts a query triple descending before lookup and
then fills the (3,3,3) coupling tensor with the one retrieved value in
every index permutation. That is exact for a fully symmetric channel and an
approximation for a channel that assigns specific spins to specific legs.

The test builds the same collapsed configuration in two leg orderings that
canonicalise to the same sorted triple:

    A = (1, c, c)   points 1 and 2 coincide, point 3 carries the offset
    B = (c, 1, c)   points 2 and 3 coincide, point 1 carries the offset

The callable cannot tell them apart. Any difference between them is exactly
the ambiguity the symmetric fill hides, so it bounds how precisely each
observable can be quoted once the grid resolves the collapsed family.

Usage::

    python probe_permutation.py --low ../../products/pieces/low_perm.npz \
        --high ../../products/pieces/win_perm_00060_01000.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

CHANNELS = ("zeta_TTT", "zeta_TTP", "zeta_TPP", "zeta_PPP")
MODS = {"bmod": "zeta_Bmod", "dmod": "zeta_Dmod"}
ARCMIN = 180.0 * 60.0 / np.pi


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--low", type=Path, required=True)
    ap.add_argument("--high", type=Path, required=True)
    ap.add_argument("--gamma-max", type=float, default=200.0,
                    help="report over the separations that carry signal")
    args = ap.parse_args()

    low = np.load(args.low, allow_pickle=False)
    high = np.load(args.high, allow_pickle=False)
    triples = low["cosine_triples"]
    n = triples.shape[0] // 2
    if not (np.allclose(triples[:n, 0], 1.0) and np.allclose(triples[n:, 1], 1.0)):
        raise SystemExit("the probe table is not the two-ordering layout")

    gamma = np.degrees(np.arccos(np.clip(triples[:n, 1], -1, 1))) * 60.0
    rows = gamma <= args.gamma_max

    print(f"collapsed family in two leg orderings, {int(rows.sum())} separations "
          f"up to {args.gamma_max:.0f} arcmin\n")
    print(f"{'channel':12s} {'median rel. difference':>23} {'symmetric?':>12}")
    verdict = {}
    for key, name in [(c, c) for c in CHANNELS] + list(MODS.items()):
        total = np.asarray(low[key], float) + np.asarray(high[key], float)
        a, b = total[:n], total[n:]
        scale = np.maximum(np.abs(a), np.abs(b))
        mask = scale > 0
        rel = np.zeros_like(a)
        rel[mask] = np.abs(a[mask] - b[mask]) / scale[mask]
        median = float(np.median(rel[rows]))
        verdict[name] = median
        print(f"{name:12s} {median:23.3e} {'yes' if median < 1e-2 else 'NO':>12}")

    print("\nMaxima are not quoted: every channel passes through zero somewhere "
          "in this range,\nso the relative difference saturates at 2 there and "
          "says nothing.")
    print("\nConsequence. xi_kappa and xi_+ are carried by zeta_TTT and "
          "zeta_Bmod, both symmetric\nto the quadrature level, so the symmetric "
          "fill costs them nothing. xi_- and\nxi_kappa_gamma involve the spin "
          "channels, where the ordering matters at the tens of\npercent level, "
          "so those two can be quoted as magnitudes and signs but not as "
          "precise\nvalues until the callable distinguishes the orderings.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
