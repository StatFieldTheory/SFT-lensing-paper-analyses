"""How much does linear interpolation between the 16 radial shells cost?

The callable interpolates linearly in lambda between the tabulated shells.
Building the vertex ON the midpoints between those shells gives the true
value there, so the deployed grid's interpolant can be compared against
truth directly, with no other change.

Usage::

    python probe_radial_interp.py --low products/pieces/low_collapsed_mid.npz \
        --high products/pieces/win_collapsed_mid_00060_01000.npz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

ARCMIN = 180.0 * 60.0 / np.pi
GAMMAS = (0.5, 5.30, 21.88, 114.27)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--low", type=Path, required=True)
    ap.add_argument("--high", type=Path, required=True)
    ap.add_argument("--channel", default="zeta_TTT")
    args = ap.parse_args()

    low = np.load(args.low, allow_pickle=False)
    high = np.load(args.high, allow_pickle=False)
    zeta = np.asarray(low[args.channel], float) + np.asarray(high[args.channel], float)
    lam = np.asarray(low["lambda_shells"], float)
    triples = low["cosine_triples"]

    # The fine grid is the 16 deployed shells with midpoints inserted, so the
    # deployed shells are the even indices and the midpoints the odd ones.
    deployed = zeta[:, 0::2]
    lam_deployed = lam[0::2]
    truth = zeta[:, 1::2]
    lam_mid = lam[1::2]
    interpolated = 0.5 * (deployed[:, :-1] + deployed[:, 1:])

    print(f"channel {args.channel}: linear interpolation between the deployed "
          f"shells against the true midpoint value")
    print(f"  {len(lam_deployed)} deployed shells, {len(lam_mid)} midpoints\n")
    header = "  midpoint lambda [Mpc]  " + "".join(
        f"{g:>12.2f}'" for g in GAMMAS)
    print(header)
    rows = []
    for g in GAMMAS:
        c = np.cos(g / ARCMIN)
        target = np.array([1.0, c, c])
        i = int(np.argmin(np.linalg.norm(triples - target, axis=1)))
        rows.append(i)
    for j, lm in enumerate(lam_mid):
        cells = []
        for i in rows:
            t = truth[i, j]
            cells.append(interpolated[i, j] / t if t != 0 else np.nan)
        print(f"  {lm:21.1f}  " + "".join(f"{v:13.3f}" for v in cells))

    ratio = np.where(truth[rows] != 0, interpolated[rows] / truth[rows], np.nan)
    finite = np.isfinite(ratio)
    print(f"\n  median interpolant/truth = {np.median(ratio[finite]):.3f}")
    print(f"  range                    = [{np.nanmin(ratio):.3f}, "
          f"{np.nanmax(ratio):.3f}]")
    print("\nThe outermost shells carry nearly all the lensing weight, so the "
          "bias that\nmatters is the one at large lambda, not the largest one "
          "in the table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
