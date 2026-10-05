"""Historical single-placement channel diagnostic, not a current FK validation.

The algebra below retains the historical K111=(PPP+3Dmod)/4 fill at one
placement. The current v2 callable instead reconstructs K111 using Dmod
at three cyclic leg placements. This script does not implement or validate
that reconstruction or the folded xi_plus/xi_minus and E/B spectra.

Printed channel combinations and raw ratios are descriptive diagnostics
only. They do not establish equal E/B power in an accepted FK fold.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

ARCMIN = 180.0 * 60.0 / np.pi


def verify_combinations() -> None:
    """Check only the historical single-placement algebra on random values."""
    rng = np.random.default_rng(0)
    tpp, bmod, ppp, dmod = rng.normal(size=4)
    k011 = 0.5 * (tpp + bmod)
    k022 = 0.5 * (bmod - tpp)
    k111 = 0.25 * (ppp + 3.0 * dmod)
    k122 = 0.25 * (dmod - ppp)
    checks = {
        "K011 + K022 = zeta_Bmod": (k011 + k022, bmod),
        "K011 - K022 = zeta_TPP": (k011 - k022, tpp),
        "K111 + K122 = zeta_Dmod": (k111 + k122, dmod),
        "K111 - K122 = (zeta_PPP + zeta_Dmod)/2": (k111 - k122,
                                                   0.5 * (ppp + dmod)),
    }
    print("Historical single-placement algebra check (not current callable validation):")
    for label, (lhs, rhs) in checks.items():
        print(f"  {label:42s} residual {abs(lhs - rhs):.3e}")
    print()


def collapsed_slot0_rows(triples):
    """Select slot0 with exact 12-decimal keys, refusing duplicate rows."""
    t = np.asarray(triples, dtype=float)
    if t.ndim != 2 or t.shape[1] != 3 or not np.isfinite(t).all():
        raise ValueError("Expected finite cosine triples (n,3)")
    keys = np.round(t, 12)
    if len(np.unique(keys, axis=0)) != len(keys):
        raise ValueError("Duplicate 12-decimal geometry keys")
    mask = (keys[:, 0] == 1.0) & (keys[:, 1] == keys[:, 2])
    idx = np.flatnonzero(mask)
    if not idx.size:
        raise ValueError("No collapsed slot0 rows")
    if len(np.unique(keys[idx, 1])) != len(idx):
        raise ValueError("Duplicate slot0 separations")
    return idx


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--shell", type=int, default=-1,
                    help="radial shell index; default the source shell")
    args = ap.parse_args()

    verify_combinations()

    d = np.load(args.table, allow_pickle=False)
    triples = np.asarray(d["cosine_triples"], float)
    lam = np.asarray(d["lambda_shells_Mpc"], float)
    idx = collapsed_slot0_rows(triples)
    gamma = np.degrees(np.arccos(np.clip(triples[idx, 1], -1, 1))) * 60.0
    order = np.argsort(gamma)
    idx, gamma = idx[order], gamma[order]

    ch = {name: np.asarray(d[f"zeta_{name}"], float)[idx, args.shell]
          for name in ("TTT", "TTP", "TPP", "PPP", "Bmod", "Dmod")}
    sum_spin = ch["Bmod"]
    diff_spin = ch["TPP"]
    sum_pure = ch["Dmod"]
    diff_pure = 0.5 * (ch["PPP"] + ch["Dmod"])

    print(f"Raw-channel slot0 diagnostic on shell lambda = {lam[args.shell]:.1f} Mpc\n")
    print(f"{'gamma':>9} {'raw Bmod':>13} {'raw TPP':>13} {'|TPP/Bmod|':>12}"
          f" {'raw Dmod':>13} {'historical (PPP+Dmod)/2':>19}")
    for i, g in enumerate(gamma):
        if g > 400:
            continue
        ratio = abs(diff_spin[i] / sum_spin[i]) if sum_spin[i] else np.nan
        print(f"{g:9.2f} {sum_spin[i]:+13.4e} {diff_spin[i]:+13.4e} "
              f"{ratio:12.4f} {sum_pure[i]:+13.4e} {diff_pure[i]:+19.4e}")

    small = gamma < 1.0
    print(f"\n  |TPP/Bmod| at the smallest separation ({gamma[0]:.2f}'): "
          f"{abs(diff_spin[0]/sum_spin[0]):.3e}")
    mid = (gamma > 10) & (gamma < 60)
    if mid.any():
        print(f"  |TPP/Bmod| median over 10' to 60': "
              f"{np.median(np.abs(diff_spin[mid]/sum_spin[mid])):.3f}")
    print("\nRaw-channel ratios are descriptive only; the historical single-placement "
          "algebra does not validate current FK xi_plus/xi_minus or equal E/B power.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
