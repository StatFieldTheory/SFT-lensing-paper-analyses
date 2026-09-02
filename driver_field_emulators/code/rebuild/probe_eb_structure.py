"""What actually carries E minus B in the FK term.

The paper's claim is that the FK power splits equally between the two
shear polarizations, Delta C_EE = Delta C_BB, because the difference
<Phi_00 Psi_+ Psi_+> - <Phi_00 Psi_x Psi_x> vanishes for a statistically
isotropic driving field.

That difference is a specific vertex channel. In the production callable's
real-basis reconstruction the coupling tensor is filled as

    K011 = (zeta_TPP + zeta_Bmod) / 2
    K022 = (zeta_Bmod - zeta_TPP) / 2
    K111 = (zeta_PPP + 3 zeta_Dmod) / 4
    K122 = (zeta_Dmod - zeta_PPP) / 4

so the combinations that matter are read off directly rather than derived:
the SUM K011 + K022 and K111 + K122 feed xi_+ (that is C_EE + C_BB), and
the DIFFERENCE K011 - K022 and K111 - K122 feed xi_- (that is
C_EE - C_BB). This script verifies those combinations numerically, then
measures the channels that survive in the difference, as a function of
separation.

If the difference channel vanishes only at zero separation, the equal split
is a zero-separation identity and not a property of the plotted range.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

ARCMIN = 180.0 * 60.0 / np.pi


def verify_combinations() -> None:
    """Check the sum and difference combinations on random channel values."""
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
    print("Combination check (exact arithmetic on the callable's own fill):")
    for label, (lhs, rhs) in checks.items():
        print(f"  {label:42s} residual {abs(lhs - rhs):.3e}")
    print()


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
    collapsed = np.isclose(triples[:, 0], 1.0) & np.isclose(triples[:, 1], triples[:, 2])
    idx = np.flatnonzero(collapsed)
    gamma = np.degrees(np.arccos(np.clip(triples[idx, 1], -1, 1))) * 60.0
    order = np.argsort(gamma)
    idx, gamma = idx[order], gamma[order]

    ch = {name: np.asarray(d[f"zeta_{name}"], float)[idx, args.shell]
          for name in ("TTT", "TTP", "TPP", "PPP", "Bmod", "Dmod")}
    sum_spin = ch["Bmod"]
    diff_spin = ch["TPP"]
    sum_pure = ch["Dmod"]
    diff_pure = 0.5 * (ch["PPP"] + ch["Dmod"])

    print(f"collapsed family on shell lambda = {lam[args.shell]:.1f} Mpc\n")
    print(f"{'gamma':>9} {'E+B: Bmod':>13} {'E-B: TPP':>13} {'|TPP/Bmod|':>12}"
          f" {'E+B: Dmod':>13} {'E-B: (PPP+Dmod)/2':>19}")
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
    print("\nThe difference channel is the one the equal-split claim needs to "
          "vanish.\nIf it is small only as gamma -> 0, the claim is a "
          "zero-separation identity.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
