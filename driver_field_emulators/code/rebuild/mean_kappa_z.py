"""Deterministic second-order mean convergence <kappa>(z_s), no Monte-Carlo.

The FF moment tends, at large separation, to the constant disconnected part
<kappa>^2, the square of the second-order mean convergence. The multi-z
figure needs that plateau per source redshift, and the talk-era cache never
computed it (it scaled the z=5 FF by a per-gamma FK ratio, which is noise at
the separations where the plateau shows). This script propagates the exact
discrete recursions of ``simulate_ff_crn`` in expectation:

    V_k = resp_k^2 V_{k-1} + Sigma2(lam_k) dlam          (Gaussian covariance)
    m_k = resp_k m_{k-1} + <F(s,s)>(V_{k-1}) dlam        (mean of the F arm)
    <kappa> = sum_k w_k m_k[0]

with <F(s,s)> = (-(V00+V11+V22), -2 V01, -2 V02) from the explicit F vertex.
The white-noise (equal-time) collapse bias is common to all z, so plateau
RATIOS anchored at the production z=5 value are accurate to the residual
z-dependence of that bias.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[3]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"


def mean_kappa(lam_source: float, n_lambda: int = 2000, lam_min: float = 406.0):
    sys.path.insert(0, str(_MC))
    import background as _bg
    import driver_stats as _ds
    bg = _bg.Background(lam_min=120.0, lam_max=lam_source + 5.0)
    lam = np.linspace(lam_min, lam_source, n_lambda)
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), float)
    resp = np.empty(n_lambda); resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    builder = _ds.Sigma2Builder(background=bg, apply_c0=False)

    V = np.zeros((3, 3)); m = np.zeros(3); kap = np.zeros(3)
    for k in range(n_lambda):
        fv = np.array([-(V[0, 0] + V[1, 1] + V[2, 2]),
                       -2.0 * V[0, 1], -2.0 * V[0, 2]])
        m = resp[k] * m + fv * dlam
        S2 = builder.matrix(1.0, float(lam[k]))
        V = resp[k] ** 2 * V + S2 * dlam
        w = 0.5 * dlam if (k == 0 or k == n_lambda - 1) else dlam
        kap += w * m
    return float(kap[0])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lam-sources", type=float, nargs="+",
                    default=[1822.721, 2094.894, 2216.793, 2266.582,
                             2297.289, 2313.0288751857356])
    ap.add_argument("--labels", type=str, nargs="+",
                    default=["z=1.0", "z=1.7", "z=2.5", "z=3.2", "z=4.0", "z=5.0"])
    ap.add_argument("--n-lambda", type=int, default=2000)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    rows = []
    print(f"{'label':>8} {'lam_s':>10} {'<kappa>':>12} {'<kappa>^2':>12}")
    for lab, ls in zip(a.labels, a.lam_sources):
        mk = mean_kappa(ls, a.n_lambda)
        rows.append((ls, mk, mk * mk))
        print(f"{lab:>8} {ls:10.2f} {mk:+12.5e} {mk*mk:12.5e}", flush=True)
    P5 = rows[-1][2]
    print("\nplateau ratios P(z)/P(z=5):",
          "  ".join(f"{r[2]/P5:.4f}" for r in rows))
    print(f"production z=5 plateau: 5.785e-07; computed z=5: {P5:.4e} "
          f"(ratio {P5/5.785e-7:.3f})")
    if a.out:
        arr = np.array(rows)
        np.savez(a.out, lam_source=arr[:, 0], mean_kappa=arr[:, 1],
                 plateau=arr[:, 2], n_lambda=a.n_lambda)
        print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
