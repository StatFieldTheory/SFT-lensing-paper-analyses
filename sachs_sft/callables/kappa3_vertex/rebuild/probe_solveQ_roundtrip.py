"""Does the Monte-Carlo's deformation tensor Q actually reproduce its input cumulant?

The FK Monte-Carlo does not read the three-point vertex directly.  It solves
for a deformation tensor Q from

    Zeta6_abc = Q_amn M2_mb M2_nc + Q_bmn M2_ma M2_nc + Q_cmn M2_ma M2_nb,

which is Q ~ Zeta6 . M2^-1 . M2^-1, and ``solve_Q`` regularises that inverse by
ZEROING every eigendirection of M2 whose eigenvalue falls below
``rcond * d_max`` (rcond = 1e-8).  The stated justification is that for
near-aligned rays the small-eigenvalue directions carry ~0 cumulant too, so
dropping them costs nothing.  That is a threshold dispatch, and the place to
test a threshold is AT the threshold: at intermediate separation the
antisymmetric ray combinations have small-but-not-negligible variance, and the
eigenvalues sweep across the floor.

This probe is deterministic -- no realisations.  For each separation it
rebuilds M2 and Zeta6 exactly as ``_field_precompute`` does, solves for Q, maps
it forward again with ``cum3_from_Q``, and reports how much of the input
cumulant survives the round trip, together with the eigenvalue spectrum
relative to the floor.

If the round-trip residual is small at 1' and large near 17', the estimator is
being fed a different vertex from the one the analytic fold uses, and the
Monte-Carlo/analytic disagreement is explained without invoking new physics.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[5]
_MC = _REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
_FIX = Path(__file__).resolve().parents[1] / "equal_time_limber_cut15360_permaware"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.6, 6.7, 17.3, 44.4, 114.3, 293.9, 957.2])
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--n-lambda", type=int, default=1000,
                    help="the estimator's grid; the probe subsamples it")
    ap.add_argument("--n-probe", type=int, default=60,
                    help="how many lambda nodes to actually test")
    ap.add_argument("--lam-min", type=float, default=406.0)
    ap.add_argument("--rcond", type=float, default=1e-8)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    sys.path.insert(0, str(_MC))
    sys.path.insert(0, str(_FIX))
    import driver_stats as ds
    import perm_aware_kappa3_callable as pa
    pa.TABLE_PATH = a.table.resolve()
    pa._CACHE = None
    ds._k3 = pa
    import background as _bg
    import sachs_mc_core as mc

    bg = _bg.Background(lam_min=120.0, lam_max=_bg.LAM_SOURCE_BASELINE + 5.0)
    lam_full = np.linspace(a.lam_min, _bg.LAM_SOURCE_BASELINE, a.n_lambda)
    idx = np.unique(np.linspace(0, a.n_lambda - 1, a.n_probe).astype(int))
    builder = ds.Sigma2Builder(background=bg, apply_c0=False)
    sig = a.sigma

    print(f"[probe] sigma_lambda={sig}, {len(idx)} of {a.n_lambda} lambda nodes, "
          f"rcond={a.rcond:g}, table={a.table.name}\n")
    # The Frobenius norm of the (6,6,6) cumulant is dominated by the COINCIDENT
    # blocks (all three legs on one ray): those are large and do not decay with
    # separation.  The FK estimator contracts the CROSS blocks, which are much
    # smaller.  A round-trip error that is 3% of the total norm can therefore be
    # several hundred percent of the block the observable actually uses, so the
    # residual is reported per block: 3+0 = all three legs on ray 1, 2+1 = two on
    # ray 1, and so on.
    ray = np.array([0, 0, 0, 1, 1, 1])
    blocks = {}
    for a3 in range(6):
        for b3 in range(6):
            for c3 in range(6):
                blocks.setdefault(int(ray[a3] + ray[b3] + ray[c3]), []).append(
                    (a3, b3, c3))
    masks = {}
    for nb, cells in blocks.items():
        m = np.zeros((6, 6, 6), dtype=bool)
        for c in cells:
            m[c] = True
        masks[nb] = m

    print(f"{'gamma':>8} {'resid med':>10} {'modes cut':>10} "
          f"{'d_min/d_max':>12} | relative round-trip error per block")
    print(f"{'':>8} {'':>10} {'':>10} {'':>12} | "
          f"{'3+0':>10} {'2+1':>10} {'1+2':>10} {'0+3':>10}")

    rows = []
    for gamma in a.gammas:
        cg = float(np.cos(np.deg2rad(gamma / 60.0)))
        resid, cut, cond, lost, blockres = [], [], [], [], []
        for k in idx:
            lk = float(lam_full[k])
            S2 = mc._assemble_6x6(builder, cg, lk)
            Vk = S2 / (2.0 * sig)
            Z6 = ds.zeta6(cg, lk) / (2.0 * sig) ** 2
            Q = ds.solve_Q(Vk, Z6, rcond=a.rcond)
            Zb = ds.cum3_from_Q(Q, Vk)
            nz = np.linalg.norm(Z6)
            resid.append(np.linalg.norm(Zb - Z6) / max(nz, 1e-300))
            d = np.linalg.eigvalsh(Vk)
            cut.append(int(np.sum(d <= a.rcond * d.max())))
            cond.append(float(d.min() / d.max()))
            # how much of the cumulant NORM sits in the directions being zeroed
            good = d > (a.rcond * d.max())
            U = np.linalg.eigh(Vk)[1]
            Zp = np.einsum("ap,bq,cr,abc->pqr", U, U, U, Z6, optimize=True)
            keep = good[:, None, None] & good[None, :, None] & good[None, None, :]
            lost.append(np.linalg.norm(np.where(keep, 0.0, Zp)) / max(nz, 1e-300))
            blockres.append({
                nb: float(np.linalg.norm((Zb - Z6)[m])
                          / max(np.linalg.norm(Z6[m]), 1e-300))
                for nb, m in masks.items()})
        resid = np.array(resid); lost = np.array(lost)
        bl = {nb: float(np.median([r[nb] for r in blockres])) for nb in masks}
        print(f"{gamma:8.1f} {np.median(resid):10.3e} {np.mean(cut):10.2f} "
              f"{np.median(cond):12.3e} | "
              + " ".join(f"{bl[nb]:10.3e}" for nb in (0, 1, 2, 3)))
        rows.append((gamma, np.median(resid), resid.max(), np.mean(cut),
                     np.median(cond), np.median(lost),
                     bl[0], bl[1], bl[2], bl[3]))

    if a.out:
        arr = np.array(rows)
        np.savez(a.out, gamma=arr[:, 0], resid_med=arr[:, 1], resid_max=arr[:, 2],
                 modes_cut=arr[:, 3], cond=arr[:, 4], norm_lost=arr[:, 5],
                 blk30=arr[:, 6], blk21=arr[:, 7], blk12=arr[:, 8],
                 blk03=arr[:, 9],
                 sigma=sig, rcond=a.rcond)
        print(f"\n[probe] -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
