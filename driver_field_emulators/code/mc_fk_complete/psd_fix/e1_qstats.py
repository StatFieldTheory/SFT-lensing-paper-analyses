"""Downstream effect of the Sigma2 repair, part 1.

(1) deformation-tensor magnitude |Q| per node, and
(2) how often ``sachs_mc_core._cholesky_psd`` has to take its eigen-floor branch,

for the stock ``Sigma2Builder`` (baseline) against ``R1Knots`` (repaired), with
``apply_c0`` matched.  The builder is swapped by patching ``grid.builder`` after
``build_grid``, which is the single place ``node_stats`` reaches the covariance.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np

HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, HERE)
sys.path.insert(0, HERE + "/psd_fix")

import fk_expect_exact as ex
import sigma2_repaired as R

_DS, _CORE = ex._DS, ex._CORE
_CHOL_ORIG = _CORE._cholesky_psd


class CholCounter:
    """Runtime wrapper around ``_cholesky_psd``; the original file is untouched."""

    def __init__(self):
        self.calls = 0
        self.fallbacks = 0
        self.min_eig = []
        self.fallback_idx = []

    def __call__(self, M):
        self.calls += 1
        try:
            return np.linalg.cholesky(M)
        except np.linalg.LinAlgError:
            self.fallbacks += 1
            self.fallback_idx.append(self.calls - 1)
            w, V = np.linalg.eigh(M)
            self.min_eig.append(float(w.min()))
            w = np.maximum(w, 0.0)
            return V * np.sqrt(w)[None, :]


def qstats(gam, sig, n_lambda, apply_c0, tag, builder_cls, calibrate="smeared"):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=n_lambda, apply_anchor=apply_c0)
    if builder_cls is not None:
        grid.builder = builder_cls(background=grid.bg, apply_c0=apply_c0)
    t0 = time.time()
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calibrate)
    # --- Cholesky fallback bookkeeping, exactly the production call pattern ---
    cnt = CholCounter()
    _CORE._cholesky_psd = cnt
    try:
        _ = np.array([_CORE._cholesky_psd(V[k]) for k in range(grid.lam.size)])
    finally:
        _CORE._cholesky_psd = _CHOL_ORIG
    # --- per-node deformation magnitude --------------------------------------
    qn = np.max(np.abs(Q).reshape(Q.shape[0], -1), axis=1)     # (N,)
    med = float(np.median(qn))
    n10 = int((qn > 10.0 * med).sum())
    n100 = int((qn > 100.0 * med).sum())
    eig = np.array([np.linalg.eigvalsh(V[k]) for k in range(V.shape[0])])
    nneg = int((eig.min(axis=1) < 0).sum())
    top = np.argsort(qn)[::-1][:5]
    print(f"  {tag:<10} max|Q|={qn.max():.4e}  median|Q|={med:.4e}  "
          f"ratio={qn.max()/med:8.1f}  n(>10x med)={n10:>4d}  n(>100x med)={n100:>4d}")
    print(f"  {'':<10} V non-PSD nodes={nneg:>4d}/{V.shape[0]}   "
          f"chol fallbacks={cnt.fallbacks:>4d}/{cnt.calls}"
          f"   min eig at fallback={min(cnt.min_eig) if cnt.min_eig else 0.0:+.3e}")
    print(f"  {'':<10} worst-5 nodes lam={np.round(grid.lam[top],1).tolist()} "
          f"|Q|={np.array2string(qn[top], precision=3)}  [{time.time()-t0:.0f}s]")
    return dict(qn=qn, lam=grid.lam, nneg=nneg, nfb=cnt.fallbacks,
                fb_idx=np.array(cnt.fallback_idx), maxQ=float(qn.max()), med=med,
                n10=n10, min_eig=np.array(cnt.min_eig))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+", default=[1.0, 17.3])
    ap.add_argument("--sigma", type=float, default=8.0)
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--apply-c0", type=int, default=0)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--out", default=HERE + "/psd_fix/_e1_qstats.npz")
    a = ap.parse_args()
    store = {}
    for gam in a.gammas:
        print(f"=== gamma = {gam}'  sigma_lambda = {a.sigma}  N = {a.n_lambda}  "
              f"apply_c0 = {bool(a.apply_c0)}  calib = {a.calibrate} ===")
        for tag, cls in (("baseline", None), ("R1Knots", R.R1Knots)):
            d = qstats(gam, a.sigma, a.n_lambda, bool(a.apply_c0), tag, cls, a.calibrate)
            for k, v in d.items():
                store[f"{tag}_{gam}_{k}"] = v
        b, r = store[f"baseline_{gam}_maxQ"], store[f"R1Knots_{gam}_maxQ"]
        print(f"  --> max|Q| baseline/repaired = {b/r:.1f}x   "
              f"(baseline {b:.3e} -> repaired {r:.3e})\n")
    np.savez(a.out, **store)
    print("wrote", a.out)
