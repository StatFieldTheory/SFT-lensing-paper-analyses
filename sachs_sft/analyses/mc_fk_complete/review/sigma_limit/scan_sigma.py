"""sigma_lambda scan of the exact FK estimator expectation, with the
placement-resolved split and the parameter-free O(sigma) contact prediction.

For each (N, sigma) it reports, at full precision:

  ref3   = sum dlam W H F:X3      -- sigma->0 limit of T1 + T1^T  (Q on the
                                     observable leg)
  ref12  = sum dlam W H F:(X1+X2) -- sigma->0 limit of T2 + T2^T  (Q on an
                                     F-vertex input leg)
  ref    = ref3 + ref12           -- must equal fk_reference(zeta_injected)
  T1s, T2s, tot                   -- the exact finite-sigma expectation
  C3, C12 = sum dlam W^2 F:X_i    -- the contact (kink) integrals

The prediction derived in review/DERIVATION is

  T1s(sigma) = ref3  - mu2 * C3   + O(sigma^2),   mu2 -> (3/4) sigma
  T2s(sigma) = ref12 - mu1 * C12  + O(sigma^2),   mu1 -> (1/2) sigma

with mu2 = E[max(U,V)], mu1 = E[max(U,0)] for the (lattice) AR(1) smearing
kernel; the continuum Laplace values are 3/4 and 1/2 times sigma.
"""
from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/"
                   "SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete")
import fk_expect_exact as ex   # noqa: E402


def lattice_moments(rho: float, kmax: int = 20000):
    """E[max(K,L)] and E[max(K,0)] in LATTICE UNITS for p_k ~ rho^|k|."""
    k = np.arange(-kmax, kmax + 1)
    p = rho ** np.abs(k)
    p = p / p.sum()
    mu1 = float((np.maximum(k, 0) * p).sum())
    # E[max(K,L)] = E|K-L|/2 for symmetric K,L
    conv = np.convolve(p, p[::-1])                 # distribution of K-L
    d = np.arange(-2 * kmax, 2 * kmax + 1)
    mu2 = float(0.5 * (np.abs(d) * conv).sum())
    return mu2, mu1


def placement_pieces(grid, Q, B):
    """X1, X2, X3 with index order [p, b, c, B]; X1+X2+X3 = zeta_injected."""
    X3 = np.einsum("pBij,pib,pjc->pbcB", Q, B, B, optimize=True)
    X1 = np.einsum("pbij,pic,pjB->pbcB", Q, B, B, optimize=True)
    X2 = np.einsum("pcij,pib,pjB->pbcB", Q, B, B, optimize=True)
    return X1, X2, X3


def fold(kern, X):
    M = np.einsum("p,Abc,pbcB->AB", kern, ex.F6, X, optimize=True)
    return M + M.T


def one(gam, N, sig, calib="smeared"):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    t0 = time.time()
    grid = ex.build_grid(n_lambda=N)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
    B = ex.calib_matrix(grid, V, A, rho, sig, calib)
    zi = ex.zeta_injected(B, Q)
    X1, X2, X3 = placement_pieces(grid, Q, B)
    # consistency: X1+X2+X3 must equal the injected cumulant
    resid = float(np.abs(X1 + X2 + X3 - np.transpose(zi, (0, 1, 2, 3))).max()
                  / np.abs(zi).max())
    kern_ref = grid.dlam * grid.Wd * grid.Hd
    kern_con = grid.dlam * grid.Wd * grid.Wd
    ref3 = fold(kern_ref, X3)[0, 3]
    ref12 = fold(kern_ref, X1 + X2)[0, 3]
    C3 = fold(kern_con, X3)[0, 3]
    C12 = fold(kern_con, X1 + X2)[0, 3]
    ref_full = ex.fk_reference(grid, zi)[0, 3]
    ref_true = ex.fk_reference(grid, np.array(
        [ex._DS.zeta6(cosg, float(l)) for l in grid.lam]))[0, 3]
    blk = 128 if N <= 2000 else 64
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho, block=blk)
    t1s = (T1 + T1.T)[0, 3]
    t2s = (T2 + T2.T)[0, 3]
    mu2, mu1 = lattice_moments(rho)
    return dict(gamma=gam, N=N, sigma=sig, dlam=grid.dlam, rho=rho,
                ratio_dl=sig / grid.dlam, resid=resid,
                ref3=ref3, ref12=ref12, ref=ref3 + ref12, ref_full=ref_full,
                ref_true=ref_true, C3=C3, C12=C12, t1s=t1s, t2s=t2s,
                tot=t1s + t2s,
                mu2_lat=mu2 * grid.dlam, mu1_lat=mu1 * grid.dlam,
                sec=time.time() - t0)


def main():
    gam = float(sys.argv[1])
    spec = sys.argv[2]           # "251:32,500:16,..."
    out = sys.argv[3]
    calib = sys.argv[4] if len(sys.argv) > 4 else "smeared"
    rows = []
    hdr = (f"{'N':>6} {'dlam':>8} {'sigma':>7} {'s/dl':>6} "
           f"{'tot':>15} {'ref':>15} {'ratio':>11} "
           f"{'T1s/ref3':>10} {'T2s/ref12':>10} {'resid':>9} {'s':>5}")
    print(f"# gamma={gam}'  calib={calib}")
    print(hdr, flush=True)
    for item in spec.split(","):
        Ns, ss = item.split(":")
        r = one(gam, int(Ns), float(ss), calib)
        rows.append(r)
        print(f"{r['N']:>6} {r['dlam']:>8.4f} {r['sigma']:>7.2f} "
              f"{r['ratio_dl']:>6.2f} {r['tot']:>15.8e} {r['ref']:>15.8e} "
              f"{r['tot']/r['ref']:>11.7f} {r['t1s']/r['ref3']:>10.6f} "
              f"{r['t2s']/r['ref12']:>10.6f} {r['resid']:>9.1e} "
              f"{r['sec']:>5.0f}", flush=True)
    keys = rows[0].keys()
    np.savez(out, **{k: np.array([r[k] for r in rows]) for k in keys})
    print(f"-> {out}")


if __name__ == "__main__":
    main()
