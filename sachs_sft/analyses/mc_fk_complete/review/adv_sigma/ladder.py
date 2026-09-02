"""Adversarial probe: dense sigma_lambda ladder of the EXACT (noise-free)
expectation, on a user-chosen grid rule.

  --fixed-N   : keep n_lambda fixed (this is the PRODUCTION configuration,
                N = 1000, from which the quoted 0.997/1.002 was extrapolated)
  --fixed-ratio R : choose N so that sigma/dlam = R

Writes an npz with exact, ref_inj (the estimator's own local white-noise
reference built from the INJECTED cumulant), ref_true (same, tabulated zeta)
and the paper fold value.
"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fk_expect_exact as ex

LAM_MIN, LAM_F = 406.0, None


def fold_value(gam_arcmin):
    import _bootstrap
    d = np.load(_bootstrap.FOLD, allow_pickle=True)
    g = d["gamma_arcmin"] if "gamma_arcmin" in d.files else d["gamma"]
    xi = d["xi"]
    # order 2, a=b=0
    v = xi[2, 0, 0] if xi.ndim == 4 else None
    return g, v


def run(gam, sig, N, calib="smeared"):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    t0 = time.time()
    grid = ex.build_grid(n_lambda=N)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calib)
    M = ex.calib_matrix(grid, V, A, rho, sig, calib)
    zi = ex.zeta_injected(M, Q)
    zt = np.array([ex._DS.zeta6(cosg, float(l)) for l in grid.lam])
    ref_inj = ex.fk_reference(grid, zi)[0, 3]
    ref_true = ex.fk_reference(grid, zt)[0, 3]
    blk = 128 if N <= 2000 else 64
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho, block=blk)
    tot = ((T1 + T2) + (T1 + T2).T)[0, 3]
    t1o = (T1 + T1.T)[0, 3]
    return dict(gamma=gam, sigma=sig, N=N, dlam=grid.dlam, ratio=sig / grid.dlam,
                exact=tot, t1=t1o, ref_inj=ref_inj, ref_true=ref_true,
                sec=time.time() - t0)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--sigmas", type=float, nargs="+", required=True)
    ap.add_argument("--fixed-N", type=int, default=0)
    ap.add_argument("--fixed-ratio", type=float, default=0.0)
    ap.add_argument("--calib", default="smeared")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    span = 2313.029 - 406.0
    for s in a.sigmas:
        if a.fixed_N:
            N = a.fixed_N
        else:
            N = int(round(span / (s / a.fixed_ratio))) + 1
        r = run(a.gamma, s, N, a.calib)
        rows.append(r)
        print(f"{r['N']:>6d} dlam={r['dlam']:>7.4f} sig={s:>7.3f} r={r['ratio']:>6.2f} "
              f"exact={r['exact']:.8e} ref_inj={r['ref_inj']:.8e} "
              f"ex/refinj={r['exact']/r['ref_inj']:.8f} [{r['sec']:.0f}s]", flush=True)
    np.savez(a.out, **{k: np.array([r[k] for r in rows]) for k in rows[0]})
    print("->", a.out)
