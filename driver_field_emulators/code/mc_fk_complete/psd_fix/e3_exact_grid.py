"""Exact expectation of the placement-complete estimator, baseline vs R1Knots,
on the paper's gamma grid at both regulator values, plus the paper's
paired sigma_lambda -> 0 extrapolation ``2 x(sigma=4) - x(sigma=8)``.
"""
from __future__ import annotations
import argparse, sys, time
import numpy as np

HERE = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete"
sys.path.insert(0, HERE)
sys.path.insert(0, HERE + "/psd_fix")

import fk_expect_exact as ex
import sigma2_repaired as R
from run_gate2 import analytic_fk

BUILDERS = {"base": None, "R1": R.R1Knots}


def run(gam, sig, n_lambda, apply_c0, calibrate, builder_cls):
    cosg = float(np.cos(np.deg2rad(gam / 60.0)))
    grid = ex.build_grid(n_lambda=n_lambda, apply_anchor=apply_c0)
    if builder_cls is not None:
        grid.builder = builder_cls(background=grid.bg, apply_c0=apply_c0)
    V, A, Q, rho = ex.node_stats(grid, cosg, sig, calibrate=calibrate)
    T1, T2 = ex.expectation(grid, cosg, sig, V=V, A=A, Q=Q, rho=rho)
    tot = (T1 + T2) + (T1 + T2).T
    vr = T1 + T1.T
    return float(tot[0, 3]), float(vr[0, 3])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 2.6, 5.0, 6.7, 17.3])
    ap.add_argument("--sigmas", type=float, nargs="+", default=[8.0, 4.0])
    ap.add_argument("--n-lambda", type=int, default=1000)
    ap.add_argument("--apply-c0", type=int, default=0)
    ap.add_argument("--calibrate", default="smeared")
    ap.add_argument("--out", default=HERE + "/psd_fix/_e3_exact_grid.npz")
    a = ap.parse_args()

    print(f"N={a.n_lambda} apply_c0={bool(a.apply_c0)} calib={a.calibrate}")
    res = {}
    for gam in a.gammas:
        fold = float(analytic_fk([gam])[0])
        line = [f"g={gam:>6.2f} fold={fold:.6e}"]
        for sig in a.sigmas:
            for tag, cls in BUILDERS.items():
                t0 = time.time()
                tot, vr = run(gam, sig, a.n_lambda, bool(a.apply_c0),
                              a.calibrate, cls)
                res[(gam, sig, tag)] = tot
                line.append(f"  s{sig:g}/{tag}={tot:.7e}({tot/fold:.5f})"
                            f"[{time.time()-t0:.0f}s]")
        print(" ".join(line), flush=True)
        res[(gam, "fold", "-")] = fold

    print()
    hdr = (f"{'gamma':>7} {'fold':>13} | {'s8 base':>12} {'s8 R1':>12} {'R1/base':>9} | "
           f"{'s4 base':>12} {'s4 R1':>12} {'R1/base':>9} | "
           f"{'ext base':>12} {'/fold':>8} {'ext R1':>12} {'/fold':>8} {'R1/base':>9}")
    print(hdr); print("-" * len(hdr))
    rows = []
    for gam in a.gammas:
        fold = res[(gam, "fold", "-")]
        b8, r8 = res[(gam, 8.0, "base")], res[(gam, 8.0, "R1")]
        b4, r4 = res[(gam, 4.0, "base")], res[(gam, 4.0, "R1")]
        eb, er = 2 * b4 - b8, 2 * r4 - r8
        print(f"{gam:>7.2f} {fold:>13.6e} | {b8:>12.6e} {r8:>12.6e} {r8/b8:>9.6f} | "
              f"{b4:>12.6e} {r4:>12.6e} {r4/b4:>9.6f} | "
              f"{eb:>12.6e} {eb/fold:>8.5f} {er:>12.6e} {er/fold:>8.5f} {er/eb:>9.6f}")
        rows.append((gam, fold, b8, r8, b4, r4, eb, er))
    np.savez(a.out, rows=np.array(rows, float))
    print("wrote", a.out)
