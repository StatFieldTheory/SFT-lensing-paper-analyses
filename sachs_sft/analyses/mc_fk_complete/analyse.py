"""Assemble the gate results: tables, the sigma_lambda -> 0 extrapolation, plot."""
from __future__ import annotations
import numpy as np
from pathlib import Path


def load(p):
    d = np.load(p, allow_pickle=True)
    return {k: d[k] for k in d.files}


def extrapolate(sig, ratio):
    """Richardson on the linear-in-sigma deficit, plus a least-squares line."""
    o = np.argsort(sig)
    s, r = np.asarray(sig)[o], np.asarray(ratio)[o]
    rich = 2 * r[0] - r[1] if s.size > 1 else np.nan     # two smallest sigma
    A = np.vstack([np.ones_like(s), s]).T
    c, m = np.linalg.lstsq(A, r, rcond=None)[0]
    return rich, c, m


def main():
    print("=" * 78)
    print("GATE 3 -- sigma_lambda dependence at gamma = 1 arcmin")
    print("=" * 78)
    d = load("_gate3_sigma.npz")
    print(f"{'sigma':>7} {'MC':>12} {'seedSEM':>10} {'MC/exact':>9} "
          f"{'exact/ref':>10} {'MC/analytic':>12} {'vr/analytic':>12}")
    for i in np.argsort(-d["sigma"]):
        print(f"{d['sigma'][i]:>7.1f} {d['mc'][i]:>12.5e} {d['sem'][i]:>10.2e} "
              f"{d['mc'][i]/d['exact'][i]:>9.4f} "
              f"{d['exact'][i]/d['ref_inj'][i]:>10.4f} "
              f"{d['mc'][i]/d['analytic'][i]:>12.4f} "
              f"{d['vr'][i]/d['analytic'][i]:>12.4f}")
    for label, num in (("exact / local reference", d["exact"] / d["ref_inj"]),
                       ("MC   / analytic FK", d["mc"] / d["analytic"])):
        rich, c, m = extrapolate(d["sigma"], num)
        print(f"\n  {label}:  linear fit  {c:.4f} + {m:.5f} * sigma"
              f"   |  Richardson (two smallest sigma) -> {rich:.4f}")
    mc_err = d["sem"] / d["analytic"]
    i0 = int(np.argmin(d["sigma"]))
    print(f"  MC seed SEM at the smallest sigma: +/- {mc_err[i0]:.4f}")

    p4 = Path("_gate4_gamma.npz")
    if p4.exists():
        print()
        print("=" * 78)
        print("GATE 4 -- gamma sweep")
        print("=" * 78)
        g = load(p4)
        for sig in sorted(set(g["sigma"]), reverse=True):
            m = g["sigma"] == sig
            print(f"\n  sigma_lambda = {sig}")
            print(f"  {'gamma':>7} {'MC':>12} {'seedSEM':>10} {'analytic':>12} "
                  f"{'MC/ana':>8} {'exact/ana':>10} {'MC/exact':>9} {'vr/ana':>7}")
            for i in np.where(m)[0][np.argsort(g["gamma"][m])]:
                print(f"  {g['gamma'][i]:>7.2f} {g['mc'][i]:>12.5e} "
                      f"{g['sem'][i]:>10.2e} {g['analytic'][i]:>12.5e} "
                      f"{g['mc'][i]/g['analytic'][i]:>8.4f} "
                      f"{g['exact'][i]/g['analytic'][i]:>10.4f} "
                      f"{g['mc'][i]/g['exact'][i]:>9.4f} "
                      f"{g['vr'][i]/g['analytic'][i]:>7.4f}")
        # per-gamma sigma->0 extrapolation from the two sigmas available
        print(f"\n  sigma->0 extrapolation of MC/analytic (linear in sigma, "
              f"from sigma = {sorted(set(g['sigma']))}):")
        print(f"  {'gamma':>7} {'MC/ana(s0)':>11} {'+/-':>8} {'exact/ana(s0)':>14}")
        out = []
        for gam in sorted(set(g["gamma"])):
            m = g["gamma"] == gam
            s = g["sigma"][m]; r = g["mc"][m] / g["analytic"][m]
            e = g["exact"][m] / g["analytic"][m]
            er = g["sem"][m] / g["analytic"][m]
            o = np.argsort(s)
            s, r, e, er = s[o], r[o], e[o], er[o]
            # linear in sigma through the two points -> value at sigma = 0
            w = s[1] / (s[1] - s[0])
            r0 = r[0] * w - r[1] * (w - 1.0)
            e0 = e[0] * w - e[1] * (w - 1.0)
            err0 = np.hypot(er[0] * w, er[1] * (w - 1.0))
            out.append((gam, r0, err0, e0))
            print(f"  {gam:>7.2f} {r0:>11.4f} {err0:>8.4f} {e0:>14.4f}")
        arr = np.array(out)
        print(f"\n  median MC/analytic at sigma->0 over 0.5'-17.3': "
              f"{np.median(arr[:,1]):.4f}   range [{arr[:,1].min():.4f}, "
              f"{arr[:,1].max():.4f}]")
        print(f"  median exact/analytic at sigma->0:                "
              f"{np.median(arr[:,3]):.4f}   range [{arr[:,3].min():.4f}, "
              f"{arr[:,3].max():.4f}]")


if __name__ == "__main__":
    main()
