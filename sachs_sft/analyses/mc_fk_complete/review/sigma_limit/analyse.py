"""Fit the sigma_lambda exponent and test the parameter-free O(sigma) prediction."""
from __future__ import annotations
import sys
import numpy as np

np.set_printoptions(precision=6, suppress=False)


def load(p):
    d = np.load(p)
    return {k: d[k] for k in d.files}


def show(path, title):
    d = load(path)
    n = d["sigma"].size
    print("=" * 100)
    print(title)
    print("=" * 100)
    print(f"{'N':>6} {'sigma':>7} {'s/dl':>6} {'ratio':>11} {'deficit':>11} "
          f"{'def/sig':>10} {'d1/sig':>10} {'d2/sig':>10} "
          f"{'pred1/sig':>10} {'pred2/sig':>10} {'p1/m1':>7} {'p2/m2':>7}")
    for i in range(n):
        ratio = d["tot"][i] / d["ref"][i]
        dfc = 1.0 - ratio
        d1 = 1.0 - d["t1s"][i] / d["ref3"][i]
        d2 = 1.0 - d["t2s"][i] / d["ref12"][i]
        sig = d["sigma"][i]
        # parameter-free contact prediction, lattice-exact moments
        p1 = d["mu2_lat"][i] * d["C3"][i] / d["ref3"][i]
        p2 = d["mu1_lat"][i] * d["C12"][i] / d["ref12"][i]
        print(f"{int(d['N'][i]):>6} {sig:>7.2f} {d['ratio_dl'][i]:>6.2f} "
              f"{ratio:>11.7f} {dfc:>11.7f} {dfc/sig:>10.6f} "
              f"{d1/sig:>10.6f} {d2/sig:>10.6f} "
              f"{p1/sig:>10.6f} {p2/sig:>10.6f} "
              f"{p1/d1:>7.4f} {p2/d2:>7.4f}")
    return d


def exponents(d, label):
    sig = d["sigma"]
    dfc = 1.0 - d["tot"] / d["ref"]
    o = np.argsort(-sig)
    sig, dfc = sig[o], dfc[o]
    print(f"\n  fitted exponent q in deficit ~ sigma^q  [{label}]")
    for i in range(len(sig) - 1):
        q = np.log(dfc[i] / dfc[i + 1]) / np.log(sig[i] / sig[i + 1])
        print(f"    sigma {sig[i]:>6.2f} -> {sig[i+1]:>6.2f} :  q = {q:.4f}")
    A = np.vstack([np.ones_like(sig), np.log(sig)]).T
    c = np.linalg.lstsq(A, np.log(dfc), rcond=None)[0]
    print(f"    global log-log slope over the whole range: q = {c[1]:.4f}")
    # 3-parameter fit ratio = r0 + a sigma + b sigma^2  -> is r0 = 1 ?
    r = d["tot"] / d["ref"]
    for deg, nm in ((1, "linear"), (2, "quadratic")):
        M = np.vander(d["sigma"], deg + 1, increasing=True)
        co = np.linalg.lstsq(M, r, rcond=None)[0]
        print(f"    {nm} fit of ratio(sigma) over all points: intercept = "
              f"{co[0]:.6f}   (theory: exactly 1)")


for path, title in [(a, a) for a in sys.argv[1:]]:
    d = show(path, title)
    exponents(d, title)
