"""Assemble the verdict: exponent, parameter-free slope, intercept budget.

Theory (see kink_moments.wl):

  T1 (deformation on the far observable leg): BOTH F-legs get smeared, the kink
      argument is Max[U,V],  E[Max[U,V]] = 3/4   ->  deficit_1 = (3/4) sigma C3
  T2 (deformation on an F-vertex input leg):  ONE  F-leg gets smeared, the kink
      argument is Max[U,0],  E[Max[U,0]] = 1/2   ->  deficit_2 = (1/2) sigma C12

  C3  = sum dlam W^2 F:X3 ,   C12 = sum dlam W^2 F:(X1+X2)   (contact integrals)
  ref3, ref12 = sum dlam W H F:X_i                            (sigma->0 limits)
"""
from __future__ import annotations
import numpy as np


def L(p):
    d = np.load(p)
    return {k: d[k] for k in d.files}


def series(d, mask=None):
    m = np.ones(d["sigma"].size, bool) if mask is None else mask
    o = np.argsort(-d["sigma"][m])
    return {k: d[k][m][o] for k in d}


def richardson(dA, dB, p=2.0):
    """lattice-error Richardson: dA at ratio r, dB at ratio 2r, error ~ r^-p."""
    return dB + (dB - dA) / (2.0 ** p - 1.0)


def block(gam, fix, prod):
    a = L(fix)                                     # ratio 4.19
    b = L(prod)                                    # ratio 8.39 subset + prod grid
    hi = series(b, b["ratio_dl"] > 8.0)            # r = 8.38-8.39
    lo = series(a)                                 # r = 4.19-4.20
    pg = series(b, np.isclose(b["N"], 1000))       # production grid N=1000
    print("=" * 104)
    print(f"gamma = {gam} arcmin")
    print("=" * 104)

    print("\n[A] fixed sigma/dlam = 4.19 (the requested scan)")
    print(f"  {'sigma':>7} {'N':>6} {'ratio':>11} {'deficit':>11} {'def/sig':>11} "
          f"{'local exponent q':>18}")
    dl = 1.0 - lo["tot"] / lo["ref"]
    for i in range(dl.size):
        q = "" if i == 0 else (
            f"{np.log(dl[i-1]/dl[i])/np.log(lo['sigma'][i-1]/lo['sigma'][i]):.4f}")
        print(f"  {lo['sigma'][i]:>7.2f} {int(lo['N'][i]):>6} "
              f"{lo['tot'][i]/lo['ref'][i]:>11.7f} {dl[i]:>11.7f} "
              f"{dl[i]/lo['sigma'][i]:>11.7f} {q:>18}")
    A = np.vstack([np.ones_like(lo["sigma"]), np.log(lo["sigma"])]).T
    q = np.linalg.lstsq(A, np.log(dl), rcond=None)[0][1]
    print(f"  global power-law exponent over sigma = 32 -> 2 :  q = {q:.4f}")

    print("\n[B] same sigma at TWICE the resolution (sigma/dlam = 8.39): the "
          "lattice contribution")
    dh = 1.0 - hi["tot"] / hi["ref"]
    # MATCH BY SIGMA, not by position
    sel = np.array([int(np.where(np.isclose(lo["sigma"], s))[0][0])
                    for s in hi["sigma"]])
    dl_m = dl[sel]
    dc = richardson(dl_m, dh)
    print(f"  {'sigma':>7} {'d(r=4.19)':>11} {'d(r=8.39)':>11} "
          f"{'d(extrap)':>11} {'lattice bias at r=4.19':>24}")
    for i in range(dh.size):
        print(f"  {hi['sigma'][i]:>7.2f} {dl_m[i]:>11.7f} {dh[i]:>11.7f} "
              f"{dc[i]:>11.7f} {dl_m[i]/dc[i]-1:>23.4%}")
    print(f"\n  lattice-converged deficit/sigma: "
          + "  ".join(f"{hi['sigma'][i]:.0f}:{dc[i]/hi['sigma'][i]:.6f}"
                      for i in range(dh.size)))
    c1 = 2 * dc[-1] / hi["sigma"][-1] - dc[-2] / hi["sigma"][-2]
    print(f"  -> continuum slope c1 = lim deficit/sigma = {c1:.6f} "
          f"(2-point Richardson in sigma from sigma = 4, 2)")

    print("\n[C] parameter-free prediction of the slope (no fitted parameter)")
    print(f"  {'sigma':>7} {'meas d1/sig':>12} {'pred 3/4':>10} {'p/m':>8} "
          f"{'meas d2/sig':>12} {'pred 1/2':>10} {'p/m':>8} "
          f"{'tot pred/meas':>14}")
    for i in range(dh.size):
        s = hi["sigma"][i]
        m1 = 1.0 - hi["t1s"][i] / hi["ref3"][i]
        m2 = 1.0 - hi["t2s"][i] / hi["ref12"][i]
        p1 = 0.75 * s * hi["C3"][i] / hi["ref3"][i]
        p2 = 0.50 * s * hi["C12"][i] / hi["ref12"][i]
        pt = (0.75 * s * hi["C3"][i] + 0.50 * s * hi["C12"][i]) / hi["ref"][i]
        print(f"  {s:>7.2f} {m1/s:>12.6f} {p1/s:>10.6f} {p1/m1:>8.4f} "
              f"{m2/s:>12.6f} {p2/s:>10.6f} {p2/m2:>8.4f} {pt/dh[i]:>14.4f}")

    print("\n[D] discrimination against the alternatives (at the smallest sigma)")
    i = dh.size - 1
    s = hi["sigma"][i]
    m1 = 1.0 - hi["t1s"][i] / hi["ref3"][i]
    m2 = 1.0 - hi["t2s"][i] / hi["ref12"][i]
    r1, r2 = hi["C3"][i] / hi["ref3"][i], hi["C12"][i] / hi["ref12"][i]
    for name, a1, a2 in [("DERIVED  (3/4 , 1/2)", 0.75, 0.50),
                         ("both 3/4  (wrong)", 0.75, 0.75),
                         ("both 1/2  (wrong)", 0.50, 0.50),
                         ("both 1    (wrong)", 1.00, 1.00),
                         ("no O(sigma) term  (quadratic)", 0.0, 0.0)]:
        print(f"    {name:<32} pred/meas  T1 = {a1*s*r1/m1:>7.4f}   "
              f"T2 = {a2*s*r2/m2:>7.4f}")

    print("\n[E] the intercept: is there anything left at sigma = 0?")
    for nm, dd in (("r = 4.19 series", lo), ("r = 8.39 series", hi)):
        r = dd["tot"] / dd["ref"]
        for deg, lab in ((1, "linear"), (2, "quadratic")):
            M = np.vander(dd["sigma"], deg + 1, increasing=True)
            co = np.linalg.lstsq(M, r, rcond=None)[0]
            print(f"    {nm:<16} {lab:>9} fit, all points : intercept = "
                  f"{co[0]:.6f}   (theory 1 exactly)")
    print("\n[F] the extrapolation the paper actually quotes "
          "(production grid N = 1000, sigma = 8 and 4, 2-point linear)")
    rp = pg["tot"] / pg["ref"]
    sp = pg["sigma"]
    i8 = int(np.where(np.isclose(sp, 8))[0][0])
    i4 = int(np.where(np.isclose(sp, 4))[0][0])
    i2 = int(np.where(np.isclose(sp, 2))[0][0])
    ext84 = 2 * rp[i4] - rp[i8]
    ext42 = 2 * rp[i2] - rp[i4]
    print(f"    ratio(8) = {rp[i8]:.6f}   ratio(4) = {rp[i4]:.6f}  "
          f" ratio(2) = {rp[i2]:.6f}")
    print(f"    2 r(4) - r(8) = {ext84:.6f}      bias = {ext84-1:+.5%}")
    print(f"    2 r(2) - r(4) = {ext42:.6f}      bias = {ext42-1:+.5%}")
    # the same extrapolation with the lattice bias removed
    rc = 1.0 - dc
    print(f"    lattice-converged: 2 r(4)-r(8) = {2*rc[2]-rc[1]:.6f}   "
          f"bias = {2*rc[2]-rc[1]-1:+.5%}   (pure sigma^2 curvature)")
    print(f"                       2 r(2)-r(4) = {2*rc[3]-rc[2]:.6f}   "
          f"bias = {2*rc[3]-rc[2]-1:+.5%}")
    print("\n[G] what the sigma->0 intercept is measured AGAINST")
    print(f"    ref (injected zeta, N=8000)  = {hi['ref'][-1]:.6e}")
    print(f"    ref (tabulated zeta, N=8000) = {hi['ref_true'][-1]:.6e}   "
          f"injection fidelity = {hi['ref'][-1]/hi['ref_true'][-1]:.4f}")
    return dc, hi


for gam, f, p in [(1.0, "review/sigma_limit/_fixratio_g1.npz", "review/sigma_limit/_prod_g1.npz"),
                  (5.0, "review/sigma_limit/_fixratio_g5.npz", "review/sigma_limit/_prod_g5.npz")]:
    block(gam, f, p)
    print()


print("=" * 104)
print("[H] normalisation-free discriminator: the RATIO of the two placement slopes")
print("    (independent of C3, C12 overall scale; tests 3/4 : 1/2 = 3:2 directly)")
print("=" * 104)
for gam, p in [(1.0, "review/sigma_limit/_prod_g1.npz"),
               (5.0, "review/sigma_limit/_prod_g5.npz")]:
    d = L(p)
    hi = series(d, d["ratio_dl"] > 8.0)
    i = hi["sigma"].size - 1
    s = hi["sigma"][i]
    m1 = (1.0 - hi["t1s"][i] / hi["ref3"][i]) / s
    m2 = (1.0 - hi["t2s"][i] / hi["ref12"][i]) / s
    r1, r2 = hi["C3"][i] / hi["ref3"][i], hi["C12"][i] / hi["ref12"][i]
    print(f"  gamma={gam}'  sigma={s}:  measured slope ratio d1/d2 = {m1/m2:.4f}")
    print(f"      DERIVED (3/4 vs 1/2)      predicts {0.75*r1/(0.50*r2):.4f}"
          f"   -> off by {abs(0.75*r1/(0.50*r2)/(m1/m2)-1):.2%}")
    print(f"      EQUAL coefficients        predicts {r1/r2:.4f}"
          f"   -> off by {abs(r1/r2/(m1/m2)-1):.2%}")
    print(f"      REVERSED (1/2 vs 3/4)     predicts {0.50*r1/(0.75*r2):.4f}"
          f"   -> off by {abs(0.50*r1/(0.75*r2)/(m1/m2)-1):.2%}")
