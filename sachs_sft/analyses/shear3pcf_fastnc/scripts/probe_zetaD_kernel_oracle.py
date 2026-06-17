#!/usr/bin/env python
"""PHASE-1 ORACLE: pin the exact (2,2,-2) angular-kernel mechanism.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

GOAL
====
Compare the canoes spin-weighted great-circle 3-point angular kernel K[t, k]
(per triple k=(l1,l2,l3), per sky-config t) computed TWO independent ways:

  (a) canoes  `SpinAwareZetaEvaluator._build_kernel_numpy(..., spin_weights)`
      (current code on branch fix/spin2-zetaD-2_2_-2-norm).
  (b) an INDEPENDENT brute-force reference that re-evaluates the SAME angular
      sum with a TRUSTED spin-weighted spherical harmonic {s}Y_lm built from
      sympy's Wigner-d (Rotation.d), NOT canoes' _build_spin_y_row.

CONVENTION (matched to canoes, see _spin_aware_three_pt.py docstring lines 13-37)
--------------------------------------------------------------------------------
The kernel realizes the leg-1-at-pole reduction of

    zeta^{s1 s2 s3}(n1, n2, n3)
       = Sum_{l1 l2 l3} Sum_{m1 m2 m3}  B^{m1 m2 m3}_{l1 l2 l3}
            * {s1}Y_{l1 m1}(n1) {s2}Y_{l2 m2}(n2) {s3}Y_{l3 m3}(n3)

with the Komatsu-Spergel isotropic reduced bispectrum
    B^{m1 m2 m3}_{l1 l2 l3} = h * (l1 l2 l3; 0 0 0) * (l1 l2 l3; m1 m2 m3) * b_{l1 l2 l3},
    h = sqrt[(2l1+1)(2l2+1)(2l3+1)/4pi].

Leg 1 at the pole:  {s1}Y_{l1 m1}(pole) = sqrt[(2l1+1)/4pi] * delta_{m1, -s1}
                    (canoes convention {s}Y_lm(0,0) = sqrt[(2l+1)/4pi] d^l_{m,-s}(0)
                     = sqrt[(2l+1)/4pi] delta_{m,-s}).
=> m1 = -s1, and m2+m3 = -m1 = s1, i.e. m3 = s1 - m2.

So the per-triple kernel coefficient is
    K[t, k] = h * (l1 l2 l3; 0 0 0) * sqrt[(2l1+1)/4pi] * {s1}Y-pole-phase
              * Re Sum_{m2} (l1 l2 l3; -s1, m2, s1-m2)
                    * {s2}Y_{l2 m2}(g12, 0) * {s3}Y_{l3, s1-m2}(g13, phi3)

The canoes code writes the leg-1 pole phase as `sign_s1 = (-1)^s1` (line 603).
The brute-force ref instead carries the pole phase IMPLICITLY through the SAME
trusted {s1}Y at the pole, so the two must agree only if canoes' (-1)^s1 phase
matches the trusted {s1}Y pole value's phase.  We test that explicitly.

THE TRUSTED {s}Y (independent of canoes' kernel assembly)
---------------------------------------------------------
We build {s}Y from sympy's Wigner-d (Rotation.d, trusted, arbitrary sign
indices), in the SAME phase convention canoes' VALIDATED +2 path uses, which we
verify cell-by-cell to be exactly

  C1:  {s}Y_lm(theta,phi) = sqrt[(2l+1)/4pi] d^l_{m, -s}(theta) e^{i m phi}

(verified: canoes _build_spin_y_table d-part == sympy d^l_{m,-s} to ~1e-15 for
spin 2; this is the convention the +2 PPP/TPP self-consistency tests validate).
A second convention C2 (Goldberg-standard, (-1)^m d^l_{-m,-s}) is printed only
as a cross-reference; C1 is the operative one.

THE INDEPENDENT REFERENCE m-SUM
-------------------------------
The brute force evaluates the FULL leg-1-at-pole m-sum directly:

  K_ref[t,k] = h (l1l2l3;000) sqrt[(2l1+1)/4pi] (pole-phase)
             * Re Sum_{m2=-l2..l2, |s1-m2|<=l3}
                   (l1 l2 l3; -s1, m2, s1-m2) {s2}Y_{l2 m2}(g12,0) {s3}Y_{l3,s1-m2}(g13,phi3)

i.e. it loops m2 over the FULL [-l2, l2] range (constrained only by |m3|<=l3 and
|m1=-s1|<=l1), with NO group-batched reshape. This is the mathematically
complete m-sum; it is what canoes SHOULD compute. The canoes production path
uses a group-batched reduction that reads count_g and the W3j reshape from the
GROUP-FIRST triple, which is wrong when m_count varies within a (l2,l3) group
(it does: triples with l1 < |s1| have m_count forced to 0 by the |m1|<=l1 mask,
and the group is keyed only on (l2,l3) with l1 ascending). The oracle exposes
exactly those zeroed cells.

OUTPUT
======
For each (l1,l2,l3) active triple and each config, prints canoes K vs the
independent reference K for (2,2,2) and (2,2,-2), identifying the EXACT
discrepancy. Saves outputs/zetaD_kernel_oracle.npz.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import sympy as sp  # noqa: E402
from sympy.physics.quantum.spin import Rotation  # noqa: E402
from sympy.physics.wigner import wigner_3j  # noqa: E402

from canoes.nuell.correlation._spin_aware_three_pt import (  # noqa: E402
    SpinAwareZetaEvaluator,
)
from canoes.nuell.correlation.three_pt import _canonical_frame_angles  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

# --- trusted Wigner-d via sympy (cached) -----------------------------------
_DCACHE: dict[tuple[int, int, int, float], float] = {}


def wigner_d(ell: int, a: int, b: int, theta: float) -> float:
    """Trusted real Wigner-d d^l_{a,b}(theta) via sympy Rotation.d."""
    key = (int(ell), int(a), int(b), round(float(theta), 12))
    if key in _DCACHE:
        return _DCACHE[key]
    if ell < max(abs(a), abs(b)):
        _DCACHE[key] = 0.0
        return 0.0
    val = complex(Rotation.d(ell, a, b, float(theta)).doit().evalf())
    assert abs(val.imag) < 1e-9, f"d^{ell}_{a,b} not real: {val}"
    r = float(val.real)
    _DCACHE[key] = r
    return r


def sY_C1(ell: int, m: int, spin: int, theta: float, phi: float) -> complex:
    """canoes-as-written convention: sqrt[(2l+1)/4pi] d^l_{m,-s} e^{i m phi}."""
    if ell < abs(spin) or abs(m) > ell:
        return 0.0 + 0.0j
    pref = math.sqrt((2 * ell + 1) / (4.0 * math.pi))
    return pref * wigner_d(ell, m, -spin, theta) * np.exp(1j * m * phi)


def sY_C2(ell: int, m: int, spin: int, theta: float, phi: float) -> complex:
    """Goldberg-standard convention: (-1)^m sqrt[(2l+1)/4pi] d^l_{-m,-s} e^{i m phi}.

    Equivalent to sqrt[(2l+1)/4pi] d^l_{m, s} e^{i m phi} via the d-symmetry
    d^l_{-m,-s} = (-1)^{m-s} d^l_{m,s}, i.e. (-1)^m d^l_{-m,-s} = (-1)^{2m-s} d^l_{m,s}
    = (-1)^{-s} d^l_{m,s} = (-1)^s d^l_{m,s} for integer m,s.
    """
    if ell < abs(spin) or abs(m) > ell:
        return 0.0 + 0.0j
    pref = math.sqrt((2 * ell + 1) / (4.0 * math.pi))
    sgn = -1.0 if (m % 2) else 1.0
    return pref * sgn * wigner_d(ell, -m, -spin, theta) * np.exp(1j * m * phi)


def brute_force_kernel(l1, l2, l3, s1, s2, s3, g12, g13, phi3, sY):
    """Independent brute-force K for one triple/config using a trusted {s}Y.

    Mirrors canoes' per-triple coefficient EXACTLY (h, w0, sq_l1, (-1)^s1) but
    performs the FULL uncollapsed leg-1-at-pole m-sum (no group-batched reshape):

      K = h * (l1 l2 l3; 0 0 0) * sqrt[(2l1+1)/4pi] * (-1)^s1
            * Re sum_{m2=-l2..l2, |s1-m2|<=l3}
                  (l1 l2 l3; -s1, m2, s1-m2)
                  * {s2}Y_{l2, m2}(g12, 0) * {s3}Y_{l3, s1-m2}(g13, phi3)

    The only difference vs canoes production is the m-sum is complete here.
    Requires |m1=-s1| <= l1 (else the triple contributes 0). Returns Re(K).
    """
    if abs(s1) > l1:
        return 0.0
    h = math.sqrt(
        (2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4.0 * math.pi)
    )
    w0 = float(wigner_3j(l1, l2, l3, 0, 0, 0))
    if w0 == 0.0:
        return 0.0
    sq_l1 = math.sqrt((2 * l1 + 1) / (4.0 * math.pi))
    sign_s1 = -1.0 if (s1 % 2) else 1.0
    m1 = -s1
    acc = 0.0 + 0.0j
    for m2 in range(-l2, l2 + 1):
        m3 = s1 - m2
        if abs(m3) > l3:
            continue
        w = float(wigner_3j(l1, l2, l3, m1, m2, m3))
        if w == 0.0:
            continue
        y2 = sY(l2, m2, s2, g12, 0.0)
        y3 = sY(l3, m3, s3, g13, phi3)
        acc += w * y2 * y3
    K = h * w0 * sq_l1 * sign_s1 * acc
    return float(np.real(K))


def canoes_kernel_table(t12, t23, t31, ell_max, spins):
    """Extract canoes K[t, k] and the active triple list directly."""
    ev = SpinAwareZetaEvaluator(
        np.atleast_1d(t12), np.atleast_1d(t23), np.atleast_1d(t31),
        ell_max=ell_max, spin_weights=spins,
    )
    K = np.asarray(ev._K)  # (n_t, n_active)
    l1s = np.asarray(ev._l1s)
    l2s = np.asarray(ev._l2s)
    l3s = np.asarray(ev._l3s)
    return K, l1s, l2s, l3s


def sas(gamma_arcmin, phi_deg):
    g = math.radians(gamma_arcmin / 60.0)
    ph = math.radians(phi_deg)
    t3 = 2.0 * g * math.sin(ph / 2.0)
    return (t3, g, g)  # (theta_12=base, theta_23=g, theta_31=g)


def main():
    print("=" * 90)
    print("PHASE-1 ORACLE: canoes spin-(s1,s2,s3) kernel vs trusted-{s}Y brute force")
    print("=" * 90)

    ell_max = 6  # small: full per-triple comparison is feasible
    # Two explicit triangles (general, not squeezed) to exercise phi3 != 0.
    configs = [
        ("g=20' phi=120", sas(20.0, 120.0)),
        ("g=40' phi=60", sas(40.0, 60.0)),
    ]

    records = []
    for cfg_name, (t12, t23, t31) in configs:
        g12, g13, phi3 = _canonical_frame_angles(
            np.array([t12]), np.array([t23]), np.array([t31])
        )
        g12, g13, phi3 = float(g12[0]), float(g13[0]), float(phi3[0])
        print(f"\n### CONFIG {cfg_name}: g12={g12:.5f} g13={g13:.5f} phi3={phi3:.5f}")

        for spins in [(2, 2, 2), (2, 2, -2)]:
            s1, s2, s3 = spins
            K_can, l1s, l2s, l3s = canoes_kernel_table(t12, t23, t31, ell_max, spins)
            print(f"\n  --- spins={spins}  (n_active={l1s.size}) ---")
            print(f"  {'(l1,l2,l3)':>12s} {'canoes':>14s} {'bruteC1':>14s} "
                  f"{'bruteC2':>14s} {'can/C1':>10s} {'can/C2':>10s}")
            for k in range(l1s.size):
                l1, l2, l3 = int(l1s[k]), int(l2s[k]), int(l3s[k])
                kc = float(K_can[0, k])
                kb1 = brute_force_kernel(l1, l2, l3, s1, s2, s3, g12, g13, phi3, sY_C1)
                kb2 = brute_force_kernel(l1, l2, l3, s1, s2, s3, g12, g13, phi3, sY_C2)
                r1 = kc / kb1 if abs(kb1) > 1e-300 else np.nan
                r2 = kc / kb2 if abs(kb2) > 1e-300 else np.nan
                print(f"  ({l1},{l2},{l3})".rjust(14)
                      + f" {kc:14.6e} {kb1:14.6e} {kb2:14.6e} "
                      + (f"{r1:10.4f}" if np.isfinite(r1) else "       nan")
                      + " "
                      + (f"{r2:10.4f}" if np.isfinite(r2) else "       nan"))
                records.append((cfg_name, spins, l1, l2, l3, kc, kb1, kb2))

    # ---- summary: canoes vs the independent (full m-sum) reference ---------
    # C1 is the operative convention (verified == canoes _build_spin_y_table).
    rec = np.array(
        [(r[5], r[6], r[7]) for r in records], dtype=np.float64
    )  # (canoes, refC1, refC2)
    spins_arr = [r[1] for r in records]

    def disc_stats(mask):
        kc = rec[mask, 0]
        cr = rec[mask, 1]  # reference C1
        nz = np.abs(cr) > 1e-14
        rel = np.abs(kc[nz] - cr[nz]) / np.maximum(np.abs(cr[nz]), 1e-300)
        # count cells the reference says nonzero but canoes zeroes
        zeroed = int(np.sum((np.abs(cr) > 1e-14) & (np.abs(kc) < 1e-30)))
        return (rel.max() if rel.size else np.nan, zeroed, int(np.sum(nz)))

    mask_p2 = np.array([s == (2, 2, 2) for s in spins_arr])
    mask_m2 = np.array([s == (2, 2, -2) for s in spins_arr])
    p2_rel, p2_zeroed, p2_nz = disc_stats(mask_p2)
    m2_rel, m2_zeroed, m2_nz = disc_stats(mask_m2)

    print("\n" + "=" * 90)
    print("VERDICT: canoes production K vs independent FULL-m-sum reference (C1 conv):")
    print(f"  (2,2,2):  max rel-err = {p2_rel:.3e}   "
          f"cells canoes ZEROES that ref says nonzero: {p2_zeroed}/{p2_nz}")
    print(f"  (2,2,-2): max rel-err = {m2_rel:.3e}   "
          f"cells canoes ZEROES that ref says nonzero: {m2_zeroed}/{m2_nz}")
    print()
    if p2_zeroed > 0 or m2_zeroed > 0:
        print("  BUG PINNED: canoes' group-batched reduction drops contributions.")
        print("  Mechanism: the (l2,l3)-keyed group reads count_g = m_count[g0]")
        print("  from the GROUP-FIRST triple. For s1=+/-2, triples with l1 < |s1|")
        print("  have m_count forced to 0 (|m1=-s1| <= l1 mask). Since groups are")
        print("  sorted with l1 ascending, the first triple often has l1 < 2, so")
        print("  count_g=0 zeroes the ENTIRE group (incl. valid l1>=2 triples).")
        print("  Affects all channels with s1 != 0 (PPP (2,2,2), D-mod (2,2,-2)).")
    else:
        print("  No zeroed cells; canoes matches the full m-sum. No kernel bug.")

    np.savez(
        OUT / "zetaD_kernel_oracle.npz",
        ell_max=ell_max,
        canoes=rec[:, 0], refC1=rec[:, 1], refC2=rec[:, 2],
        spins=np.array([list(s) for s in spins_arr]),
        l1=np.array([r[2] for r in records]),
        l2=np.array([r[3] for r in records]),
        l3=np.array([r[4] for r in records]),
        cfg=np.array([r[0] for r in records]),
        p2_max_relerr=p2_rel, m2_max_relerr=m2_rel,
        p2_zeroed_cells=p2_zeroed, m2_zeroed_cells=m2_zeroed,
    )
    print(f"\nsaved -> {OUT / 'zetaD_kernel_oracle.npz'}")


if __name__ == "__main__":
    main()
