#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""PHASE-1 ORACLE for the HIGH / flat-sky Born-Limber spin-2 orientation kernel.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

WHAT THIS PINS
==============
The deployed cosmic-shear zeta_D (spins (2,2,-2)) is 4-17x TOO LARGE vs the
trusted B_kappa-phase reference (= fastnc to ~1%), GROWING with opening angle phi.
The prior LOW-kernel fix (SpinAwareZetaEvaluator) did NOT move it: zeta_D is
HIGH-dominated by ~8 orders of magnitude, and the HIGH path uses an INDEPENDENT
flat-sky orientation integral `_kappa3_limber_alpha_kernel` (kappa3.py:1747).

This probe drills into `_kappa3_limber_alpha_kernel` and compares it, term by
term, to the analytic flat-sky spin-2 orientation integral that the trusted
reference uses, at the kernel level (fixed triangle, fixed (u,v,phi)) AND folded
to the shear Gamma.  No factor is fit to anything; the reference orientation
integral is DERIVED from the same definition as stage1_bkappa_phase_reference.py.

THE ANALYTIC ORIENTATION INTEGRAL (derived; see scripts/docs/zetaD_HIGH_fix.md)
------------------------------------------------------------------------------
Flat-sky natural component for spins (s1,s2,s3) (canoes uses s_j in {0,+-2}; the
reference uses spin-weight 2*sigma_j with sigma_j in {+-1}, so s_j == 2 sigma_j):

  Gamma = INT d^2 l2 d^2 l3 /(2pi)^4 B(l1,l2,l3)
            * exp[i (l2.d2 + l3.d3)]
            * exp[i (s1 phi_{l1} + s2 phi_{l2} + s3 phi_{l3})]      (l1 = -(l2+l3))

Parametrize l2 = u e^{i a}, l3 = v e^{i(a+phi)}, l1 = -(u + v e^{i phi}) e^{i a}.
Integrating over the absolute orientation a in [0,2pi):

  I(s1,s2,s3; u,v,phi) = INT_0^{2pi} da exp[i Re(e^{-i a} Q)] exp[i(S a + s1 psi + s3 phi)]
                       = 2 pi i^S J_S(|Q|) e^{i S beta} e^{i(s1 psi + s3 phi)}

  Q   = u r2 + v e^{-i phi} r3   (r2 = d2, r3 = d3 as complex 2D vectors)
  S   = s1 + s2 + s3
  beta= arg(Q)
  psi = arg(l1_rel),  l1_rel = -(u + v e^{i phi})

This I is COMPLEX.  Gamma = INT (real radial measure) * B * I, so Gamma is complex.

THE SUSPECT
-----------
canoes `_kappa3_limber_alpha_kernel` returns `np.real(I)` per cell for spin_sum!=0
(line 1785), i.e. it folds Re(I) over the radial measure, producing a REAL zeta_D,
then |Gamma| = |fold(Re(I)...)|.  The correct |Gamma| = |fold(I_complex...)|.
This probe measures, cell-by-cell and folded, canoes-real vs analytic-complex vs
a brute-force a-integral (both the canoes-Q version and the reference-definition
version), for (2,2,-2) at the PRIMARY (10',60) and WORST (10',120) configs and a
phi sweep, to pin which term is wrong and whether it reproduces the 4-17x growth.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
from scipy.special import jv

os.environ.setdefault("JAX_PLATFORMS", "cpu")

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

ARCMIN2RAD = math.pi / 10800.0

# Import the EXACT deployed canoes kernel + geometry helper (current branch).
import sys
_CANOES = Path("/Users/zzhang/projects/canoes")
if str(_CANOES / "src") not in sys.path:
    sys.path.insert(0, str(_CANOES / "src"))
from canoes.sachs.kappa3 import (  # noqa: E402
    _kappa3_limber_alpha_kernel,
    _kappa3_flat_triangle_vertices,
)


# ---------------------------------------------------------------------------
# 1. Analytic complex orientation integral I (the DERIVED truth)
# ---------------------------------------------------------------------------
def analytic_I(spins, u, v, phi, w, r2, r3):
    """Full COMPLEX I = 2pi i^S J_S(|Q|) e^{iS beta} e^{i(s1 psi + s3 phi)}."""
    s1, s2, s3 = spins
    S = int(s1 + s2 + s3)
    Q = u * r2 + v * np.exp(-1j * phi) * r3
    qabs = np.abs(Q)
    beta = np.where(qabs > 0.0, np.angle(Q), 0.0)
    l1_rel = -(u + v * np.exp(1j * phi))
    psi = np.where(w > 0.0, np.angle(l1_rel), 0.0)
    return (2.0 * math.pi * (1j ** S) * np.exp(1j * S * beta)
            * jv(S, qabs) * np.exp(1j * (s1 * psi + s3 * phi)))


# ---------------------------------------------------------------------------
# 2. Brute-force a-integral with the CANOES Q-geometry (settles analytic algebra)
# ---------------------------------------------------------------------------
def brute_I_canoesgeom(spins, u, v, phi, r2, r3, n_a=4096):
    """INT_0^{2pi} da exp[i Re(e^{-i a} Q)] exp[i(S a + s1 psi + s3 phi)] (complex).

    Built from canoes' OWN (Q, l1_rel) algebra, so it isolates whether the
    closed-form Bessel reduction is right (vs the numerical a-integral)."""
    s1, s2, s3 = spins
    S = int(s1 + s2 + s3)
    a = np.linspace(0.0, 2.0 * np.pi, n_a, endpoint=False)
    da = 2.0 * np.pi / n_a
    Q = u * r2 + v * np.exp(-1j * phi) * r3
    l1_rel = -(u + v * np.exp(1j * phi))
    psi = np.angle(l1_rel) if abs(l1_rel) > 0 else 0.0
    pos = np.real(np.exp(-1j * a) * Q)                      # l2.d2 + l3.d3
    integ = np.exp(1j * pos) * np.exp(1j * (S * a + s1 * psi + s3 * phi))
    return np.sum(integ) * da


# ---------------------------------------------------------------------------
# 3. Brute-force a-integral from the REFERENCE DEFINITION (trusted, independent)
# ---------------------------------------------------------------------------
def ref_triangle_geometry(gamma_rad, phi_rad):
    """SAS isoceles vertices + centroid alpha angles, IDENTICAL to
    stage1_bkappa_phase_reference.triangle_geometry."""
    g = gamma_rad
    X3 = np.array([0.0, 0.0])
    X1 = g * np.array([np.cos(phi_rad / 2.0),  np.sin(phi_rad / 2.0)])
    X2 = g * np.array([np.cos(phi_rad / 2.0), -np.sin(phi_rad / 2.0)])
    cen = (X1 + X2 + X3) / 3.0
    a1 = np.arctan2((cen - X1)[1], (cen - X1)[0])
    a2 = np.arctan2((cen - X2)[1], (cen - X2)[0])
    a3 = np.arctan2((cen - X3)[1], (cen - X3)[0])
    d2 = X2 - X1
    d3 = X3 - X1
    return d2, d3, a1, a2, a3


def ref_I(sigma_triplet, u, v, phi, d2, d3, a1, a2, a3, n_a=4096):
    """The reference's orientation integral for spin-weight-2 components.

    sigma_triplet = (s1,s2,s3) in {+-1} (spin weight = 2 sigma).  This rebuilds
    the EXACT integrand of stage1_bkappa_phase_reference.natural_components but
    for a SINGLE (u,v,phi) (B and the radial measure factored out), integrating
    only over the absolute orientation a:

      I_ref = INT da exp[i(l2.d2 + l3.d3)]
                     exp[2 i (s1 (phi_l1 - a1) + s2 (phi_l2 - a2) + s3 (phi_l3 - a3))]

    with l2 = u e^{i a}, l3 = v e^{i(a+phi)}, l1 = -(l2+l3).  Centroid alphas a_j
    are pure constants (the projection is a magnitude-1 phase)."""
    s1, s2, s3 = sigma_triplet
    a = np.linspace(0.0, 2.0 * np.pi, n_a, endpoint=False)
    da = 2.0 * np.pi / n_a
    l2x = u * np.cos(a);            l2y = u * np.sin(a)
    l3x = v * np.cos(a + phi);      l3y = v * np.sin(a + phi)
    l1x = -(l2x + l3x);             l1y = -(l2y + l3y)
    phl1 = np.arctan2(l1y, l1x)
    phl2 = a
    phl3 = a + phi
    pos = l2x * d2[0] + l2y * d2[1] + l3x * d3[0] + l3y * d3[1]
    spin = 2.0 * (s1 * (phl1 - a1) + s2 * (phl2 - a2) + s3 * (phl3 - a3))
    integ = np.exp(1j * pos) * np.exp(1j * spin)
    return np.sum(integ) * da


# ---------------------------------------------------------------------------
# 4. Geometry bridges between canoes (r2,r3) and the reference (d2,d3)
# ---------------------------------------------------------------------------
def sas_thetas(gamma_arcmin, phi_deg):
    """canoes cosine-triple -> (theta_12, theta_23, theta_31) via the SAS map
    used in stage1_ours_v2_FIXED._sas_cosine_triples and _kappa3_*_high:
      cos_arr[:,0]=cos(t3), [:,1]=cos(g), [:,2]=cos(g);  t3=2 g sin(phi/2).
      theta_12=arccos(cos_arr[0])=t3, theta_23=g, theta_31=g."""
    g = gamma_arcmin * ARCMIN2RAD
    ph = math.radians(phi_deg)
    t3 = 2.0 * g * math.sin(ph / 2.0)
    return t3, g, g  # theta_12, theta_23, theta_31


def main():
    print("=" * 96)
    print("PHASE-1 HIGH/Limber ORACLE: canoes _kappa3_limber_alpha_kernel (2,2,-2)")
    print("  vs analytic complex I, brute(canoes Q), brute(reference definition)")
    print("=" * 96)

    # canoes spins (2,2,-2) == reference sigma (+1,+1,-1) == Gamma^3 (s=+,+,-)
    can_spins = (2, 2, -2)
    ref_sigma = (1, 1, -1)

    configs = [
        ("PRIMARY g=10' phi=60", 10.0, 60.0),
        ("WORST   g=10' phi=120", 10.0, 120.0),
        ("        g=10' phi=30", 10.0, 30.0),
        ("        g=50' phi=60", 50.0, 60.0),
    ]
    # representative (u,v,phi_ell) cells from the HIGH Gauss grid (ell in (60,1000])
    cells = [
        (200.0, 200.0, math.radians(60.0)),
        (200.0, 200.0, math.radians(120.0)),
        (500.0, 300.0, math.radians(90.0)),
        (800.0, 800.0, math.radians(50.0)),
    ]

    records = []
    for cfg_name, gam, phid in configs:
        t12, t23, t31 = sas_thetas(gam, phid)
        r2, r3 = _kappa3_flat_triangle_vertices(t12, t23, t31)
        d2, d3, a1, a2, a3 = ref_triangle_geometry(gam * ARCMIN2RAD, math.radians(phid))
        print(f"\n### {cfg_name}")
        print(f"    canoes r2={r2:.4e} r3={r3:.4e}  | ref d2={d2} d3={d3}")
        print(f"    {'(u,v,phi_l)':>22s} {'can_real':>13s} "
              f"{'Re(anaI)':>13s} {'|anaI|':>12s} {'|brute_cQ|':>12s} {'|ref_I|':>12s} "
              f"{'can/|ref|':>10s}")
        for (u, v, phil) in cells:
            ua = np.array([u]); va = np.array([v]); pa = np.array([phil])
            wa = np.sqrt(np.maximum(u * u + v * v + 2 * u * v * np.cos(phil), 0.0))
            can = float(_kappa3_limber_alpha_kernel(
                spins=can_spins, u=ua, v=va, phi=pa, w=np.array([wa]),
                r2=r2, r3=r3)[0])
            aI = analytic_I(can_spins, u, v, phil, wa, r2, r3)
            bI = brute_I_canoesgeom(can_spins, u, v, phil, r2, r3)
            rI = ref_I(ref_sigma, u, v, phil, d2, d3, a1, a2, a3)
            ratio = can / abs(rI) if abs(rI) > 1e-300 else float("nan")
            print(f"    ({u:.0f},{v:.0f},{math.degrees(phil):.0f})".rjust(24)
                  + f" {can:13.4e} {aI.real:13.4e} {abs(aI):12.4e} "
                  + f"{abs(bI):12.4e} {abs(rI):12.4e} {ratio:10.3f}")
            records.append(dict(cfg=cfg_name, gam=gam, phi=phid, u=u, v=v,
                                phil=phil, can=can, aI=aI, bI=bI, rI=rI))

    # ---- KEY: does Re(I) match the reference's COMPLEX I in magnitude? -------
    print("\n" + "=" * 96)
    print("KERNEL-LEVEL VERDICT (per cell):")
    print("  (A) analytic complex I  == brute(canoes Q)?  (closed-form algebra check)")
    print("  (B) |ref_I| (reference def) == |analytic I| (canoes-geom def)?  "
          "(geometry+def agreement)")
    print("  (C) canoes returns Re(I); is Re(I) the right magnitude vs |I|?")
    rel_AB = []
    rel_geom = []
    for r in records:
        # (A) closed form vs brute, canoes geometry
        if abs(r["bI"]) > 1e-300:
            rel_AB.append(abs(r["aI"] - r["bI"]) / abs(r["bI"]))
        # (B) reference-definition |I| vs canoes-geom analytic |I|
        if abs(r["rI"]) > 1e-300:
            rel_geom.append(abs(abs(r["aI"]) - abs(r["rI"])) / abs(r["rI"]))
    print(f"  (A) max rel-err analytic-vs-brute(canoesQ): {max(rel_AB):.2e}")
    print(f"  (B) max rel-err |ref_I| vs |analytic I|:     {max(rel_geom):.2e}")
    print("  (C) per-cell Re(I)/|I| and Re(I)/|ref_I|:")
    for r in records[:8]:
        reI = r["aI"].real
        print(f"      {r['cfg'][:18]:18s} u={r['u']:.0f} v={r['v']:.0f} "
              f"phi_l={math.degrees(r['phil']):.0f}: "
              f"Re(I)/|I|={reI/abs(r['aI']):+.4f}  "
              f"Re(I)/|ref_I|={reI/abs(r['rI']):+.4f}  "
              f"|I|/|ref_I|={abs(r['aI'])/abs(r['rI']):.4f}")

    np.savez(
        OUT / "HIGH_limber_oracle.npz",
        cfg=np.array([r["cfg"] for r in records]),
        gam=np.array([r["gam"] for r in records]),
        phi=np.array([r["phi"] for r in records]),
        u=np.array([r["u"] for r in records]),
        v=np.array([r["v"] for r in records]),
        phil=np.array([r["phil"] for r in records]),
        can_real=np.array([r["can"] for r in records]),
        analytic_I=np.array([r["aI"] for r in records], dtype=np.complex128),
        brute_canoesQ=np.array([r["bI"] for r in records], dtype=np.complex128),
        ref_I=np.array([r["rI"] for r in records], dtype=np.complex128),
        note=("HIGH/Limber spin-2 orientation kernel oracle for (2,2,-2). "
              "can_real = deployed _kappa3_limber_alpha_kernel (np.real(I)). "
              "analytic_I = derived complex 2pi i^S J_S(|Q|) e^{iS beta} "
              "e^{i(s1 psi + s3 phi)}. brute_canoesQ = numeric a-integral on "
              "canoes Q-geometry. ref_I = numeric a-integral from the "
              "stage1_bkappa_phase_reference DEFINITION (centroid spin phases). "
              "All COMPLEX; canoes keeps only Re."),
    )
    print(f"\nsaved -> {OUT / 'HIGH_limber_oracle.npz'}")


if __name__ == "__main__":
    main()
