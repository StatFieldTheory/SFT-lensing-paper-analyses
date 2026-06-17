#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""HIGH/Limber fold comparison: canoes orientation kernel vs reference, SAME B.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Isolates the orientation kernel from everything else.  Builds, on ONE common
flat-sky polar grid (u, v, phi_ell) and ONE common B_kappa(l1,l2,l3) Limber
integral, the natural component Gamma^3 (spins (2,2,-2)) folded TWO ways:

  (1) CANOES: Gamma3 = INT (u du)(v dv) dphi / (2pi)^4 * B * canoes_kernel
      with canoes_kernel = np.real( _kappa3_limber_alpha_kernel((2,2,-2)) ).
      (phi runs 0..2pi so the relative-angle integral is complete; this mirrors
      the deployed compute_kappa3_mod_sigma3_high orientation sum.)

  (2) REFERENCE: Gamma3 = INT (u du)(v dv) dphi / (2pi)^4 * B * ref_orient
      with ref_orient = the COMPLEX a-integral of the reference DEFINITION
      (e^{2 i sigma_j (phi_lj - a_j)} centroid phases), divided out of the
      centroid constant so only the orientation content remains.

Both use the SAME B (the reference's fastnc-identical Limber B_kappa) and the
SAME (u,v,phi) measure.  The ONLY difference is the orientation kernel.  So the
ratio |Gamma3_canoes| / |Gamma3_ref| at fixed (gamma,phi_open) is exactly the
HIGH-kernel error, isolated from radial/B/normalization.

This reproduces the 4-17x growing-with-phi factor (or not) using the kernel alone.
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

os.environ.setdefault("JAX_PLATFORMS", "cpu")

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"

# reference B_kappa (fastnc-identical Limber) + geometry
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import (  # noqa: E402
    BkappaLimber, BkappaInterp, triangle_geometry, ARCMIN2RAD,
)

sys.path.insert(0, str(Path("/Users/zzhang/projects/canoes/src")))
from canoes.sachs.kappa3 import (  # noqa: E402
    _kappa3_limber_alpha_kernel, _kappa3_flat_triangle_vertices,
)


def sas_thetas(gamma_arcmin, phi_deg):
    g = gamma_arcmin * ARCMIN2RAD
    ph = math.radians(phi_deg)
    t3 = 2.0 * g * math.sin(ph / 2.0)
    return t3, g, g


def ref_orient_grid(sigma, U, V, PHI, d2, d3, a1, a2, a3, n_a=512):
    """Complex reference orientation integral I_ref(u,v,phi) on a grid.

    INT_0^{2pi} da exp[i(l2.d2+l3.d3)]
                  exp[2 i (s1(phi_l1-a1)+s2(phi_l2-a2)+s3(phi_l3-a3))],
    l2=u e^{i a}, l3=v e^{i(a+phi)}, l1=-(l2+l3).  Vectorised over the (u,v,phi)
    grid; loops a."""
    s1, s2, s3 = sigma
    bshape = np.broadcast_shapes(U.shape, V.shape, PHI.shape)
    acc = np.zeros(bshape, dtype=np.complex128)
    da = 2.0 * np.pi / n_a
    for ia in range(n_a):
        a = ia * da
        l2x = U * math.cos(a);          l2y = U * math.sin(a)
        l3x = V * np.cos(a + PHI);      l3y = V * np.sin(a + PHI)
        l1x = -(l2x + l3x);             l1y = -(l2y + l3y)
        phl1 = np.arctan2(l1y, l1x)
        pos = l2x * d2[0] + l2y * d2[1] + l3x * d3[0] + l3y * d3[1]
        spin = 2.0 * (s1 * (phl1 - a1) + s2 * ((a) - a2) + s3 * ((a + PHI) - a3))
        acc += np.exp(1j * pos) * np.exp(1j * spin)
    return acc * da


def canoes_orient_grid(spins, U, V, PHI, W, r2, r3):
    """canoes np.real orientation kernel on the (u,v,phi) grid (deployed)."""
    return _kappa3_limber_alpha_kernel(spins=spins, u=U, v=V, phi=PHI, w=W,
                                       r2=r2, r3=r3)


def canoes_orient_grid_complex(spins, U, V, PHI, W, r2, r3):
    """Same algebra but KEEP the complex I (analytic), to test the np.real choice."""
    s1, s2, s3 = spins
    S = int(s1 + s2 + s3)
    Q = U * r2 + V * np.exp(-1j * PHI) * r3
    qabs = np.abs(Q)
    beta = np.where(qabs > 0.0, np.angle(Q), 0.0)
    l1_rel = -(U + V * np.exp(1j * PHI))
    psi = np.where(W > 0.0, np.angle(l1_rel), 0.0)
    return (2.0 * math.pi * (1j ** S) * np.exp(1j * S * beta)
            * jv(S, qabs) * np.exp(1j * (s1 * psi + s3 * PHI)))


def main():
    print("=" * 96)
    print("HIGH/Limber FOLD compare: canoes orientation kernel vs reference, SAME B")
    print("=" * 96)

    bk_exact = BkappaLimber(n_chi=120)
    bk = BkappaInterp(bk_exact, n_l=96, n_c=64)

    # flat-sky polar grid (matches the deployed high tail ell range loosely; we
    # integrate the WHOLE plane so the kernel test is clean, not just ell>60).
    n_l, n_phi = 160, 96
    lr = np.geomspace(1.0, 2.0e4, n_l)
    dlnl = np.gradient(np.log(lr))
    wl = lr * lr * dlnl
    pa = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    dpa = 2.0 * np.pi / n_phi
    norm = 1.0 / (2.0 * np.pi) ** 4

    U = lr[:, None, None]
    V = lr[None, :, None]
    PHI = pa[None, None, :]
    W = np.sqrt(np.maximum(U * U + V * V + 2.0 * U * V * np.cos(PHI), 0.0))
    meas = (wl[:, None, None] * wl[None, :, None] * dpa) * norm

    # B_kappa on the grid (depends only on the three magnitudes)
    Bgrid = bk(W, np.broadcast_to(U, W.shape), np.broadcast_to(V, W.shape))

    can_spins = (2, 2, -2)
    ref_sigma = (1, 1, -1)

    configs = [
        ("PRIMARY 10' 60", 10.0, 60.0),
        ("        10' 30", 10.0, 30.0),
        ("WORST   10' 120", 10.0, 120.0),
        ("        50' 60", 50.0, 60.0),
        ("        10' 90", 10.0, 90.0),
        ("squeezed 10' 10", 10.0, 10.0),
        ("        150' 120", 150.0, 120.0),
    ]
    rows = []
    print(f"\n{'config':>18s} {'|G3_can_real|':>14s} {'|G3_can_cplx|':>14s} "
          f"{'|G3_ref|':>12s} {'real/ref':>9s} {'cplx/ref':>9s}")
    for name, gam, phid in configs:
        t12, t23, t31 = sas_thetas(gam, phid)
        r2, r3 = _kappa3_flat_triangle_vertices(t12, t23, t31)
        _X1, _X2, _X3, d2, d3, a1, a2, a3 = triangle_geometry(
            gam * ARCMIN2RAD, math.radians(phid))

        K_can_real = canoes_orient_grid(can_spins, U, V, PHI, W, r2, r3)
        K_can_cplx = canoes_orient_grid_complex(can_spins, U, V, PHI, W, r2, r3)
        K_ref = ref_orient_grid(ref_sigma, U, V, PHI, d2, d3, a1, a2, a3, n_a=256)

        G_can_real = np.sum(meas * Bgrid * K_can_real)
        G_can_cplx = np.sum(meas * Bgrid * K_can_cplx)
        G_ref = np.sum(meas * Bgrid * K_ref)
        rr = abs(G_can_real) / abs(G_ref) if abs(G_ref) > 0 else float("nan")
        rc = abs(G_can_cplx) / abs(G_ref) if abs(G_ref) > 0 else float("nan")
        print(f"{name:>18s} {abs(G_can_real):14.4e} {abs(G_can_cplx):14.4e} "
              f"{abs(G_ref):12.4e} {rr:9.3f} {rc:9.3f}")
        rows.append((name, gam, phid, G_can_real, G_can_cplx, G_ref))

    np.savez(
        OUT / "HIGH_fold_compare.npz",
        names=np.array([r[0] for r in rows]),
        gam=np.array([r[1] for r in rows]),
        phi=np.array([r[2] for r in rows]),
        G_can_real=np.array([r[3] for r in rows], dtype=np.complex128),
        G_can_cplx=np.array([r[4] for r in rows], dtype=np.complex128),
        G_ref=np.array([r[5] for r in rows], dtype=np.complex128),
        note=("Gamma^3 folded over a COMMON (u,v,phi) grid and COMMON B_kappa; "
              "only the orientation kernel differs. can_real=np.real(I) "
              "(deployed), can_cplx=complex I, ref=reference-definition "
              "orientation integral. Ratio = isolated HIGH-kernel error."),
    )
    print(f"\nsaved -> {OUT / 'HIGH_fold_compare.npz'}")


if __name__ == "__main__":
    main()
