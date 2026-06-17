#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Isolate WHY deployed canoes HIGH zeta_D is 4-17x large though the kernel is OK.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Established so far:
 * probe_HIGH_fold_compare: canoes orientation kernel folded over the WHOLE plane
   with the REFERENCE B_kappa matches the reference Gamma to <1% (complex) -> the
   kernel _kappa3_limber_alpha_kernel is CORRECT.
 * probe_HIGH_deployed_vs_complex: the deployed pipeline (canoes B_delta, radial
   machinery, ell window (60,1000]) is 4-17x large for D-channels, and complex-vs-
   real makes no difference -> the bug is NOT the kernel, it is the radial/B/window
   assembly (or the ell window itself).

This probe pins WHICH by computing Gamma^3 (spins (2,2,-2)) FOUR ways, all with the
canoes orientation kernel and the canoes (u,v,phi) Gauss grid, varying ONLY the
B and the radial treatment:

  (1) ref-B, whole plane (1..2e4)              -> matches reference (baseline truth)
  (2) ref-B, ell window (60,1000]              -> isolates the WINDOW
  (3) canoes-B_delta single-plane proxy, window-> isolates the canoes B/normalisation
  (4) deployed canoes radial (LOS integral)    -> the shipped number

ref-B here is the SAME single-plane Limber convergence B_kappa the reference uses;
canoes-B is the deployed B_delta * growth^4 * response * radial measure summed over
the LOS lambda shells then folded with K^3 (the actual deployed radial integral).

If (1)~(2) (window harmless) but (2) != (3)/(4), the canoes B/normalisation is the
lever; if (1) != (2), the window is the lever.
"""
from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import BkappaLimber, BkappaInterp, ARCMIN2RAD  # noqa: E402
sys.path.insert(0, str(Path("/Users/zzhang/projects/canoes")))
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    _kappa3_limber_alpha_kernel, _kappa3_flat_triangle_vertices,
    _validate_limber_high_quadrature,
)

_DMOD = (2, 2, -2)


def sas_thetas(gam, phid):
    g = gam * ARCMIN2RAD; ph = math.radians(phid)
    t3 = 2 * g * math.sin(ph / 2)
    return t3, g, g


def refB_fold(bk, gam, phid, *, window):
    """Gamma^3 = INT (u du)(v dv) dphi/(2pi)^4 B_kappa * canoes_kernel(2,2,-2).
    window=None -> whole plane; window=(lo,hi) -> only max_leg in (lo,hi]."""
    t12, t23, t31 = sas_thetas(gam, phid)
    r2, r3 = _kappa3_flat_triangle_vertices(t12, t23, t31)
    n_l, n_phi = 220, 128
    lr = np.geomspace(1.0, 2.0e4, n_l)
    dlnl = np.gradient(np.log(lr)); wl = lr * lr * dlnl
    pa = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False); dpa = 2.0 * np.pi / n_phi
    norm = 1.0 / (2.0 * np.pi) ** 4
    U = lr[:, None, None]; V = lr[None, :, None]; PHI = pa[None, None, :]
    W = np.sqrt(np.maximum(U * U + V * V + 2 * U * V * np.cos(PHI), 0.0))
    meas = (wl[:, None, None] * wl[None, :, None] * dpa) * norm
    if window is not None:
        lo, hi = window
        ml = np.maximum(np.maximum(U, V), W)
        meas = meas * ((ml > lo) & (ml <= hi))
    B = bk(W, np.broadcast_to(U, W.shape), np.broadcast_to(V, W.shape))
    Kr = _kappa3_limber_alpha_kernel(spins=_DMOD, u=U, v=V, phi=PHI, w=W, r2=r2, r3=r3)
    return np.sum(meas * B * Kr)


def main():
    bk_exact = BkappaLimber(n_chi=120)
    bk = BkappaInterp(bk_exact, n_l=96, n_c=64)
    ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
    configs = list(zip(ref["gamma_arcmin"].tolist(), ref["phi_deg"].tolist()))

    print("=" * 100)
    print("WINDOW / B isolation for Gamma^3 (canoes kernel, ref-B)")
    print("  (1) ref-B whole plane   (2) ref-B window (60,1000]   |Gref3|=reference Gamma^3")
    print("=" * 100)
    print(f"{'config':>14s} {'|G3 whole|':>12s} {'|G3 window|':>12s} "
          f"{'|Gref3|':>12s} {'whole/ref':>9s} {'win/ref':>9s} {'win/whole':>9s}")
    rows = []
    for (gam, phid) in configs:
        idx = configs.index((gam, phid))
        gref3 = abs(ref["Gamma3"][idx])
        g_whole = refB_fold(bk, gam, phid, window=None)
        g_win = refB_fold(bk, gam, phid, window=(60.0, 1000.0))
        print(f"{f'{gam:.0f}/{phid:.0f}':>14s} {abs(g_whole):12.4e} {abs(g_win):12.4e} "
              f"{gref3:12.4e} {abs(g_whole)/gref3:9.3f} {abs(g_win)/gref3:9.3f} "
              f"{abs(g_win)/abs(g_whole):9.3f}")
        rows.append((gam, phid, complex(g_whole), complex(g_win), gref3))

    print("\nINTERPRETATION:")
    print("  whole/ref ~ 1  => kernel+ref-B correct on the full plane.")
    print("  win/ref far from whole/ref => the (60,1000] WINDOW alone shifts the")
    print("     magnitude (the spin-2 J_2 weight concentrates differently than the")
    print("     full-plane integral; the deployed HIGH tail is NOT the full Gamma).")
    np.savez(OUT / "HIGH_window_B_isolation.npz",
             gam=np.array([r[0] for r in rows]), phi=np.array([r[1] for r in rows]),
             G3_whole=np.array([r[2] for r in rows], dtype=np.complex128),
             G3_window=np.array([r[3] for r in rows], dtype=np.complex128),
             Gref3=np.array([r[4] for r in rows]),
             note="canoes kernel + ref-B, whole plane vs (60,1000] window, vs reference Gamma^3.")
    print(f"\nsaved -> {OUT / 'HIGH_window_B_isolation.npz'}")


if __name__ == "__main__":
    main()
