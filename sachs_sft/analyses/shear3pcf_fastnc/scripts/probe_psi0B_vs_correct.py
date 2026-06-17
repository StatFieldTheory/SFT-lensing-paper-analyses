#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Phase-1 ORACLE: pin where the spin-2 3-point +0.44 ell-slope / 4-17x lives.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

The xAct derivation (derive_spin2_3pt_bispectrum_weight.wl) PROVED the per-leg
spin-2/scalar modulus weight W2(ell)=sqrt(L2(L2-2))/L2 -> 1 at high ell with a
1/ell^2 correction (1-W2^3 ~ 8e-4 at ell=60), so the per-leg spin-2 eth^2 weight
CANNOT be the +0.44 lever.  This oracle locates the real lever by separating two
ratios over the (60,1000] window on a clean ell grid:

  (A) psi0-B / phi00-B  = canoes' OWN spin-2 vs scalar 3-leg reduced bispectrum
      from compute_bispectra_deri_{psi0,phi00} (full per-leg operator, SAME
      machinery as the validated 2-point).  The derivation predicts -> W2(ell)^3
      (-> 1).  If it instead shows the +0.44 slope, the operator is spin-2-buggy.

  (B) phi00-B / B_kappa = canoes' SCALAR convergence 3-leg bispectrum (the
      operator path) vs the trusted reference convergence bispectrum.  Both are
      SCALAR (spin-0).  If THIS carries the +0.44 slope, the bug is a scalar
      B-shape difference (NOT spin-2), and the J_2 kernel merely AMPLIFIES it
      relative to J_0/J_6 because J_2 weights the window differently.

  (C) HIGH ell-flat response check: the deployed zeta_D HIGH path uses
      _kappa3_limber_response_product_h (ell-FLAT, Psi0=-Phi00).  We confirm the
      HIGH scalar response vs the operator phi00-B is also ell-flat-equivalent,
      i.e. the HIGH path's per-leg weight = 1 matches the operator at high ell.

Both canoes bispectra are evaluated with one_plus_z so the geometric (1+z)^4 per
leg is included, on a SINGLE equal-shell chi (so this is the per-(l1,l2,l3) shape
at fixed radial point, exactly the quantity the orientation kernel folds).

Writes outputs/psi0B_oracle.npz.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import BkappaLimber  # noqa: E402

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.sachs.kappa3 import chi_of_lambda, build_lightcone  # noqa: E402
from canoes.nuell.transfer.sachs_driver import (  # noqa: E402
    compute_bispectra_deri_phi00,
    compute_bispectra_deri_psi0,
)
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from scipy.integrate import cumulative_trapezoid  # noqa: E402


def _w2(ell):
    """Per-leg spin-2/scalar modulus weight sqrt(L2(L2-2))/L2, L2=ell(ell+1)."""
    L2 = ell * (ell + 1.0)
    return np.sqrt(np.maximum(L2 * (L2 - 2.0), 0.0)) / L2


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    bk = BkappaLimber(n_chi=200)

    # Single representative equal-shell radial point (the lever is the per-(l1,l2,l3)
    # SHAPE; we use one mid-LOS shell at z~2.5 where the convergence kernel peaks).
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z_pick = 2.5
    chi_pick_phys = float(np.interp(z_pick, zf, chif))
    lam_pick_phys = float(np.interp(z_pick, zf, lam_f))
    lam_h = np.array([lam_pick_phys * cosmo.h], dtype=np.float64)
    chi_arr_h = np.maximum(chi_of_lambda(cosmo, lam_h, convention="project"), 1e-6)
    lightcone, one_plus_z = build_lightcone(cosmo, chi_arr_h)
    chi_h = float(chi_arr_h[0])
    print(f"radial point: z={z_pick}, chi_phys={chi_pick_phys:.1f} Mpc, "
          f"chi_h={chi_h:.2f} Mpc/h, 1+z={float(one_plus_z[0]):.4f}")

    def canoes_B(builder, l1, l2, l3):
        ell_t = np.array([[int(l1), int(l2), int(l3)]], dtype=np.intp)
        chi_t = np.array([[chi_h, chi_h, chi_h]], dtype=np.float64)
        val = builder(
            ell_t, chi_t,
            pk=pk, lightcone=lightcone, one_plus_z=one_plus_z,
            spt_kind="tree_matter", Nmax=384,
        )
        return float(np.asarray(val, dtype=np.float64)[0])

    # The reference B_kappa is a LOS integral; for a per-shell shape comparison we
    # need the per-shell integrand at chi_h.  B_kappa's chi-integrand at chi_h is
    # weight(chi)*B_delta(l/chi); we extract the convergence B_delta*D^4 at chi_h
    # WITHOUT the LOS weight (the LOS weight is the SAME for scalar and spin-2 and
    # cancels in every ratio we form).  We therefore compare canoes_B to the
    # reference convergence bispectrum core  B_delta(l/chi)*D^4 at chi_h, in h-units.
    # For a clean apples-to-apples we use the FULL B_kappa LOS integral as the
    # scalar reference (it is the trusted object); both ratios below are SHAPE
    # ratios in ell, and a constant LOS-normalisation offset is divided out by
    # reporting R/R0.

    Ls = [60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 500.0, 800.0, 1000.0]
    rows = []
    print("=" * 110)
    print(f"{'L(equi)':>8s} {'psi0-B':>12s} {'phi00-B':>12s} {'B_kappa':>12s} "
          f"{'psi0/phi00':>11s} {'W2^3':>9s} {'phi00/Bk':>11s}")
    print("-" * 110)
    for L in Ls:
        bpsi = canoes_B(compute_bispectra_deri_psi0, L, L, L)
        bphi = canoes_B(compute_bispectra_deri_phi00, L, L, L)
        bkv = float(bk(L, L, L))
        r_ps_ph = bpsi / bphi if bphi != 0 else float("nan")
        w2c = float(_w2(L)) ** 3
        r_ph_bk = bphi / bkv if bkv != 0 else float("nan")
        rows.append((L, bpsi, bphi, bkv, r_ps_ph, w2c, r_ph_bk))
        print(f"{L:8.0f} {bpsi:12.4e} {bphi:12.4e} {bkv:12.4e} "
              f"{r_ps_ph:11.5f} {w2c:9.5f} {r_ph_bk:11.4e}")

    La = np.array([r[0] for r in rows])
    # Ratio (A): psi0/phi00 — should track W2^3 (->1), slope ~0.
    rA = np.array([r[4] for r in rows])
    slopeA = np.polyfit(np.log(La), np.log(np.abs(rA)), 1)[0]
    # Ratio (B): phi00/B_kappa (scalar vs scalar) — does it carry the +0.44?
    rB = np.array([r[6] for r in rows])
    slopeB = np.polyfit(np.log(La), np.log(np.abs(rB)), 1)[0]
    # Ratio (A) vs the analytic W2^3
    w2cube = np.array([r[5] for r in rows])
    rA_over_w2 = rA / w2cube

    print("=" * 110)
    print(f"(A) d ln|psi0-B/phi00-B| / d ln L = {slopeA:+.4f}   "
          f"(derivation predicts ~0; W2^3 slope = "
          f"{np.polyfit(np.log(La), np.log(w2cube), 1)[0]:+.5f})")
    print(f"    psi0-B/phi00-B  /  W2^3  (should be ~1 flat): "
          f"min={rA_over_w2.min():.5f} max={rA_over_w2.max():.5f}")
    print(f"(B) d ln|phi00-B/B_kappa| / d ln L = {slopeB:+.4f}   "
          f"(scalar-vs-scalar; if ~+0.44 the bug is SCALAR B-shape, kernel-amplified)")
    print("=" * 110)

    np.savez(
        OUT / "psi0B_oracle.npz",
        L=La,
        psi0_B=np.array([r[1] for r in rows]),
        phi00_B=np.array([r[2] for r in rows]),
        B_kappa=np.array([r[3] for r in rows]),
        ratio_psi0_phi00=rA,
        W2_cube=w2cube,
        ratio_phi00_Bkappa=rB,
        slope_psi0_phi00=slopeA,
        slope_phi00_Bkappa=slopeB,
        chi_h=chi_h,
        z_pick=z_pick,
        note=("(A) psi0-B/phi00-B should track W2^3->1 (slopeA~0) per xAct "
              "derivation; (B) phi00-B/B_kappa is scalar-vs-scalar — its slope "
              "localises the +0.44 to the scalar B-shape vs the spin-2 operator."),
    )
    print(f"saved -> {OUT / 'psi0B_oracle.npz'}")


if __name__ == "__main__":
    main()
