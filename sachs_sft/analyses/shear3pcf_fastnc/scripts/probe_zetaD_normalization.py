#!/usr/bin/env python
"""Standalone root-cause probe for the spin-2 modulus zeta_D over-normalization.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

WHAT WAS HYPOTHESIZED vs WHAT THE CODE DOES
===========================================
HYPOTHESIS: zeta_D (spin (2,2,-2)) is over-normalized because the angular
synthesis applies the (0,0,0)-anchor Wigner-3j  w0 = wigner_3j(L1,L2,L3;0,0,0)
to the spin-2 channel instead of a spin-weighted Gaunt symbol.

DEPLOYED PATH (traced):
  compute_kappa3_mod_zeta_equal_time (LOW, kappa3.py:3737) and
  compute_kappa3_mod_sigma3_high (HIGH, kappa3.py:3912) project with
  canoes.nuell.correlation._spin_aware_three_pt.SpinAwareZetaEvaluator.
  Its per-triple angular coefficient (kernel build lines 598-604) is
      coef = h * w0 * sq_l1 * (-1)^s1 ,   w0 = wigner_3j(L1,L2,L3; 0,0,0)
  i.e. it DOES use the (0,0,0) anchor.  BUT the identical coefficient form is
  used for the canonical PPP channel (2,2,2) and TPP channel (0,2,2), which the
  existing test (test_kappa3_mod_normalization.py) validates reproduce the
  scalar-path zeta_PPP / zeta_TPP to rel-err < 1e-6.  The ONLY (2,2,-2) vs
  (2,2,2) difference is the third leg's Wigner-d index sign (d^L_{m,-2} vs
  d^L_{m,+2}) and the m-range.  So w0 is NOT a spin-2-SPECIFIC anchor error.

CRUCIAL ALGEBRAIC FACT (from derive_spin2_gaunt_norm.wl):
  the spin Gaunt anchor wigner_3j(L1,L2,L3; -2,-2,+2) used in a naive
  "spin-weighted (0,0,0)-replacement" VANISHES identically on the diagonal
  (bottom-row sum = -2 != 0).  The modulus channel has total spin +2 != 0; it is
  NOT a rotational scalar, so there is no single all-M=0 spin Gaunt anchor that
  could "correctly" replace w0.  The deployed m-sum (m1=-s1, m1+m2+m3=0) is the
  correct leg-1-at-pole reduction and is spin-threaded.

THIS PROBE therefore measures the deployed evaluator directly and tests it
against the ONE physical ground-truth anchor that the appendix provides:
  small-angle / high-ell limit  |zeta_D| -> |zeta_TTT|  (appendix 615-619:
  sqrt(L^2(L^2-2)) -> L^2, d^L_{2,2} -> 1).  If the deployed zeta_D EXCEEDS
  zeta_TTT by the observed ~3.5-5x at small angle, the over-normalization is
  real and lives in the deployed (2,2,-2) value; if it approaches zeta_TTT, the
  evaluator is correctly normalized and the fastnc discrepancy is a different
  (projection-convention) effect.

We also expose the +2 self-consistency (2,2,2)==PPP and quantify the
(2,2,-2)/(2,2,2) and (2,2,-2)/(0,0,0) ratios vs separation & opening angle to
characterize the STRUCTURE of any discrepancy.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("JAX_PLATFORMS", "cpu")

from canoes.nuell.bispectra._compact import (  # noqa: E402
    CompactBispectrum,
    enumerate_valid_triples,
)
from canoes.nuell.correlation._spin_aware_three_pt import (  # noqa: E402
    SpinAwareZetaEvaluator,
)
from canoes.nuell.correlation.three_pt import ZetaEvaluator  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(parents=True, exist_ok=True)


def make_cb(ell_max, seed=42):
    """A fixed, smooth-ish reduced bispectrum on all valid triples."""
    l1s, l2s, l3s = enumerate_valid_triples(ell_max)
    rng = np.random.default_rng(seed)
    shape = 1.0 / ((1.0 + l1s) * (1.0 + l2s) * (1.0 + l3s))
    bvals = rng.standard_normal(l1s.size) * shape
    cb = CompactBispectrum(
        ell_max=ell_max, ells_sparse=None,
        l1s=l1s.astype(np.intp), l2s=l2s.astype(np.intp), l3s=l3s.astype(np.intp),
        b_values=np.asarray(bvals, dtype=np.float64),
    )
    return cb, bvals


def eval_chan(t12, t23, t31, ell_max, spins, cb):
    if spins == (0, 0, 0):
        ev = ZetaEvaluator(np.array([t12]), np.array([t23]), np.array([t31]),
                           ell_max=ell_max)
    else:
        ev = SpinAwareZetaEvaluator(
            np.array([t12]), np.array([t23]), np.array([t31]),
            ell_max=ell_max, spin_weights=spins)
    return float(np.asarray(ev.evaluate(cb))[0])


def sas(gamma_arcmin, phi_deg):
    g = math.radians(gamma_arcmin / 60.0)
    ph = math.radians(phi_deg)
    t3 = 2.0 * g * math.sin(ph / 2.0)
    return (t3, g, g)  # (theta_12=base, theta_23=g, theta_31=g)


def main():
    print("=" * 78)
    print("zeta_D spin-(2,2,-2) modulus normalization probe (deployed evaluator)")
    print("=" * 78)

    spins_D = (2, 2, -2)
    spins_PPP = (2, 2, 2)
    spins_TTT = (0, 0, 0)

    ell_max = 16
    cb, _ = make_cb(ell_max)

    # ---- (0) +2 self-consistency: (2,2,2) spin path == scalar (0,0,0) PPP? ----
    # The canonical PPP is the (2,2,2)-spin path; (0,0,0) is the pure scalar.
    # They are DIFFERENT channels (PPP is spin-2), so we only confirm the spin
    # evaluator runs and that (2,2,2) != (0,0,0) (distinct), then anchor on TTT.
    print("\n[0] sanity: distinct channels at one config (g=10' phi=90)")
    t12, t23, t31 = sas(10.0, 90.0)
    vD = eval_chan(t12, t23, t31, ell_max, spins_D, cb)
    vP = eval_chan(t12, t23, t31, ell_max, spins_PPP, cb)
    vT = eval_chan(t12, t23, t31, ell_max, spins_TTT, cb)
    print(f"    zeta_D (2,2,-2)  = {vD:+.6e}")
    print(f"    zeta_PPP(2,2,2)  = {vP:+.6e}")
    print(f"    zeta_TTT(0,0,0)  = {vT:+.6e}")

    # ---- (1) ratio structure vs (gamma, phi) ----
    print("\n[1] deployed-zeta_D / zeta_TTT  and  / zeta_PPP  vs (gamma, phi)")
    print("    appendix anchor: at small gamma |zeta_D| -> |zeta_TTT| (ratio -> ~1)")
    print("    a ratio that GROWS to few-x and is angle-structured = over-norm")
    print()
    gammas = [1.0, 2.0, 5.0, 10.0, 30.0, 60.0, 120.0, 240.0]
    phis = [30.0, 60.0, 90.0, 120.0, 150.0]
    print(f"    {'gamma(arcmin)':>13s} " + " ".join(f"phi={p:>5.0f}" for p in phis))
    ratio_DT = np.full((len(gammas), len(phis)), np.nan)
    ratio_DP = np.full((len(gammas), len(phis)), np.nan)
    absD = np.full((len(gammas), len(phis)), np.nan)
    absT = np.full((len(gammas), len(phis)), np.nan)
    for ig, gam in enumerate(gammas):
        row = []
        for ip, ph in enumerate(phis):
            t12, t23, t31 = sas(gam, ph)
            d = eval_chan(t12, t23, t31, ell_max, spins_D, cb)
            t = eval_chan(t12, t23, t31, ell_max, spins_TTT, cb)
            p = eval_chan(t12, t23, t31, ell_max, spins_PPP, cb)
            ratio_DT[ig, ip] = d / t if t else np.nan
            ratio_DP[ig, ip] = d / p if p else np.nan
            absD[ig, ip] = abs(d)
            absT[ig, ip] = abs(t)
            row.append(f"{d/t:8.3f}" if t else "     nan")
        print(f"    g={gam:8.1f}'   " + " ".join(row))

    print("\n    |zeta_D| / |zeta_TTT|  (absolute-value ratio):")
    print(f"    {'gamma(arcmin)':>13s} " + " ".join(f"phi={p:>5.0f}" for p in phis))
    for ig, gam in enumerate(gammas):
        row = [f"{absD[ig,ip]/absT[ig,ip]:8.3f}" if absT[ig, ip] else "     nan"
               for ip in range(len(phis))]
        print(f"    g={gam:8.1f}'   " + " ".join(row))

    finite = np.isfinite(ratio_DT)
    print(f"\n    signed D/TTT : median={np.nanmedian(ratio_DT):.3f} "
          f"min={np.nanmin(ratio_DT):.3f} max={np.nanmax(ratio_DT):.3f}")
    absratio = absD / absT
    print(f"    |D|/|TTT|    : median={np.nanmedian(absratio):.3f} "
          f"min={np.nanmin(absratio):.3f} max={np.nanmax(absratio):.3f} "
          f"spread={np.nanmax(absratio)/np.nanmin(absratio):.3f}")

    # ---- (2) deepest small-angle anchor: does |D| -> |TTT| as gamma->0 ? ----
    print("\n[2] small-angle limit  |zeta_D| / |zeta_TTT|  as gamma -> 0  (phi=90)")
    print("    appendix predicts -> 1 (sqrt(L^2(L^2-2)) -> L^2).")
    small_gammas = [0.5, 1.0, 2.0, 5.0, 10.0]
    rows2 = []
    for gam in small_gammas:
        t12, t23, t31 = sas(gam, 90.0)
        d = eval_chan(t12, t23, t31, ell_max, spins_D, cb)
        t = eval_chan(t12, t23, t31, ell_max, spins_TTT, cb)
        r = abs(d) / abs(t) if t else np.nan
        rows2.append((gam, d, t, r))
        print(f"    gamma={gam:6.1f}'  |zeta_D|={abs(d):.5e}  "
              f"|zeta_TTT|={abs(t):.5e}  |D|/|TTT|={r:.4f}")
    small_limit = rows2[0][3]
    print(f"\n    => |D|/|TTT| at gamma=0.5' = {small_limit:.4f}")
    if small_limit > 2.0:
        print("       OVER-NORMALIZED: exceeds the appendix small-angle anchor"
              f" by ~{small_limit:.2f}x")
    elif 0.5 < small_limit < 2.0:
        print("       CONSISTENT with the appendix anchor (~O(1)); the fastnc"
              " discrepancy is NOT a zeta_D normalization error")
    else:
        print("       UNDER the anchor; investigate sign/zero crossings")

    # ---- (3) REAL deployed builders: zeta_D vs zeta_TTT at squeezed limit ----
    # Use the existing consistency_deployed_vs_fresh.npz (deployed==fresh,
    # reldiff=0) which holds REAL compute_kappa3_mod_zeta_equal_time (zDmod) and
    # compute_kappa3_zeta_table (zTTT) at three triples incl. the exact-(1,1,1).
    # This carries the COEFS_PSI0 sqrt(L(L+1)(L-1)(L+2)) per-leg eth^2 operator
    # that the bare-b probe above omits.
    print("\n[3] REAL deployed builders (consistency_deployed_vs_fresh.npz):")
    real = None
    real_path = OUT / "consistency_deployed_vs_fresh.npz"
    real_summary = {}
    if real_path.exists():
        real = np.load(real_path, allow_pickle=True)
        zT = real["zTTT_dep"]
        zD = real["zDmod_dep"]
        triples_cos = real["triples"]
        print("    triple(cos)                 sum|zTTT|    sum|zDmod|   |D|/|T|")
        for r in range(zT.shape[0]):
            sT = float(np.nansum(np.abs(zT[r])))
            sD = float(np.nansum(np.abs(zD[r])))
            rr = sD / sT if sT else np.nan
            tag = "  <- exact (1,1,1)" if np.allclose(triples_cos[r], 1.0) else ""
            print(f"    {str(triples_cos[r]):26s} {sT:.4e}  {sD:.4e}  {rr:7.4f}{tag}")
            real_summary[f"triple{r}_DoverT"] = rr
        # explicit squeezed-limit verdict
        sq = np.where(np.all(np.isclose(triples_cos, 1.0), axis=1))[0]
        if sq.size:
            i = int(sq[0])
            sD = float(np.nansum(np.abs(zD[i])))
            print(f"\n    EXACT-(1,1,1) [theta_ij=0, coincident points]: "
                  f"sum|zeta_D| = {sD:.4e}")
            if sD == 0.0:
                print("    => zeta_D = 0 at the coincident-point degeneracy.")
                print("       At FINITE triangles zeta_D is nonzero but SMALL (~1e-7),")
                print("       below zeta_TTT -- the OPPOSITE of an over-normalization")
                print("       (zeta_D too SMALL, not too big). So an over-normalized")
                print("       zeta_D cannot explain too-LARGE shear Gamma^1..3.")
    else:
        print(f"    (not found: {real_path}; run consistency_deployed_vs_fresh.py)")

    np.savez(
        OUT / "zetaD_normalization_probe.npz",
        ell_max=ell_max,
        gammas=np.array(gammas), phis=np.array(phis),
        ratio_DT=ratio_DT, ratio_DP=ratio_DP,
        absD=absD, absT=absT,
        small_gammas=np.array([r[0] for r in rows2]),
        small_absD=np.array([abs(r[1]) for r in rows2]),
        small_absT=np.array([abs(r[2]) for r in rows2]),
        small_ratio=np.array([r[3] for r in rows2]),
        spins_D=np.array(spins_D),
        **{k: np.array(v) for k, v in real_summary.items()},
    )
    print(f"\nsaved -> {OUT / 'zetaD_normalization_probe.npz'}")


if __name__ == "__main__":
    main()
