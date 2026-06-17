#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Root-cause diagnostic for the +0.44 ell-slope of W_can / B_kappa.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

GOAL (diagnostic B, task 1+3):
  (1) Re-measure the canoes-HIGH B_delta / B_kappa-ref log-log slope in ell, both
      from the saved decomp and by re-evaluating the canoes HIGH B_delta kernel
      (compute_kappa3_sigma3_high path, tree SPT F2 P P) directly across ell.
  (3) Decide whether +0.44 is a MISSING/EXTRA radial factor (power of ell, k, 1+z,
      F2 arg bug, window/Limber mismatch) OR a genuine tree-Limber-vs-fastnc
      modeling difference.

METHOD: the two LOS integrands both fold the SAME b_delta(l/chi) (same SPT F2,
same P(k), same k=l/chi).  The ratio of the chi-integrals is therefore governed
ENTIRELY by the chi-WEIGHTING W(chi) each route applies to b_delta(l/chi):

    Gamma_can(l) ~ INT dchi  W_can(chi) b_delta(l/chi)
    Gamma_ref(l) ~ INT dchi  W_ref(chi) b_delta(l/chi)

If W_can(chi) propto W_ref(chi) (same chi-shape), the ratio is l-INDEPENDENT (flat
slope, no +0.44).  The +0.44 means W_can(chi) has a DIFFERENT chi-shape -- it
weights high-chi (low-k for fixed l) shells more, so as l grows the integral
samples b_delta at a different effective chi and the ratio drifts.

We (a) extract W_can(chi), W_ref(chi); (b) show their SHAPE ratio vs chi;
(c) attribute the chi-shape mismatch factor-by-factor; (d) verify that replacing
W_can's chi-shape with W_ref's chi-shape KILLS the +0.44 slope (boundary check).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"
sys.path.insert(0, str(ANALYSIS))
from stage1_bkappa_phase_reference import BkappaLimber, H0_HUNIT, OMEGA_M  # noqa: E402

sys.path.insert(0, "/Users/zzhang/projects/canoes")
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    chi_of_lambda, build_lightcone, _f2_kernel_numpy,
    _kappa3_limber_response_product_h, _lambda_jacobian,
    _kappa3_radial_density_conversion, _H0_H_PER_MPC,
)
from scipy.integrate import cumulative_trapezoid  # noqa: E402


def build_weights(cosmo, n_chi=120):
    """Return (chi_h, dchi-fold weights) for canoes-HIGH and reference, on a
    COMMON chi grid, both as pure chi-density LOS weights W(chi) so that
    Gamma(l) = INT dchi W(chi) b_delta(l/chi)."""
    # background
    zf = np.linspace(0, 5, 8000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    chi_s = float(chif[-1])
    chi_s_h = chi_s * cosmo.h

    # common z grid (canoes-style equal-z shells)
    z = np.linspace(5.0 / n_chi, 5.0, n_chi)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * cosmo.h
    a = 1.0 / (1.0 + z)
    opz = 1.0 + z
    chi_h = chi_phys * cosmo.h  # Mpc/h  (k = l/chi_h, h/Mpc)

    # canoes lightcone (growth + opz on canoes' projected chi grid)
    chi_arr = np.maximum(chi_of_lambda(cosmo, lam_h, convention="project"), 1e-6)
    lc, opz_c = build_lightcone(cosmo, chi_arr)
    growth = np.asarray(lc.M[0] / opz_c, dtype=np.float64)  # D(z), D(0)=1
    poisson_amp = -1.5 * cosmo.Omega_m * (_H0_H_PER_MPC ** 2) * opz_c
    djdc = _lambda_jacobian(opz_c, convention="project")
    rf, _, _ = _kappa3_radial_density_conversion(djdc, "lambda")  # (1+z)^-4
    resp_PPP = np.array([
        _kappa3_limber_response_product_h("PPP", float(poisson_amp[i]),
                                          float(opz_c[i]))
        for i in range(len(opz_c))
    ])
    h6 = cosmo.h ** 6
    # convergence lensing kernel (physical), as in stage1_ours_v2
    Kg = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s

    # ---- canoes-HIGH chi-density LOS weight (b_delta factored out) -------------
    # The deployed fold is Gamma ~ INT dlambda Kg^3 * zeta,
    #   zeta = growth^4 * b_delta / chi_h^4 * resp_PPP * rf * h6.
    # INT dlambda = INT a^2 dchi_phys.  We want a pure dchi_phys density times
    # b_delta, so multiply by (dlambda/dchi_phys)=a^2 and absorb growth^4 (b_delta
    # already carries D^4? NO -- canoes b_delta here is the z=0 P(k) SPT bispectrum
    # WITHOUT growth; growth^4 is the explicit factor).  So:
    W_can = (Kg ** 3) * (a ** 2) * (growth ** 4) / (chi_arr ** 4) * resp_PPP * rf * h6

    # ---- reference B_kappa chi-density LOS weight (b_delta factored out) -------
    # B_kappa = INT dchi g^3 (1/chi) (1+z)^3 D^4 b_delta,  g=1.5 H0^2 Om(1-chi/chis)
    g_ref = 1.5 * (H0_HUNIT ** 2) * OMEGA_M * (1.0 - chi_h / chi_s_h)
    W_ref = (g_ref ** 3) * (1.0 / chi_h) * (opz ** 3) * (growth ** 4)

    return dict(z=z, chi_phys=chi_phys, chi_h=chi_h, chi_arr=chi_arr,
                opz=opz, a=a, growth=growth, Kg=Kg, resp_PPP=resp_PPP, rf=rf,
                W_can=W_can, W_ref=W_ref, g_ref=g_ref, h6=h6, chi_s_h=chi_s_h)


def gamma_of_l(W, chi_phys_for_k, pk, f2, Ls):
    """Gamma(l) = INT dchi W(chi) b_delta(l/chi) at equilateral (cos=-0.5).
    chi_phys_for_k must be in the SAME units as k=l/chi expects (h/Mpc => chi_h)."""
    out = np.zeros(len(Ls))
    for j, L in enumerate(Ls):
        k = L / chi_phys_for_k
        p = pk(k)
        f = f2(k, k, np.full_like(k, -0.5))
        bdelta = 2.0 * (3.0 * f * p * p)  # 3 identical cyc terms, equilateral
        out[j] = np.trapezoid(W * bdelta, chi_phys_for_k)
    return out


def slope(x, y):
    return np.polyfit(np.log(x), np.log(np.abs(y)), 1)[0]


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    f2 = lambda ka, kb, c: _f2_kernel_numpy(np.asarray(ka), np.asarray(kb),
                                            np.asarray(c))
    W = build_weights(cosmo, n_chi=160)
    chi_h = W["chi_h"]
    Ls = np.array([60., 80., 100., 150., 200., 300., 500., 800., 1000.,
                   2000., 5000.])

    # both routes integrate over chi_h with k=l/chi_h
    G_can = gamma_of_l(W["W_can"], chi_h, pk, f2, Ls)
    G_ref = gamma_of_l(W["W_ref"], chi_h, pk, f2, Ls)
    R = G_can / G_ref
    s_full = slope(Ls, R)
    s_lo = slope(Ls[Ls <= 200], R[Ls <= 200])
    s_hi = slope(Ls[(Ls >= 200) & (Ls <= 1000)], R[(Ls >= 200) & (Ls <= 1000)])
    s_doc = slope(Ls[Ls <= 1000], R[Ls <= 1000])

    print("=" * 90)
    print("DIRECT canoes-HIGH B_delta / B_kappa-ref ratio vs ell (equilateral, "
          "same b_delta both sides)")
    print("=" * 90)
    print(f"{'L':>7s} {'G_can':>13s} {'G_ref':>13s} {'R=can/ref':>12s}")
    for j, L in enumerate(Ls):
        print(f"{L:7.0f} {G_can[j]:13.4e} {G_ref[j]:13.4e} {R[j]:12.4e}")
    print(f"\n  slope(R) full [60..5000]   = {s_full:+.4f}")
    print(f"  slope(R) [60..1000] (==doc)= {s_doc:+.4f}   (decomp doc = +0.4395)")
    print(f"  slope(R) low [60..200]     = {s_lo:+.4f}")
    print(f"  slope(R) high [200..1000]  = {s_hi:+.4f}")

    # ---- chi-shape mismatch: normalize each W to unit peak, show ratio vs chi ---
    Wc = W["W_can"] / np.nanmax(W["W_can"])
    Wr = W["W_ref"] / np.nanmax(W["W_ref"])
    shape_ratio = (W["W_can"] / W["W_ref"])
    shape_ratio = shape_ratio / shape_ratio[0]  # normalize at near-side
    print("\n" + "=" * 90)
    print("chi-SHAPE of the two LOS weights (normalized at near-side chi). If "
          "W_can/W_ref were flat,\nthe ratio would be ell-independent (no +0.44).")
    print("=" * 90)
    print(f"{'z':>5s} {'chi_h':>9s} {'W_can(norm)':>12s} {'W_ref(norm)':>12s} "
          f"{'(Wc/Wr)/near':>13s}")
    for i in range(0, 160, 16):
        print(f"{W['z'][i]:5.2f} {chi_h[i]:9.1f} {Wc[i]:12.4e} {Wr[i]:12.4e} "
              f"{shape_ratio[i]:13.4e}")

    # ---- factor-by-factor attribution of the W_can/W_ref chi-shape -------------
    # W_can/W_ref = [Kg^3 a^2 growth^4/chi_arr^4 resp rf h6] / [g^3 (1/chi_h)(1+z)^3 growth^4]
    # growth^4 cancels. Remaining chi-shape pieces:
    a = W["a"]; opz = W["opz"]; Kg = W["Kg"]; resp = W["resp_PPP"]; rf = W["rf"]
    g_ref = W["g_ref"]; chi_arr = W["chi_arr"]
    print("\n" + "=" * 90)
    print("Factor-by-factor chi-slope d ln(factor)/d ln(chi_h)  (drives the shape "
          "mismatch)")
    print("=" * 90)
    factors = {
        "Kg^3 (canoes lens kernel, phys)": Kg ** 3,
        "a^2 (dlambda/dchi)": a ** 2,
        "1/chi_arr^4 (canoes)": 1.0 / chi_arr ** 4,
        "resp_PPP (=A^3 (1+z)^15)": resp,
        "rf=(1+z)^-4": rf,
        "--- canoes W_can total ---": W["W_can"],
        "g_ref^3 (ref lens kernel)": g_ref ** 3,
        "1/chi_h (ref)": 1.0 / chi_h,
        "(1+z)^3 (ref)": opz ** 3,
        "--- ref W_ref total ---": W["W_ref"],
        ">>> W_can / W_ref (shape) <<<": W["W_can"] / W["W_ref"],
    }
    lnchi = np.log(chi_h)
    for name, fac in factors.items():
        fac = np.asarray(fac, dtype=float)
        m = np.isfinite(fac) & (fac > 0)
        s = np.polyfit(lnchi[m], np.log(fac[m]), 1)[0]
        print(f"  d ln/d ln chi = {s:+8.3f}   {name}")

    # ---- BOUNDARY CHECK: swap canoes chi-shape -> ref chi-shape, slope vanishes -
    # Build a hybrid weight: ref chi-shape but canoes overall normalization.
    W_hybrid = W["W_ref"] * (np.nansum(W["W_can"]) / np.nansum(W["W_ref"]))
    G_hyb = gamma_of_l(W_hybrid, chi_h, pk, f2, Ls)
    R_hyb = G_hyb / G_ref
    s_hyb = slope(Ls[Ls <= 1000], R_hyb[Ls <= 1000])
    print("\n" + "=" * 90)
    print("BOUNDARY CHECK: give the canoes route the REFERENCE chi-shape "
          "(rescaled to same total).")
    print("If the +0.44 is purely the chi-shape mismatch, this slope -> 0.")
    print("=" * 90)
    print(f"  slope(R_hybrid) [60..1000] = {s_hyb:+.4f}  (target ~ 0.000)")

    np.savez(
        OUT / "ellslope_root_diagnostic.npz",
        Ls=Ls, G_can=G_can, G_ref=G_ref, R=R, slope_full=s_full,
        slope_doc=s_doc, slope_lo=s_lo, slope_hi=s_hi,
        chi_h=chi_h, z=W["z"], W_can=W["W_can"], W_ref=W["W_ref"],
        shape_ratio=shape_ratio, slope_hybrid=s_hyb, R_hybrid=R_hyb,
        note="canoes-HIGH/B_kappa-ref ell-slope root diagnostic; chi-shape attribution.",
    )
    print(f"\nsaved -> {OUT / 'ellslope_root_diagnostic.npz'}")


if __name__ == "__main__":
    main()
