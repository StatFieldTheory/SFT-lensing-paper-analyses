#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""Decompose W_can/B_kappa per-triangle into every radial factor to find the
ell-dependent +0.44 lever.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Both W_can and B_kappa are SCALAR convergence-3pt weights that multiply the SAME
orientation kernel.  The deployed shear fold is

    Gamma ~ INT dlambda K_g(lambda)^3 * zeta_HIGH(lambda; l1,l2,l3)
    zeta_HIGH = sum_{u,v,phi} measure/(2pi)^4 * W_lambda * kernel,
    W_lambda  = growth^4 * b_delta(l_i/chi) / chi^4 * resp_PPP * radial_factor * h^6.

The reference convergence weight is

    B_kappa = INT dchi  g(chi)^3 * (1/chi) * (1+z)^3 * D^4 * b_delta(l_i/chi),
    g(chi)  = (3/2) H0^2 Omega_m (1 - chi/chi_s)   [h-units].

This probe forms, at a FIXED equilateral triangle L, the per-shell ratio of the
two LOS integrands (collapsing the (u,v,phi) sum, which is shared), and the
fully-folded per-triangle ratio, factor by factor:

    canoes per-shell scalar weight  S_can(chi; L) =
        K_g^3 * growth^4 * b_delta * resp_PPP / chi^4 * radial_factor * h^6
    reference per-shell scalar weight S_ref(chi; L) =
        g^3 * (1/chi) * (1+z)^3 * D^4 * b_delta

    Ratio R(chi; L) = S_can / S_ref.

Because b_delta is IDENTICAL (same F2, same P(k), same k=l/chi), it cancels and
R is a PURE radial/normalisation ratio that does NOT depend on L at fixed chi...
UNLESS some factor secretly depends on l via chi-resampling.  We print R(chi;L)
for several L to expose any residual L-dependence, and the chi-integrated
W_can/B_kappa to reproduce the +0.44 slope and attribute it.
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


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    bk = BkappaLimber(n_chi=200)

    # canoes LOS lambda grid + convergence kernel (mirror stage1_ours_v2)
    zf = np.linspace(0, 5, 6000)
    chif = np.array([chi_of_z(z, Omega_m=cosmo.Omega_m, h=cosmo.h) for z in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(5.0 / 40, 5.0, 40)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * cosmo.h
    chi_s = float(chif[-1])
    a = 1.0 / (1.0 + z)
    Kg = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s          # convergence kernel (phys)
    chi_arr = np.maximum(chi_of_lambda(cosmo, lam_h, convention="project"), 1e-6)
    lc, opz = build_lightcone(cosmo, chi_arr)
    poisson_amp = -1.5 * cosmo.Omega_m * (_H0_H_PER_MPC ** 2) * opz
    djdc = _lambda_jacobian(opz, convention="project")
    rf, _, _ = _kappa3_radial_density_conversion(djdc, "lambda")
    growth = np.asarray(lc.M[0] / opz, dtype=np.float64)        # D(z), D(0)=1

    # reference per-shell scalar weight pieces (on canoes' chi grid, h-units)
    # g(chi) = (3/2) H0^2 Om (1 - chi/chi_s); weight = g^3 * (1/chi) * (1+z)^3 * D^4
    chi_h_arr = chi_arr                                          # Mpc/h
    chi_s_h = chi_s * cosmo.h
    g_ref = 1.5 * (H0_HUNIT ** 2) * OMEGA_M * (1.0 - chi_h_arr / chi_s_h)
    Dz = growth                                                 # same growth normalisation
    w_ref_perchi = (g_ref ** 3) * (1.0 / chi_h_arr) * (opz ** 3) * (Dz ** 4)

    # canoes per-shell scalar weight pieces (h-units, the deployed zeta integrand
    # AFTER the (u,v,phi) sum cancels b_delta; we keep b_delta in both and divide).
    resp_PPP = np.array([
        _kappa3_limber_response_product_h("PPP", float(poisson_amp[i]), float(opz[i]))
        for i in range(len(opz))
    ])
    # The deployed fold multiplies zeta by K_g^3 in lambda; zeta carries
    # growth^4 * b_delta / chi^4 * resp * rf * h^6.  Per-shell (b_delta divided out):
    # also the lambda-measure: INT dlambda = INT (a^2) dchi, but stage1_ours folds in
    # lambda directly with K_g (phys).  We compare the chi-integrated weights with the
    # SAME measure by folding canoes in lambda and ref in chi, then dividing totals.
    h6 = cosmo.h ** 6

    Ls = [60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 500.0, 800.0, 1000.0]
    print("=" * 104)
    print("Per-triangle folded ratio  W_can_fold / B_kappa  (equilateral L), and "
          "per-shell radial ratio L-independence check.")
    print("=" * 104)
    print(f"{'L':>6s} {'W_can_fold':>13s} {'B_kappa':>13s} {'R=Wc/Bk':>11s} "
          f"{'R/R0':>9s}")
    R0 = None
    rows = []
    perchi_ratio_by_L = {}
    for L in Ls:
        c12 = c13 = c23 = -0.5  # equilateral
        bdelta = np.zeros_like(chi_arr)
        for il, cv in enumerate(chi_arr):
            cf = float(cv)
            k = L / cf
            p = float(pk(k))
            f = float(_f2_kernel_numpy(np.array([k]), np.array([k]), np.array([-0.5]))[0])
            bdelta[il] = 2.0 * (f * p * p + f * p * p + f * p * p)
        # canoes per-shell zeta integrand (lambda density), b_delta included
        zeta_can = (growth ** 4) * bdelta / (chi_arr ** 4) * resp_PPP * rf * h6
        Wcan_fold = float(np.trapezoid((Kg ** 3) * zeta_can, lam_phys))
        # reference: B_kappa folded in chi (its native measure)
        bkv = float(bk(L, L, L))
        R = Wcan_fold / bkv
        if R0 is None:
            R0 = R
        rows.append((L, Wcan_fold, bkv, R))
        # per-shell radial ratio (b_delta cancels): canoes-vs-ref radial factor.
        # canoes radial: K_g^3 (in lambda fold = *a^2 dchi) * growth^4/chi^4*resp*rf*h6
        # ref radial:    g^3 * (1/chi) * (1+z)^3 * D^4
        canoes_radial = (Kg ** 3) * (a ** 2) * (growth ** 4) / (chi_arr ** 4) * resp_PPP * rf * h6
        ref_radial = w_ref_perchi
        perchi_ratio_by_L[L] = canoes_radial / ref_radial
        print(f"{L:6.0f} {Wcan_fold:13.4e} {bkv:13.4e} {R:11.4e} {R/R0:9.4f}")

    La = np.array([r[0] for r in rows])
    Ra = np.array([r[3] for r in rows])
    slope = np.polyfit(np.log(La), np.log(np.abs(Ra)), 1)[0]
    print(f"\n  log-log slope d ln|R| / d ln L = {slope:+.4f}")

    # Is the PER-SHELL radial ratio L-independent? If yes, the +slope comes ONLY
    # from b_delta(l/chi) being weighted differently across chi by the two
    # measures (a chi-reweighting), NOT from a per-leg ell factor.
    print("\nPer-shell radial ratio canoes/ref at chi index 20 (mid-LOS), across L:")
    for L in Ls:
        print(f"  L={L:6.0f}:  ratio[20] = {perchi_ratio_by_L[L][20]:.6e}")
    # they should be IDENTICAL across L (b_delta divided out) -> confirms the
    # radial normalisation is a pure (L-independent) factor; the +0.44 then comes
    # from the chi-WEIGHTING difference acting on b_delta(l/chi).
    ratios_at_20 = np.array([perchi_ratio_by_L[L][20] for L in Ls])
    print(f"  spread across L: max/min = {ratios_at_20.max()/ratios_at_20.min():.6f} "
          f"(==1 => radial factor is L-independent)")

    np.savez(
        OUT / "Wcan_vs_Bkappa_decomp.npz",
        L=La, W_can_fold=np.array([r[1] for r in rows]),
        B_kappa=np.array([r[2] for r in rows]), R=Ra, slope=slope,
        perchi_ratio_at_20=ratios_at_20,
        chi_arr=chi_arr, Kg=Kg, growth=growth, resp_PPP=resp_PPP,
        g_ref=g_ref, opz=opz,
        note="W_can folded vs B_kappa; per-shell radial ratio L-independence check.",
    )
    print(f"\nsaved -> {OUT / 'Wcan_vs_Bkappa_decomp.npz'}")


if __name__ == "__main__":
    main()
