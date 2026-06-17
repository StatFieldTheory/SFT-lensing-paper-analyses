#!/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
"""NUMERICAL closure check: apply the derived photon-energy reduction factor
R(z) = -1/[(1+z)^4 * h^2] to the REAL canoes zeta_TTT 3-pt fold, per-shell
inside the radial integral, and verify the canoes-reduced/PyCCL ratio drops
from ~10x to ~1 (within tree-Limber tolerance).

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python

WHAT IS APPLIED (exactly)
-------------------------
The derivation (scripts/docs/photon_energy_reduction.md, VERDICT) gives the
per-shell, per-source-z reduction factor for the four equal-shell kappa3
builders:

    R(z) = -1 / [ (1+z)^4 * h^2 ]

split into:
  * radial, PER SHELL : divide the equal-shell zeta by (1+z)^4 at each source
    shell z  (drops the 2 uncancelled local-E (1+z)^4 legs the equal-shell
    delta-collapse fails to return, the 3-leg analogue of the 2-pt projection
    reducing (1+z)^4 -> (1+z)^1);
  * units, CONSTANT   : divide by h^2 (h^6 "physical" footing -> h^4
    dimensionless-correct footing);
  * sign kept (R<0): flips the canoes-negative zeta to the standard-positive
    convergence-3PCF sign.

It is applied INSIDE the radial fold, per shell, as a function of z:

    Z_reduced(gamma) = INT dlambda  K(lambda,lambda_s)^3  [ R(z(lambda)) * zeta_TTT(gamma,lambda) ]

i.e. R(z) multiplies the per-shell zeta_TTT before the K^3 dlambda integral.
The radial 1/(1+z)^4 is genuinely per-shell (z-dependent); the 1/h^2 is a
constant. NOTHING else changes: K, zeta_TTT, the lambda grid are the deployed
stage0 values.

ORACLE
------
PyCCL standard tree conv-3PCF Z_std(gamma) from
outputs/scalar_vs_pyccl_tree.npz (built on byte-identical angular machinery;
the only difference vs canoes is the radial response convention, so the ratio
isolates exactly the reduction this script applies).

Output: outputs/scalar_3pt_reduced_vs_pyccl.npz
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ANALYSIS = _HERE.parent
_OUT = _ANALYSIS / "outputs"


def main() -> int:
    # --- load the REAL canoes per-shell zeta_TTT fold inputs (deployed stage0) --
    st0 = np.load(_OUT / "stage0_ours.npz", allow_pickle=True)
    gamma_arcmin = st0["gamma_arcmin"]
    z = st0["z"]                       # (40,) source-shell redshifts
    chi_phys = st0["chi_phys"]         # (40,) physical Mpc
    lam_phys = st0["lam_phys"]         # (40,) physical Mpc (affine)
    zeta_TTT = st0["zeta_TTT"]         # (n_gamma, 40) deployed equal-shell 3-cumulant
    Z_lambda_deployed = st0["Z_lambda"]  # (n_gamma,) deployed canoes fold (unreduced)
    chi_s = float(st0["chi_s"])
    h = float(st0["h"])
    z_s = float(st0["z_s"])

    # --- PyCCL oracle (byte-identical angular machinery; radial-convention only) -
    orc = np.load(_OUT / "scalar_vs_pyccl_tree.npz", allow_pickle=True)
    Z_std = orc["Z_std"]
    ga_orc = orc["gamma_arcmin"]
    assert np.allclose(ga_orc, gamma_arcmin), "gamma grid mismatch oracle vs stage0"
    # sanity: the stored real-canoes fold must equal the deployed stage0 fold
    Z_real_orc = orc["Z_real_canoes"]
    assert np.allclose(Z_real_orc, Z_lambda_deployed, rtol=1e-6), \
        "stored Z_real_canoes != deployed stage0 Z_lambda"

    # --- the convergence kernel K (physical Mpc), same as stage0_ours -----------
    a = 1.0 / (1.0 + z)
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s  # (40,)

    # ============================================================================
    # APPLY R(z) = -1/[(1+z)^4 h^2] PER SHELL, INSIDE the radial fold.
    # ============================================================================
    R_z = -1.0 / ((1.0 + z) ** 4 * h ** 2)              # (40,), per-shell
    zeta_reduced = R_z[None, :] * zeta_TTT              # (n_gamma, 40)
    integrand = (K ** 3)[None, :] * zeta_reduced        # (n_gamma, 40)
    Z_reduced = np.trapezoid(integrand, lam_phys, axis=1)  # (n_gamma,)

    # --- ALSO record the constant-per-shell variant: apply R as a single scalar
    #     evaluated at the lensing-weighted z_eff (sanity; the per-shell route is
    #     the principled one). Reported only for comparison.
    z_eff = orc["z_eff"]                                # (n_gamma,) from oracle
    R_eff = -1.0 / ((1.0 + z_eff) ** 4 * h ** 2)        # (n_gamma,)
    Z_reduced_zeff = R_eff * Z_lambda_deployed          # scalar-R applied to deployed fold

    # --- ratios -----------------------------------------------------------------
    eps = 1e-300
    ratio_reduced = Z_reduced / np.where(np.abs(Z_std) > eps, Z_std, eps)  # signed
    ratio_reduced_abs = np.abs(Z_reduced) / np.maximum(np.abs(Z_std), eps)
    ratio_unreduced_abs = np.abs(Z_lambda_deployed) / np.maximum(np.abs(Z_std), eps)
    ratio_zeff_abs = np.abs(Z_reduced_zeff) / np.maximum(np.abs(Z_std), eps)

    # sign agreement: reduced and std should now share sign (both > 0)
    sign_ok = np.sign(Z_reduced) == np.sign(Z_std)

    print("=" * 104)
    print(f"3-pt PHOTON-ENERGY REDUCTION CLOSURE  (z_s={z_s}, h={h})")
    print("Applied: R(z) = -1/[(1+z)^4 h^2] PER SHELL inside  Z = INT dlambda K^3 [R(z) zeta_TTT]")
    print("=" * 104)
    hdr = (f"{'gamma':>9} {'Z_std(pyccl)':>13} {'Z_reduced':>13} "
           f"{'unred/std':>10} {'REDUCED/std':>12} {'(signed)':>10} "
           f"{'zeff-R/std':>11} {'sign_ok':>8}")
    print(hdr)
    for i in range(gamma_arcmin.size):
        print(f"{gamma_arcmin[i]:9.3f} {Z_std[i]:13.4e} {Z_reduced[i]:13.4e} "
              f"{ratio_unreduced_abs[i]:10.3f} {ratio_reduced_abs[i]:12.3f} "
              f"{ratio_reduced[i]:10.3f} {ratio_zeff_abs[i]:11.3f} "
              f"{str(bool(sign_ok[i])):>8}")

    # --- closure summary over the clean signal regime (gamma <= ~30') -----------
    sig = gamma_arcmin <= 30.0
    print()
    print("-" * 104)
    print("CLOSURE SUMMARY (signal regime gamma <= 30', where tree-Limber + lensing")
    print("kernel are well-behaved; gamma > 60' is sign-marginal / Limber-tail noise):")
    print(f"  unreduced |canoes/PyCCL| : "
          f"min={ratio_unreduced_abs[sig].min():.3f}  "
          f"max={ratio_unreduced_abs[sig].max():.3f}  "
          f"median={np.median(ratio_unreduced_abs[sig]):.3f}")
    print(f"  REDUCED   |canoes/PyCCL| : "
          f"min={ratio_reduced_abs[sig].min():.3f}  "
          f"max={ratio_reduced_abs[sig].max():.3f}  "
          f"median={np.median(ratio_reduced_abs[sig]):.3f}")
    print(f"  REDUCED   sign agreement : {int(sign_ok[sig].sum())}/{int(sig.sum())} shells")
    med = float(np.median(ratio_reduced_abs[sig]))
    closed = (0.7 <= med <= 1.3)
    print(f"  VERDICT: median reduced ratio = {med:.3f} "
          f"-> {'CLOSES (within ~30%)' if closed else 'does NOT close'}")

    np.savez(
        _OUT / "scalar_3pt_reduced_vs_pyccl.npz",
        gamma_arcmin=gamma_arcmin,
        z=z, chi_phys=chi_phys, lam_phys=lam_phys, K=K,
        zeta_TTT=zeta_TTT, R_z=R_z, h=h, z_s=z_s,
        Z_std=Z_std,
        Z_lambda_deployed=Z_lambda_deployed,
        Z_reduced=Z_reduced,
        Z_reduced_zeff=Z_reduced_zeff,
        ratio_unreduced_abs=ratio_unreduced_abs,
        ratio_reduced_abs=ratio_reduced_abs,
        ratio_reduced_signed=ratio_reduced,
        ratio_zeff_abs=ratio_zeff_abs,
        sign_ok=sign_ok,
        z_eff=z_eff,
        signal_mask=sig,
        note="canoes zeta_TTT fold reduced by R(z)=-1/[(1+z)^4 h^2] applied per "
             "shell inside the radial fold (radial 1/(1+z)^4 per shell + units "
             "1/h^2 constant, sign kept), vs PyCCL standard tree conv-3PCF Z_std. "
             "Closure = reduced/PyCCL ~ 1 within tree-Limber tolerance.",
    )
    print(f"\nsaved -> {_OUT / 'scalar_3pt_reduced_vs_pyccl.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
