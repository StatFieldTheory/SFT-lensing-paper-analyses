"""EXTERNAL TEXTBOOK tree-level convergence 3PCF vs canoes equal-shell zeta fold.

Interpreter (ABSOLUTE; subagent Bash cannot conda activate):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

PURPOSE (2026-06-09 external check; non-same-family)
----------------------------------------------------
Build a GENUINELY INDEPENDENT textbook tree-level CONVERGENCE 3PCF and compare it
to the canoes production equal-shell zeta_TTT fold (the value behind FK kk
+1.219e-4).

Why this is NOT same-family
---------------------------
``stage0_spt_reference.py`` is SAME-FAMILY: it hand-rolls the Sachs (1+z)^4 zeta
(response = (A(a)(1+z)^4)^3 per leg) and folds it with the SAME K^3 kernel as
``stage0_ours.py``.  So it carries the SAME (1+z)^4 over-count and cannot expose
it.  Here we instead build the STANDARD textbook convergence bispectrum from the
ground up with the TEXTBOOK lensing kernel

    W(chi) = (3/2) Omega_m (H0/c)^2 (1+z) chi (chi_s - chi)/chi_s    [ONE (1+z)/leg]

and the SPT tree-level matter bispectrum B_delta = 2 F2 P P + cyc (growth D^4),
projected to the squeezed/collapsed equal-shell config via the same spin-0 J0
orientation kernel.  The convergence bispectrum is

    b^kappa_{l1 l2 l3} = INT dchi (W(chi)^3 / chi^4) B_delta(l1/chi, l2/chi, l3/chi; chi),

and the squeezed real-space zeta_TTT(gamma) at fixed source plane is the flat-sky
J0 transform of that integrand on the collapsed triangle.  This is the SAME
physics fastnc uses for the shear 3PCF (external ground truth).

THE CANOES PRODUCTION FOLD (raw, no user (1+z)^-4 recipe)
--------------------------------------------------------
    Z_canoes(gamma) = INT_0^{lambda_s} dlambda  K(lambda, lambda_s)^3  zeta_TTT^canoes(gamma, lambda)
    K(lambda, lambda_s) = a^2(chi) chi (chi_s - chi)/chi_s,   dlambda = a^2 dchi.
zeta_TTT^canoes is built fresh by the canonical equal_time_limber recipe
(LOW exact-3j ell<=60 + HIGH flat-sky Limber (60,1000], tree_phi, units=physical,
radial_measure=lambda, channels=TTT).  This is EXACTLY the production fold that
feeds FK kk -- we deliberately DROP the analysis-side (1+z)^-4 born reduction so
the ratio exposes whether the raw production object over-counts.

OUTPUT
------
For each source plane z_s in a scan, report
    ratio(gamma) = Z_canoes(gamma) / Z_textbook(gamma)
and fit log|ratio| vs log(1+z_s) to extract the (1+z) power.  A power ~ +4 (flat
in gamma) confirms the canoes 3-pt over-counts by (1+z)^4 (user (1+z)^-4 recipe
needed, production FK affected).  A power ~ 0 with ratio ~ 1 means canoes is
already correct.

Result printed as JSON to stdout (line prefixed 'JSON_RESULT: ').
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.special import jv
from scipy.integrate import cumulative_trapezoid

_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

_HERE = Path(__file__).resolve().parent
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402

from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
from canoes.cosmo.background import C_KM_S  # noqa: E402

# Source planes to scan.  Span a wide (1+z) lever arm so the power fit is clean.
Z_SOURCE_SCAN = [0.5, 1.0, 2.0, 3.0, 5.0]
_N_LAMBDA = 40        # source-shell grid per plane
_GAMMA_ARCMIN = np.unique(np.concatenate(
    [np.array([1.0]), np.geomspace(0.5, 300.0, 9)]))

# Textbook flat-sky quadrature (independent of canoes' Nmax/ell_max).
_N_ELL = 220
_N_PHI = 96
_ELL_MIN = 1e-3
_ELL_MAX = 1000.0     # match canoes LOW+HIGH cap so we compare the SAME ell band


# ---------------------------------------------------------------------------
# Background / affine grid (identical conventions to stage0_ours)
# ---------------------------------------------------------------------------
def _build_grid(cosmo, z_source: float):
    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, z_source, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])  # physical Mpc
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_source / _N_LAMBDA, z_source, _N_LAMBDA)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    lam_h = lam_phys * h
    return z, chi_phys, lam_phys, lam_h, float(chif[-1]), float(lam_f[-1])


def _growth_canoes(z: np.ndarray, cosmo) -> np.ndarray:
    """Linear growth D(z) consistent with the canoes tables (M[0]/(1+z))."""
    from canoes.sachs.kappa3 import build_lightcone  # type: ignore
    chi_h = np.array([chi_of_z(zi, Omega_m=cosmo.Omega_m, h=cosmo.h)
                      for zi in z]) * cosmo.h
    lightcone, one_plus_z = build_lightcone(cosmo, chi_h)
    D = lightcone.M[0, :] / one_plus_z
    return np.asarray(D, dtype=np.float64)


def _f2_spt(k_a, k_b, mu):
    safe_a = np.where(k_a > 0, k_a, 1.0)
    safe_b = np.where(k_b > 0, k_b, 1.0)
    ratio = safe_a / safe_b + safe_b / safe_a
    return 5.0 / 7.0 + 0.5 * mu * ratio + (2.0 / 7.0) * mu * mu


def _squeezed_triples(gamma_rad: np.ndarray) -> np.ndarray:
    cg = np.cos(gamma_rad)
    return np.stack([np.ones_like(cg), cg, cg], axis=1)


# ---------------------------------------------------------------------------
# (1) CANOES production equal-shell zeta_TTT fold (raw, no (1+z)^-4 recipe)
# ---------------------------------------------------------------------------
def canoes_fold(cosmo, pk, z, chi_phys, lam_phys, lam_h, chi_s, gamma_rad):
    from canoes.sachs import (
        compute_kappa3_sigma3_high,
        compute_kappa3_zeta_table,
        kappa3_combine_low_high,
    )
    triples = _squeezed_triples(gamma_rad)
    common_low = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_max=60, Nmax=128, spt_kind="tree_phi",
        units="physical", lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda", channels=("TTT",), parallelism="sequential",
        shell_block_size=4, triple_chunk_size=4096, n_workers=1,
        fuse_tree_phi_channels=True,
    )
    common_high = dict(
        pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000, ell_min=1e-3,
        n_ell=96, n_phi=64, channels=("TTT",), units="physical",
        lambda_convention="project", radial_discretization="sample",
        radial_measure="lambda",
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low = compute_kappa3_zeta_table(triples, lam_h, **common_low)
        high = compute_kappa3_sigma3_high(triples, lam_h, **common_high)
    combined = kappa3_combine_low_high(low, high)
    zeta_TTT = np.asarray(combined.zeta_TTT, dtype=np.float64)  # (n_gamma, n_lambda)

    a = 1.0 / (1.0 + z)
    # RAW production fold: K^3 zeta, NO (1+z)^-4 born reduction.
    K = a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s  # physical Mpc
    integrand = (K ** 3)[None, :] * zeta_TTT
    Z = np.trapezoid(integrand, lam_phys, axis=1)
    return Z, zeta_TTT


# ---------------------------------------------------------------------------
# (2) TEXTBOOK tree-level convergence 3PCF (standard lensing kernel W^3/chi^4)
# ---------------------------------------------------------------------------
def textbook_fold(cosmo, pk, z, chi_phys, D, chi_s, gamma_rad):
    """Standard tree-level convergence 3PCF, squeezed/collapsed equal-shell.

    b^kappa = INT dchi (W^3/chi^4) B_delta,  W = (3/2) Om (H0/c)^2 (1+z) g(chi),
    g(chi) = chi (chi_s-chi)/chi_s.  Project to zeta_TTT(gamma) with the SAME
    spin-0 J0 orientation kernel as the reference, so only the KERNEL differs
    (textbook W per leg vs Sachs (A(1+z)^4) per leg).  All h-units internally;
    convert to physical by h^4 (matter convergence 3PCF native dim (h/Mpc)^4).
    """
    h = cosmo.h
    Om = cosmo.Omega_m
    # H0 in h-units (h/Mpc): H0 = 100 h km/s/Mpc -> (H0/c) = 100/c per (h/Mpc).
    H0_over_c_hunits = 100.0 / C_KM_S  # = (H0/c) in (h/Mpc)
    chi_phys_h = chi_phys * h
    chi_s_h = chi_s * h
    one_plus_z = 1.0 + z
    a = 1.0 / one_plus_z

    # textbook lensing kernel W(chi), h-units (per leg ONE (1+z))
    g = chi_phys_h * (chi_s_h - chi_phys_h) / chi_s_h     # geometric efficiency [Mpc/h]
    W = 1.5 * Om * (H0_over_c_hunits ** 2) * one_plus_z * g  # dimensionless density [h/Mpc]

    ell = np.geomspace(_ELL_MIN, _ELL_MAX, _N_ELL)
    phi = np.linspace(0.0, 2.0 * np.pi, _N_PHI, endpoint=False)
    dphi = 2.0 * np.pi / _N_PHI
    u = ell[:, None, None]
    v = ell[None, :, None]
    ph = phi[None, None, :]
    cph = np.cos(ph)
    w = np.sqrt(np.maximum(u * u + v * v + 2.0 * u * v * cph, 0.0))
    safe_w = np.maximum(w, _ELL_MIN)
    cos_12 = -(u + v * cph) / safe_w
    cos_13 = -(v + u * cph) / safe_w
    cos_23 = cph

    dlnell = np.gradient(np.log(ell))
    w_u = (ell * dlnell)[:, None, None]
    w_v = (ell * dlnell)[None, :, None]
    base_meas = w_u * w_v * dphi * u * v
    prefactor = 1.0 / ((2.0 * np.pi) ** 4)

    n_g = gamma_rad.size
    n_z = z.size
    # zeta_TTT-equivalent per-shell density built the SAME way as the reference
    # (spin-0 J0 squeezed), but with the textbook W per leg in place of canoes'
    # Sachs response.  Native h-units; physical conversion h^4 by caller.
    zeta = np.zeros((n_g, n_z), dtype=np.float64)
    for iz in range(n_z):
        chi = float(chi_phys_h[iz])
        k1 = safe_w / chi
        k2 = np.maximum(u, _ELL_MIN) / chi
        k3 = np.maximum(v, _ELL_MIN) / chi
        p1 = np.asarray(pk(k1), dtype=np.float64)
        p2 = np.asarray(pk(k2), dtype=np.float64)
        p3 = np.asarray(pk(k3), dtype=np.float64)
        f12 = _f2_spt(k1, k2, cos_12)
        f13 = _f2_spt(k1, k3, cos_13)
        f23 = _f2_spt(k2, k3, cos_23)
        b_delta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
        # growth: B_delta carries D^4 (tree, three legs).  D^4 = (D^2)^2; the
        # P(k) table is z=0 so multiply by D^4 explicitly.
        radial = prefactor * (D[iz] ** 4) / chi ** 4
        for ig in range(n_g):
            ang = 2.0 * np.pi * jv(0, v * gamma_rad[ig])
            integrand = base_meas * b_delta * ang
            zeta[ig, iz] = radial * float(np.sum(integrand))

    zeta_phys = zeta * (h ** 4)  # (h/Mpc)^4 native -> physical

    # textbook fold: Z = INT dchi (W^3/chi^4) ... but we have already pulled the
    # 1/chi^4 and W-less B_delta into zeta_phys (the per-shell collapsed density).
    # The convergence 3PCF is INT dchi W(chi)^3 * [zeta_phys / (per-leg-norm)].
    # zeta_phys here is the MATTER squeezed 3PCF density (B_delta J0 / chi^4) with
    # growth; multiplying by W^3 (in physical Mpc) and dchi gives the convergence
    # 3PCF.  W is dimensionless density [h/Mpc] -> use physical chi for dchi.
    W_phys = W * h  # convert W from [h/Mpc] to [1/Mpc] so W^3 dchi[Mpc] is clean
    # NOTE: zeta_phys already carries the matter bispectrum shape; W_phys^3 is the
    # textbook lensing weight per leg.  dchi in physical Mpc.
    integrand = (W_phys ** 3)[None, :] * zeta_phys
    Z = np.trapezoid(integrand, chi_phys, axis=1)
    return Z, zeta_phys


def main() -> int:
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)

    gamma_arcmin = _GAMMA_ARCMIN
    gamma_rad = np.radians(gamma_arcmin / 60.0)
    i1 = int(np.argmin(np.abs(gamma_arcmin - 1.0)))

    rows = []
    for z_s in Z_SOURCE_SCAN:
        z, chi_phys, lam_phys, lam_h, chi_s, lam_s = _build_grid(cosmo, z_s)
        D = _growth_canoes(z, cosmo)
        Zc, _ = canoes_fold(cosmo, pk, z, chi_phys, lam_phys, lam_h, chi_s, gamma_rad)
        Zt, _ = textbook_fold(cosmo, pk, z, chi_phys, D, chi_s, gamma_rad)
        ratio = Zc / np.where(np.abs(Zt) > 0, Zt, np.nan)
        rows.append(dict(
            z_s=z_s,
            Z_canoes=Zc.tolist(),
            Z_textbook=Zt.tolist(),
            ratio=ratio.tolist(),
            ratio_at_1arcmin=float(ratio[i1]),
            absratio_median=float(np.nanmedian(np.abs(ratio))),
        ))
        print(f"[scan] z_s={z_s:.2f}  Zc(1')={Zc[i1]:+.4e}  Zt(1')={Zt[i1]:+.4e}  "
              f"ratio(1')={ratio[i1]:+.4e}  |ratio|_med={np.nanmedian(np.abs(ratio)):.4e}",
              flush=True)

    # --- fit log|ratio(1')| vs log(1+z_s): the (1+z) power -------------------
    zs = np.array([r["z_s"] for r in rows])
    r1 = np.array([abs(r["ratio_at_1arcmin"]) for r in rows])
    rmed = np.array([r["absratio_median"] for r in rows])
    x = np.log(1.0 + zs)
    p1 = np.polyfit(x, np.log(r1), 1)
    pmed = np.polyfit(x, np.log(rmed), 1)
    power_1arcmin = float(p1[0])
    power_median = float(pmed[0])
    # amplitude at z=0 (extrapolated) so we can see the constant offset:
    amp0_1arcmin = float(np.exp(p1[1]))
    amp0_median = float(np.exp(pmed[1]))

    # Reference (1+z)^4 expectation: pin to the z_s=5 ratio and see (1+z)^-4 collapse.
    i5 = int(np.argmin(np.abs(zs - 5.0)))
    ratio5 = r1[i5]
    onepz5 = 1.0 + zs[i5]
    reduced5 = ratio5 / (onepz5 ** 4)  # if ~1, the over-count is exactly (1+z)^4

    result = dict(
        gamma_arcmin=gamma_arcmin.tolist(),
        z_source_scan=zs.tolist(),
        ratio_at_1arcmin_per_zs={f"{z:.2f}": float(r) for z, r in zip(zs, r1)},
        absratio_median_per_zs={f"{z:.2f}": float(r) for z, r in zip(zs, rmed)},
        onepz_power_fit_1arcmin=power_1arcmin,
        onepz_power_fit_median=power_median,
        amp0_extrapolated_1arcmin=amp0_1arcmin,
        amp0_extrapolated_median=amp0_median,
        ratio_at_zs5_1arcmin=float(ratio5),
        ratio_div_1pz4_at_zs5=float(reduced5),
        rows=rows,
    )

    out = _HERE.parent / "outputs" / "textbook_3pt_external_check.npz"
    np.savez(out,
             gamma_arcmin=gamma_arcmin, z_source_scan=zs,
             ratio_1arcmin=r1, absratio_median=rmed,
             power_1arcmin=power_1arcmin, power_median=power_median)
    print(f"[saved] -> {out}", flush=True)
    print()
    print(f"[fit] log|ratio(1')| vs log(1+z_s): power = {power_1arcmin:+.3f}  (amp0={amp0_1arcmin:.3e})")
    print(f"[fit] log|ratio|_median   vs log(1+z_s): power = {power_median:+.3f}  (amp0={amp0_median:.3e})")
    print(f"[anchor] ratio(1', z_s=5) = {ratio5:.4e};  ratio/(1+z_s)^4 = {reduced5:.4e}  "
          f"(==1 => over-count is exactly (1+z)^4)")
    print()
    print("JSON_RESULT: " + json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
