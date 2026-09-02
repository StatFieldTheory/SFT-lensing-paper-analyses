"""End-to-end normalisation test of the kappa3 vertex against standard Limber.

This is the decisive Step-1 test of the normalisation study.  Unlike
``verify_chain.py``, which re-implements the pipeline's factors and checks
they assemble correctly, this module CALLS the production function
``canoes.sachs.kappa3.compute_kappa3_sigma3_high`` and confronts its output
with an independently coded standard-Limber convergence three-point
function.  A coding slip anywhere inside the production routine (a stray
factor of two, a wrong power of ``h``, a misnormalised orientation
integral) shows up here as a ratio away from one.

The two sides
-------------
Pipeline side, per shell, for the TTT channel:

    f_pipe(chi) = -(d lambda/d chi) K(lambda)^3 zeta_TTT^lambda(chi)

with ``zeta_TTT^lambda`` exactly as returned by the production routine
(``radial_measure='lambda'``), ``K`` the paper's affine lensing window
[cosmology.tex Eq. "K affine kernel"], and the leading minus sign the cube
of ``kappa = - int d lambda K Phi_00`` [Eq. "kappa from Dsachs1"].

Reference side, per shell:

    f_ref(chi) = W(chi)^3 / chi^4 * I(chi)

    W(chi)  = 3/2 Omega_m H0^2 (1+z) chi (chi_s - chi) / chi_s
    I(chi)  = (2pi)^-4 Re int d^2 l2 d^2 l3 W_high
                  B_delta(|l2+l3|/chi, l2/chi, l3/chi; z)
                  exp(i [l2.r2 + l3.r3])

``I`` is coded here from scratch: the absolute-orientation integral is done
by an explicit periodic quadrature over the angle alpha rather than by the
analytic ``2 pi J_0(|Q|)`` the production routine uses, so the reference
does not inherit the production alpha kernel.  The radial GL nodes, the
window and the matter bispectrum are deliberately shared, so quadrature
error cancels and only the NORMALISATION is compared.

If the chain is right, ``f_pipe / f_ref = 1`` shell by shell and triple by
triple.

Run:
  CAN=/Users/zzhang/projects/angular_statistics/canoes
  cd sachs_sft/analyses/fk_audit_2026-08
  PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python normalisation/bkappa_endtoend.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from canoes.cosmo.background import H_DEFAULT, OMEGA_M_DEFAULT
from canoes.cosmo.lightcone import build_lightcone, lambda_of_chi
from canoes.sachs.kappa3 import (
    _H0_H_PER_MPC,
    _kappa3_flat_triangle_vertices,
    _validate_limber_high_quadrature,
    compute_kappa3_sigma3_high,
)

ARCMIN = math.pi / (180.0 * 60.0)


@dataclass(frozen=True)
class SimpleCosmo:
    """Minimal ``CosmoLike``: canoes only reads ``Omega_m`` and ``h``."""

    Omega_m: float
    h: float


COSMO = SimpleCosmo(Omega_m=OMEGA_M_DEFAULT, h=H_DEFAULT)


# ------------------------------------------------------------------ geometry
def chi_of_z(cosmo: SimpleCosmo, z_target: float) -> float:
    """Invert the canoes lightcone for comoving distance in Mpc/h."""
    grid = np.linspace(1.0, 12000.0, 60001)
    _lc, opz = build_lightcone(cosmo, grid, distance_unit="Mpc/h")
    return float(np.interp(1.0 + z_target, opz, grid))


def affine_window_K(
    cosmo: SimpleCosmo, chi_h: NDArray[np.float64], chi_s_h: float
) -> NDArray[np.float64]:
    """``K(lambda, lambda_s) = a^2 chi (chi_s - chi) / chi_s``, Mpc/h.

    Verified in ``verify_chain.py`` T2 against brute-force quadrature of the
    defining form ``Dbar^2 int_lambda^{lambda_s} d lambda_1 / Dbar^2``.
    """
    chi_h = np.asarray(chi_h, dtype=np.float64)
    _lc, opz = build_lightcone(cosmo, chi_h, distance_unit="Mpc/h")
    return (chi_h * (float(chi_s_h) - chi_h) / float(chi_s_h)) / opz**2


def lensing_efficiency_W(
    cosmo: SimpleCosmo, chi_h: NDArray[np.float64], chi_s_h: float
) -> NDArray[np.float64]:
    """``W = 3/2 Omega_m H0^2 (1+z) chi (chi_s-chi)/chi_s`` in h-native units."""
    chi_h = np.asarray(chi_h, dtype=np.float64)
    _lc, opz = build_lightcone(cosmo, chi_h, distance_unit="Mpc/h")
    return (
        1.5
        * float(cosmo.Omega_m)
        * _H0_H_PER_MPC**2
        * opz
        * chi_h
        * (float(chi_s_h) - chi_h)
        / float(chi_s_h)
    )


# ----------------------------------------------------- independent reference
def reference_angular_integral(
    *,
    cosine_triple: tuple[float, float, float],
    chi_h: float,
    z: float,
    b_delta_fn,
    ell_cut: int,
    ell_high_max: int,
    ell_min: float,
    n_ell: int,
    n_phi: int,
    n_alpha: int = 512,
) -> float:
    """``I(chi)``: the flat-sky scalar 3PCF integral, orientation done by hand.

    The absolute orientation ``alpha`` is integrated with the periodic
    trapezoid rule, which is spectrally accurate for the smooth periodic
    integrand ``exp(i |Q| cos(alpha - arg Q))``.
    """
    ell_nodes, ell_w, phi_nodes, phi_w = _validate_limber_high_quadrature(
        ell_cut=ell_cut,
        ell_high_max=ell_high_max,
        ell_min=ell_min,
        n_ell=n_ell,
        n_phi=n_phi,
    )
    u = ell_nodes[:, None, None]
    v = ell_nodes[None, :, None]
    phi = phi_nodes[None, None, :]
    w = np.sqrt(np.maximum(u * u + v * v + 2.0 * u * v * np.cos(phi), 0.0))
    window = ((np.maximum(np.maximum(u, v), w) > float(ell_cut))
              & (np.maximum(np.maximum(u, v), w) <= float(ell_high_max)))
    measure = (
        ell_w[:, None, None] * ell_w[None, :, None] * phi_w[None, None, :]
        * u * v * window
    )
    safe_w = np.maximum(w, ell_min)
    safe_u = np.maximum(u, ell_min)
    safe_v = np.maximum(v, ell_min)
    b_val = np.asarray(
        b_delta_fn(safe_w / chi_h, safe_u / chi_h, safe_v / chi_h, z),
        dtype=np.float64,
    )

    t12 = math.acos(float(np.clip(cosine_triple[0], -1.0, 1.0)))
    t23 = math.acos(float(np.clip(cosine_triple[1], -1.0, 1.0)))
    t31 = math.acos(float(np.clip(cosine_triple[2], -1.0, 1.0)))
    r2, r3 = _kappa3_flat_triangle_vertices(t12, t23, t31)

    weighted = measure * b_val
    alpha = 2.0 * math.pi * np.arange(n_alpha) / n_alpha
    dalpha = 2.0 * math.pi / n_alpha
    acc = 0.0
    for a in alpha:
        # l2 = u e^{i a}, l3 = v e^{i(a+phi)};  l.r = Re(conj(l) r)
        phase = (
            u * np.real(np.exp(-1j * a) * r2)
            + v * np.real(np.exp(-1j * (a + phi)) * r3)
        )
        acc += float(np.sum(weighted * np.cos(phase)))
    return acc * dalpha / (2.0 * math.pi) ** 4


# ------------------------------------------------------------------- driver
def main() -> None:
    import b_model

    z_s = 1.0
    chi_s = chi_of_z(COSMO, z_s)
    ell_cut, ell_high_max = 60, 1000
    ell_min, n_ell, n_phi = 1e-3, 64, 32

    z_shells = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    chi_shells = np.array([chi_of_z(COSMO, float(z)) for z in z_shells])
    lam_shells = lambda_of_chi(
        COSMO, chi_shells, convention="project", distance_unit="Mpc/h"
    )

    gammas = np.array([10.0, 60.0]) * ARCMIN
    triples = []
    labels = []
    for g in gammas:
        triples.append((1.0, math.cos(g), math.cos(g)))
        labels.append(f"collapsed  gamma={g / ARCMIN:6.1f}'")
    for g in gammas:
        triples.append((math.cos(g), math.cos(g), math.cos(g)))
        labels.append(f"equilateral side={g / ARCMIN:6.1f}'")
    # a scalene open triangle: sides 20', 40', 50'
    a1, a2, a3 = 20.0 * ARCMIN, 40.0 * ARCMIN, 50.0 * ARCMIN
    triples.append((math.cos(a1), math.cos(a2), math.cos(a3)))
    labels.append("scalene    20'/40'/50'")
    cos_arr = np.array(triples, dtype=np.float64)

    b_delta_fn = b_model.make_b_delta_fn("tree")
    pk_lin = b_model.load_pk_lin(b_model.PK_TABLE_3PT)

    print(f"z_s = {z_s}, chi_s = {chi_s:.2f} Mpc/h; "
          f"ell in ({ell_cut}, {ell_high_max}], n_ell={n_ell}, n_phi={n_phi}")
    print("matter bispectrum: b_model 'tree' (canoes-identical), h-native units\n")

    out_h = compute_kappa3_sigma3_high(
        cos_arr,
        lam_shells,
        pk_matter_today=pk_lin,
        cosmo=COSMO,
        ell_cut=ell_cut,
        ell_high_max=ell_high_max,
        ell_min=ell_min,
        n_ell=n_ell,
        n_phi=n_phi,
        channels=("TTT",),
        units="h",
        radial_measure="lambda",
        b_delta_fn=b_delta_fn,
    )
    out_phys = compute_kappa3_sigma3_high(
        cos_arr,
        lam_shells,
        pk_matter_today=pk_lin,
        cosmo=COSMO,
        ell_cut=ell_cut,
        ell_high_max=ell_high_max,
        ell_min=ell_min,
        n_ell=n_ell,
        n_phi=n_phi,
        channels=("TTT",),
        units="physical",
        radial_measure="lambda",
        b_delta_fn=b_delta_fn,
    )

    chi_code = np.asarray(out_h.chi_shells, dtype=np.float64)
    _lc, opz = build_lightcone(COSMO, chi_code, distance_unit="Mpc/h")
    d_lam_d_chi = 1.0 / opz**2
    kk = affine_window_K(COSMO, chi_code, chi_s)
    ww = lensing_efficiency_W(COSMO, chi_code, chi_s)

    print(f"{'configuration':<24} {'z':>5} "
          f"{'f_pipe':>13} {'f_ref':>13} {'ratio':>12}")
    ratios = []
    for i_cos, label in enumerate(labels):
        for i_lam in range(chi_code.size):
            chi = float(chi_code[i_lam])
            z = float(opz[i_lam] - 1.0)
            zeta_lam = float(out_h.zeta_TTT[i_cos, i_lam])
            f_pipe = -d_lam_d_chi[i_lam] * kk[i_lam] ** 3 * zeta_lam
            i_ref = reference_angular_integral(
                cosine_triple=tuple(cos_arr[i_cos]),
                chi_h=chi,
                z=z,
                b_delta_fn=b_delta_fn,
                ell_cut=ell_cut,
                ell_high_max=ell_high_max,
                ell_min=ell_min,
                n_ell=n_ell,
                n_phi=n_phi,
            )
            f_ref = ww[i_lam] ** 3 / chi**4 * i_ref
            ratio = f_pipe / f_ref
            ratios.append(ratio)
            print(f"{label:<24} {z:5.2f} {f_pipe:13.5e} {f_ref:13.5e} "
                  f"{ratio:12.8f}")
    ratios_arr = np.array(ratios)
    print(f"\nratio f_pipe/f_ref: min {ratios_arr.min():.10f}  "
          f"max {ratios_arr.max():.10f}  "
          f"max |1-ratio| {np.abs(1.0 - ratios_arr).max():.3e}")

    # units cross-check: physical zeta must be h^4 times the h-native one
    r_units = out_phys.zeta_TTT / out_h.zeta_TTT
    print(f"physical/h-native zeta_TTT: {r_units.min():.10f} .. "
          f"{r_units.max():.10f}  (expected h^4 = {COSMO.h**4:.10f})")
    print(f"physical lambda_shells / h-native: "
          f"{float(np.max(out_phys.lambda_shells / out_h.lambda_shells)):.10f} "
          f"(expected 1/h = {1.0 / COSMO.h:.10f})")

    # the h-cancellation statement: the folded, dimensionless <kkk> must not
    # depend on which unit system the table is served in.
    lam_h = np.asarray(out_h.lambda_shells, dtype=np.float64)
    lam_p = np.asarray(out_phys.lambda_shells, dtype=np.float64)
    kk_p = kk / COSMO.h
    fold_h = np.trapezoid(-(kk**3) * out_h.zeta_TTT[0], lam_h)
    fold_p = np.trapezoid(-(kk_p**3) * out_phys.zeta_TTT[0], lam_p)
    print(f"folded <kkk> h-native {fold_h:.8e} vs physical {fold_p:.8e} "
          f"ratio {fold_h / fold_p:.12f}")


if __name__ == "__main__":
    main()
