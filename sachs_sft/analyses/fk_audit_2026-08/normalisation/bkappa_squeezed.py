"""Step 2: the pipeline's convergence bispectrum vs the published measurement.

Predicts B_kappa(l1, l2, l3) at z_s = 1 with the SFT pipeline's own transfer
chain and the BiHalofit matter bispectrum, in the cosmology of the
ray-tracing maps the published measurement was made on (Sato et al. 2009 /
Kayo et al. 2013: Omega_m = 0.238, Omega_b = 0.042, h = 0.732, n_s = 0.958,
sigma8 = 0.76).  The target is Figure 13 of Takahashi et al. 2020
(arXiv:1911.07886), which shows the same three shape families:

    equilateral   B_kappa(l, l, l)        plotted as l^3 B
    squeezed      B_kappa(l, l, 84)       plotted as 60 l^2 B
    squeezed      B_kappa(l, l, 510)      plotted as 300 l^2 B

Nothing textbook is assumed: every projection factor is taken from the
canoes vertex code itself.

    B_kappa = - int d lambda K(lambda)^3 * J(z) * R_TTT(z) / chi^4 * B_delta

      K        canoes / paper affine lensing window, a^2 chi (chi_s-chi)/chi_s
      J        _kappa3_radial_density_conversion(..., 'lambda') = (1+z)^-4
      R_TTT    _kappa3_limber_response_product_h('TTT', A(a), 1+z)
      leading minus sign: kappa = - int d lambda K Phi_00

``bkappa_endtoend.py`` shows that this chain, applied shell by shell to the
output of the production routine ``compute_kappa3_sigma3_high``, reproduces
that routine to 4e-15, so the assembly here IS the pipeline's projection.

Run:
  CAN=/Users/zzhang/projects/angular_statistics/canoes
  cd sachs_sft/analyses/fk_audit_2026-08
  PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python normalisation/bkappa_squeezed.py
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from canoes.cosmo.growth import growth_D
from canoes.cosmo.lightcone import build_lightcone
from canoes.sachs.kappa2 import _lambda_jacobian
from canoes.sachs.kappa3 import (
    _H0_H_PER_MPC,
    _kappa3_limber_response_product_h,
    _kappa3_radial_density_conversion,
)

DATA = Path(__file__).resolve().parent / "data"
PK_WMAP3 = DATA / "PCAMB_wmap3_sato2009_z0.txt"


@dataclass(frozen=True)
class SimCosmo:
    """Sato et al. 2009 / Kayo et al. 2013 ray-tracing cosmology."""

    Omega_m: float = 0.238
    Omega_b: float = 0.042
    h: float = 0.732
    n_s: float = 0.958
    sigma8: float = 0.76

    @property
    def Omega_de(self) -> float:
        return 1.0 - self.Omega_m


SIM = SimCosmo()


# --------------------------------------------------------------- background
def chi_of_z(cosmo, z_target: float) -> float:
    grid = np.linspace(1.0, 12000.0, 60001)
    _lc, opz = build_lightcone(cosmo, grid, distance_unit="Mpc/h")
    return float(np.interp(1.0 + z_target, opz, grid))


# ---------------------------------------------------------- matter spectrum
def load_pk_wmap3() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    arr = np.loadtxt(PK_WMAP3)
    return arr[:, 0].copy(), arr[:, 1].copy()


def make_bihalofit(cosmo: SimCosmo):
    """A fastnc BiHalofit instance for this cosmology and linear P(k)."""
    from b_model.halofit_backend import Z_SPLINE_MAX, load_fastnc_halofit

    module = load_fastnc_halofit()
    k_h, p_m = load_pk_wmap3()
    z_grid = np.linspace(0.0, Z_SPLINE_MAX, 66)
    halofit = module.Halofit()
    halofit.set_cosmology(
        {
            "Om0": cosmo.Omega_m,
            "Ode0": cosmo.Omega_de,
            "ns": cosmo.n_s,
            "w0": -1.0,
            "wa": 0.0,
            "fnu0": 0.0,
            "sigma8": cosmo.sigma8,
        }
    )
    halofit.set_lgr(
        z_grid,
        np.asarray(growth_D(z_grid, Omega_m=cosmo.Omega_m), dtype=np.float64),
    )
    halofit.set_pklin(k_h.copy(), p_m.copy())

    def b_delta(k1, k2, k3, z):
        shape = np.broadcast(k1, k2, k3).shape
        flat = [
            np.broadcast_to(np.asarray(k, dtype=np.float64), shape).ravel().copy()
            for k in (k1, k2, k3)
        ]
        z_arr = np.full(flat[0].size, float(z), dtype=np.float64)
        out = halofit.get_bihalofit(
            flat[0], flat[1], flat[2], z_arr,
            all_physical=True, which=["Bh3", "Bh1"],
        )
        return np.asarray(out, dtype=np.float64).reshape(shape)

    return b_delta


def make_tree(cosmo: SimCosmo):
    """Tree-level SPT bispectrum on the same linear P(k), for reference."""
    from canoes.cosmo.pk import _CambTablePk

    k_h, p_m = load_pk_wmap3()
    pk = _CambTablePk(k_h, p_m, n_s=cosmo.n_s, extrapolation_low_k=True)

    def f2(ka, kb, cos_ab):
        return (
            5.0 / 7.0
            + 0.5 * cos_ab * (ka / kb + kb / ka)
            + (2.0 / 7.0) * cos_ab**2
        )

    def b_delta(k1, k2, k3, z):
        k1 = np.asarray(k1, float)
        k2 = np.asarray(k2, float)
        k3 = np.asarray(k3, float)
        c12 = (k3**2 - k1**2 - k2**2) / (2.0 * k1 * k2)
        c13 = (k2**2 - k1**2 - k3**2) / (2.0 * k1 * k3)
        c23 = (k1**2 - k2**2 - k3**2) / (2.0 * k2 * k3)
        p1, p2, p3 = pk(k1), pk(k2), pk(k3)
        b0 = 2.0 * (
            f2(k1, k2, c12) * p1 * p2
            + f2(k1, k3, c13) * p1 * p3
            + f2(k2, k3, c23) * p2 * p3
        )
        return b0 * float(growth_D(float(z), Omega_m=cosmo.Omega_m)) ** 4

    return b_delta


# ------------------------------------------------------------- the pipeline
def b_kappa_pipeline(
    ell_triples: NDArray[np.float64],
    *,
    cosmo,
    z_s: float,
    b_delta_fn,
    n_chi: int = 400,
    chi_min: float | None = None,
    k_max: float = 500.0,
) -> tuple[NDArray[np.float64], float]:
    """``B_kappa`` from the canoes vertex factors plus the affine window.

    Returns ``(B_kappa, k_max_reached)``; ``B_kappa`` is dimensionless and
    ``h``-independent (the h-native chi and H0 cancel).
    """
    chi_s = chi_of_z(cosmo, z_s)
    # The Limber map k = l / chi diverges toward the observer, so the inner
    # end of the radial integral is set by the P(k) table rather than by
    # physics. The lensing weight vanishes there as chi^3, so truncating is
    # harmless, but the cut is made explicit and its cost is measured by
    # `--chi-min-scan` rather than left implicit.
    if chi_min is None:
        ell_top = float(np.max(ell_triples))
        chi_min = max(2.0, 1.02 * ell_top / k_max)
    chi = np.exp(np.linspace(np.log(chi_min), np.log(chi_s * (1 - 1e-6)), n_chi))
    _lc, opz = build_lightcone(cosmo, chi, distance_unit="Mpc/h")
    z = opz - 1.0

    poisson_amp = -1.5 * cosmo.Omega_m * _H0_H_PER_MPC**2 * opz
    resp = np.array(
        [
            _kappa3_limber_response_product_h("TTT", float(a), float(o))
            for a, o in zip(poisson_amp, opz)
        ]
    )
    d_lam_d_chi = _lambda_jacobian(opz, convention="project")
    jac, _name, _desc = _kappa3_radial_density_conversion(d_lam_d_chi, "lambda")
    k_window = (chi * (chi_s - chi) / chi_s) / opz**2       # K(lambda)

    # int d lambda = int d chi (d lambda/d chi)
    weight = -d_lam_d_chi * k_window**3 * jac * resp / chi**4

    out = np.empty(ell_triples.shape[0], dtype=np.float64)
    k_max = 0.0
    for i, (l1, l2, l3) in enumerate(ell_triples):
        k1, k2, k3 = l1 / chi, l2 / chi, l3 / chi
        k_max = max(k_max, float(max(k1.max(), k2.max(), k3.max())))
        bd = np.array(
            [
                float(b_delta_fn(np.array([k1[j]]), np.array([k2[j]]),
                                 np.array([k3[j]]), float(z[j]))[0])
                for j in range(chi.size)
            ]
        )
        out[i] = float(np.trapezoid(weight * bd, chi))
    return out, k_max


def b_kappa_textbook(
    ell_triples: NDArray[np.float64],
    *,
    cosmo,
    z_s: float,
    b_delta_fn,
    n_chi: int = 400,
    chi_min: float = 2.0,
) -> NDArray[np.float64]:
    """``int d chi W^3/chi^4 B_delta`` with the textbook W, for the audit trail."""
    chi_s = chi_of_z(cosmo, z_s)
    # The Limber map k = l / chi diverges toward the observer, so the inner
    # end of the radial integral is set by the P(k) table rather than by
    # physics. The lensing weight vanishes there as chi^3, so truncating is
    # harmless, but the cut is made explicit and its cost is measured by
    # `--chi-min-scan` rather than left implicit.
    if chi_min is None:
        ell_top = float(np.max(ell_triples))
        chi_min = max(2.0, 1.02 * ell_top / k_max)
    chi = np.exp(np.linspace(np.log(chi_min), np.log(chi_s * (1 - 1e-6)), n_chi))
    _lc, opz = build_lightcone(cosmo, chi, distance_unit="Mpc/h")
    z = opz - 1.0
    w = (1.5 * cosmo.Omega_m * _H0_H_PER_MPC**2 * opz
         * chi * (chi_s - chi) / chi_s)
    weight = w**3 / chi**4
    out = np.empty(ell_triples.shape[0], dtype=np.float64)
    for i, (l1, l2, l3) in enumerate(ell_triples):
        k1, k2, k3 = l1 / chi, l2 / chi, l3 / chi
        bd = np.array(
            [
                float(b_delta_fn(np.array([k1[j]]), np.array([k2[j]]),
                                 np.array([k3[j]]), float(z[j]))[0])
                for j in range(chi.size)
            ]
        )
        out[i] = float(np.trapezoid(weight * bd, chi))
    return out


# ------------------------------------------------------------------- driver
def main() -> None:
    z_s = 1.0
    chi_s = chi_of_z(SIM, z_s)
    print(f"Sato/Kayo cosmology: Omega_m={SIM.Omega_m}, h={SIM.h}, "
          f"n_s={SIM.n_s}, sigma8={SIM.sigma8}")
    print(f"chi(z_s=1) = {chi_s:.2f} Mpc/h = {chi_s / SIM.h:.1f} Mpc\n")

    b_nl = make_bihalofit(SIM)
    b_tr = make_tree(SIM)

    # Figure 13 bins: dlog10 l = 0.13, running over roughly 100 .. 6000
    ell = 10 ** np.arange(2.0, 3.79, 0.13)

    configs = {
        "equilateral l,l,l": (np.column_stack([ell, ell, ell]), "l^3", ell**3),
        "squeezed l,l,84": (
            np.column_stack([ell, ell, np.full_like(ell, 84.0)]),
            "60 l^2",
            60.0 * ell**2,
        ),
        "squeezed l,l,510": (
            np.column_stack([ell, ell, np.full_like(ell, 510.0)]),
            "300 l^2",
            300.0 * ell**2,
        ),
    }

    results: dict[str, dict] = {}
    for name, (triples, scale_label, scale) in configs.items():
        # squeezed triangles need l >= l_soft/2 to close
        l_soft = float(triples[0, 2])
        keep = triples[:, 0] >= max(l_soft / 2.0, 1.0)
        tri = triples[keep]
        sc = scale[keep]
        ll = ell[keep]
        b_pipe, k_max = b_kappa_pipeline(
            tri, cosmo=SIM, z_s=z_s, b_delta_fn=b_nl
        )
        b_book = b_kappa_textbook(tri, cosmo=SIM, z_s=z_s, b_delta_fn=b_nl)
        b_tree_v, _ = b_kappa_pipeline(
            tri, cosmo=SIM, z_s=z_s, b_delta_fn=b_tr
        )
        # Sensitivity of the answer to the inner radial cut.
        b_tight, _ = b_kappa_pipeline(tri, cosmo=SIM, z_s=z_s, b_delta_fn=b_nl,
                                      chi_min=None, k_max=250.0)
        drift = float(np.max(np.abs(b_tight / b_pipe - 1.0)))
        print(f"--- {name}   (Fig 13 y-axis: {scale_label} B_kappa), "
              f"k_max = {k_max:.1f} h/Mpc; halving the reachable k moves the "
              f"answer by at most {drift:.2e}")
        print(f"{'l':>8} {'B_kappa(BiHalofit)':>20} {'scaled':>12} "
              f"{'B_tree':>13} {'pipe/textbook':>15}")
        for j in range(tri.shape[0]):
            print(f"{ll[j]:8.1f} {b_pipe[j]:20.5e} {sc[j] * b_pipe[j]:12.4e} "
                  f"{b_tree_v[j]:13.4e} {b_pipe[j] / b_book[j]:15.12f}")
        print()
        results[name] = {
            "ell": ll.tolist(),
            "B_kappa_bihalofit": b_pipe.tolist(),
            "B_kappa_tree": b_tree_v.tolist(),
            "scaled_for_fig13": (sc * b_pipe).tolist(),
            "scale_label": scale_label,
            "pipeline_over_textbook": (b_pipe / b_book).tolist(),
            "k_max_h_per_Mpc": k_max,
        }

    out = DATA / "bkappa_squeezed_zs1.json"
    out.write_text(json.dumps(results, indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
