"""Convergence diagnostics of the HIGH-branch quadrature and cutoff.

The deployed vertex truncates the Limber branch at ell_high_max = 1000 with
n_ell = 96 Gauss-Legendre nodes mapped LINEARLY onto [ell_min, ell_high_max].
Two effects are entangled and must be separated before either is believed:

A. quadrature resolution at fixed cutoff (is the deployed 96-node rule
   converged for the integrand it actually integrates?)
B. genuine truncation, i.e. the content above ell = 1000, measured at
   CONSTANT node density (n_ell scaled with the cutoff) so that raising the
   cutoff does not silently coarsen the rule

Naively sweeping the cutoff at fixed n_ell conflates the two: the same 96
nodes spread over a 4x wider interval, so the apparent cutoff dependence
carries a resolution artifact.

Usage (from ``driver_field_emulators/code``)::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python sweep_ell_high_max.py
"""

from __future__ import annotations

import numpy as np

from canoes.sachs.kappa3 import compute_kappa3_sigma3_high

from b_model import FIDUCIAL, load_pk_lin, make_b_delta_fn

ARCMIN = np.pi / (180.0 * 60.0)

#: Production gamma points. gamma = 1000' is deliberately excluded from the
#: ratio tables: the band sum cancels there, so ratios divide by ~0.
GAMMAS = [0.5, 5.45, 42.16, 117.31]
#: Deployed shell grid extremes plus the lensing-kernel peak region.
LAMBDA_MPC = [396.63, 1200.0, 2326.61]


def collapsed_triples(gammas_arcmin: list[float]) -> np.ndarray:
    """The (1, cos gamma, cos gamma) rows the FK 2-point fold queries."""
    cos_gamma = np.cos(np.asarray(gammas_arcmin) * ARCMIN)
    return np.column_stack([np.ones_like(cos_gamma), cos_gamma, cos_gamma])


def zeta(model: str, cutoff: int, n_ell: int, n_phi: int = 64) -> np.ndarray:
    b_delta_fn = None if model == "internal_tree" else make_b_delta_fn(model)
    result = compute_kappa3_sigma3_high(
        collapsed_triples(GAMMAS),
        np.asarray(LAMBDA_MPC) * FIDUCIAL.h,
        pk_matter_today=load_pk_lin(),
        cosmo=FIDUCIAL,
        ell_cut=60,
        ell_high_max=cutoff,
        n_ell=n_ell,
        n_phi=n_phi,
        b_delta_fn=b_delta_fn,
    )
    return np.asarray(result.zeta_TTT)


def _table(title: str, columns: list[str], values: list[np.ndarray]) -> None:
    print(f"\n=== {title} ===")
    print("  gamma[']  lambda[Mpc]  " + "".join(f"{c:>12s}" for c in columns))
    for i, gamma in enumerate(GAMMAS):
        for j, lam in enumerate(LAMBDA_MPC):
            row = "".join(f"{v[i, j]:12.4f}" for v in values)
            print(f"  {gamma:8.2f}  {lam:11.1f}  " + row)


def main() -> None:
    for model in ("internal_tree", "bihalofit"):
        # A. quadrature convergence at the DEPLOYED cutoff
        base = zeta(model, 1000, 96)
        refined = [zeta(model, 1000, n) for n in (192, 384)]
        _table(
            f"{model}: quadrature at cutoff 1000, zeta(n_ell)/zeta(384)",
            ["n_ell=96", "n_ell=192", "n_ell=384"],
            [base / refined[-1], refined[0] / refined[-1],
             refined[-1] / refined[-1]],
        )

        # B. genuine truncation at CONSTANT node density (96 nodes per 1000)
        constant_density = {
            cutoff: zeta(model, cutoff, 96 * cutoff // 1000)
            for cutoff in (1000, 2000, 4000)
        }
        reference = constant_density[4000]
        _table(
            f"{model}: truncation at constant density, zeta(cutoff)/zeta(4000)",
            ["cut=1000", "cut=2000", "cut=4000"],
            [constant_density[c] / reference for c in (1000, 2000, 4000)],
        )


if __name__ == "__main__":
    main()
