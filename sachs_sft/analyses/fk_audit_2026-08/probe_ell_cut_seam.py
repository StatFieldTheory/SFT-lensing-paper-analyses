"""Boundary probe of the ell_cut = 60 dispatch seam.

The vertex splits at ell_cut: LOW is an exact Wigner-3j sum over all
triples with every leg <= 60, HIGH a flat-sky Limber integral over triples
whose largest leg exceeds it. The two are different approximations meeting
at a hard boundary, and after the rebuild they also carry DIFFERENT
bispectra, because LOW's FFTlog machinery needs a separable (tree) model
while HIGH takes any callable. Two questions follow:

* dispatch error: do the two methods agree on a shared band, with the same
  bispectrum? This is the classic boundary check, run first as a control.
* model mismatch: how much does keeping tree in LOW while HIGH goes
  nonlinear cost at the seam?

Method (per the project's boundary-validation rule, the dispatcher is
bypassed): compare a LOW band-difference, zeta(ell_max=60) minus
zeta(ell_max=ell0), against a HIGH window integral run over the same band,
(ell0, 60].

Usage (from ``sachs_sft/analyses/fk_audit_2026-08``)::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python probe_ell_cut_seam.py
"""

from __future__ import annotations

import numpy as np

from canoes.sachs.kappa3 import (
    compute_kappa3_sigma3_high,
    compute_kappa3_zeta_table,
)

from b_model import FIDUCIAL, load_pk_lin, make_b_delta_fn

ARCMIN = np.pi / (180.0 * 60.0)
#: Band edges: the probe compares the methods on (ELL0, ELL_CUT].
ELL0, ELL_CUT = 40, 60
#: Collapsed family points spanning the FK range, plus an open triangle.
GAMMAS = [0.5, 42.16, 1000.0]
#: Extreme shells first: the boundary rule requires the corners.
LAMBDA_MPC = [396.63, 1200.0, 2326.61]


def triples() -> np.ndarray:
    cos_gamma = np.cos(np.asarray(GAMMAS) * ARCMIN)
    collapsed = np.column_stack([np.ones_like(cos_gamma), cos_gamma, cos_gamma])
    return np.vstack([collapsed, np.array([[0.8, 0.7, 0.6]])])


def low_band(ell0: int, ell_cut: int) -> np.ndarray:
    """zeta from the exact branch restricted to the band (ell0, ell_cut]."""
    common = dict(
        pk_matter_today=load_pk_lin(),
        cosmo=FIDUCIAL,
        Nmax=128,
        spt_kind="tree_phi",
        units="physical",
        lambda_convention="project",
        radial_discretization="sample",
        radial_measure="lambda",
        parallelism="sequential",
        n_workers=1,
    )
    lam = np.asarray(LAMBDA_MPC) * FIDUCIAL.h
    upper = compute_kappa3_zeta_table(triples(), lam, ell_max=ell_cut, **common)
    lower = compute_kappa3_zeta_table(triples(), lam, ell_max=ell0, **common)
    return np.asarray(upper.zeta_TTT) - np.asarray(lower.zeta_TTT)


def high_band(ell0: int, ell_cut: int, model: str) -> np.ndarray:
    """zeta from the Limber branch windowed onto the same band."""
    b_delta_fn = None if model == "internal_tree" else make_b_delta_fn(model)
    result = compute_kappa3_sigma3_high(
        triples(),
        np.asarray(LAMBDA_MPC) * FIDUCIAL.h,
        pk_matter_today=load_pk_lin(),
        cosmo=FIDUCIAL,
        ell_cut=ell0,
        ell_high_max=ell_cut,
        n_ell=96,
        n_phi=64,
        units="physical",
        lambda_convention="project",
        radial_discretization="sample",
        radial_measure="lambda",
        b_delta_fn=b_delta_fn,
    )
    return np.asarray(result.zeta_TTT)


def report(title: str, reference: np.ndarray, other: np.ndarray) -> None:
    print(f"\n=== {title} ===")
    print("  triple            lambda[Mpc]      exact-3j        Limber      rel.err")
    labels = [f"gamma={g:.2f}'" for g in GAMMAS] + ["open(0.8,0.7,0.6)"]
    for i, label in enumerate(labels):
        for j, lam in enumerate(LAMBDA_MPC):
            a, b = reference[i, j], other[i, j]
            rel = abs(b - a) / max(abs(a), abs(b), 1e-300)
            flag = "OK  " if rel < 1e-1 else "BAD "
            print(f"  {label:17s} {lam:9.1f}  {a:13.4e} {b:13.4e}  {rel:9.3f} {flag}")


def main() -> None:
    exact = low_band(ELL0, ELL_CUT)
    print(f"band = ({ELL0}, {ELL_CUT}]  (the dispatcher is bypassed on both sides)")
    report(
        "control: same tree bispectrum on both sides (dispatch error only)",
        exact,
        high_band(ELL0, ELL_CUT, "internal_tree"),
    )
    report(
        "model mismatch: tree in LOW vs BiHalofit in HIGH",
        exact,
        high_band(ELL0, ELL_CUT, "bihalofit"),
    )
    print(
        "\nThe control measures how well the flat-sky Limber branch reproduces\n"
        "the exact curved-sky sum at the seam; the second table adds the cost\n"
        "of leaving LOW at tree level once HIGH goes nonlinear."
    )


if __name__ == "__main__":
    main()
