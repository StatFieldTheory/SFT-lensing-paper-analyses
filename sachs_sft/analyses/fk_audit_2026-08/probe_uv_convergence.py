"""Measure the UV convergence degree of the collapsed kappa3 vertex.

Sweeping ell_high_max answers "how much have I accumulated"; it cannot say
whether the sum converges, because a slowly growing partial sum and a
divergent one look alike over one decade. The differential contribution
can: evaluate the HIGH branch on successive octave bands (L, 2L] and read
off how the band contribution scales with L.

    band(L) ~ L^(-p),  p > 0  -> the sum converges, ell_high_max is a
                                 numerical parameter to be converged
                 p <= 0  -> the sum diverges, ell_high_max is a physical
                            smoothing scale that must be chosen and stated

Each band is evaluated with its own quadrature over that band alone, so
node density is uniform across bands by construction and the comparison
is not contaminated by the node-dilution effect that a fixed-n_ell cutoff
sweep suffers from.

Analytic expectation, for orientation. The collapsed vertex integrates the
squeezed bispectrum over the hard pair's transverse wavevector, so the band
contribution goes as INT_L^2L dl l B(l/chi, l/chi, l_soft/chi), and with
B_sq ~ P(k_soft) P(k_hard) that is INT dl l P(l/chi). A linear CDM spectrum
falls as k^-3 at high k, giving band ~ L^-1, convergent. A nonlinear
spectrum falls as roughly k^-1.5, giving band ~ L^+0.5, divergent. The
measurement below tests exactly this.

Usage (from ``sachs_sft/analyses/fk_audit_2026-08``)::

    CAN=/Users/zzhang/projects/angular_statistics/canoes
    PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python probe_uv_convergence.py
"""

from __future__ import annotations

import numpy as np

from canoes.sachs.kappa3 import compute_kappa3_sigma3_high

from b_model import FIDUCIAL, load_pk_lin, make_b_delta_fn

ARCMIN = np.pi / (180.0 * 60.0)

#: Octave band edges. Each band is integrated on its own quadrature.
#: Capped at 16000 because the third leg reaches 2 * ell_max, and the
#: fiducial CAMB table stops at k_max = 206 h/Mpc: at the nearest shell
#: (chi = 266 Mpc/h) that is ell = 55000, so bands beyond 16000 would ask
#: the tree branch to extrapolate its own P(k) table.
BAND_EDGES = [125, 250, 500, 1000, 2000, 4000, 8000, 16000]
#: Small separations, where the collapsed vertex is largest and the FK
#: amplitude the draft quotes is measured.
GAMMAS = [0.5, 5.45, 42.16]
#: Near, middle and far shell.
LAMBDA_MPC = [396.63, 1200.0, 2326.61]
#: Nodes per band. The band is one octave wide regardless of L, so a fixed
#: count is a fixed density in ln ell.
N_ELL, N_PHI = 96, 64


def triples() -> np.ndarray:
    cos_gamma = np.cos(np.asarray(GAMMAS) * ARCMIN)
    return np.column_stack([np.ones_like(cos_gamma), cos_gamma, cos_gamma])


def band_value(low: int, high: int, model: str) -> np.ndarray:
    b_delta_fn = None if model == "internal_tree" else make_b_delta_fn(model)
    result = compute_kappa3_sigma3_high(
        triples(),
        np.asarray(LAMBDA_MPC) * FIDUCIAL.h,
        pk_matter_today=load_pk_lin(),
        cosmo=FIDUCIAL,
        ell_cut=low,
        ell_high_max=high,
        n_ell=N_ELL,
        n_phi=N_PHI,
        units="physical",
        lambda_convention="project",
        radial_discretization="sample",
        radial_measure="lambda",
        b_delta_fn=b_delta_fn,
    )
    return np.asarray(result.zeta_TTT)


def main() -> None:
    for model in ("internal_tree", "bihalofit"):
        bands = {}
        for low, high in zip(BAND_EDGES[:-1], BAND_EDGES[1:]):
            bands[(low, high)] = band_value(low, high, model)

        print(f"\n=== {model}: octave-band contributions to zeta_TTT ===")
        for i, gamma in enumerate(GAMMAS):
            for j, lam in enumerate(LAMBDA_MPC):
                row = [bands[b][i, j] for b in bands]
                print(f"\n  gamma = {gamma}', lambda = {lam} Mpc")
                print("    band            value        |value| ratio to previous"
                      "     local slope p")
                previous = None
                for (low, high), value in zip(bands, row):
                    if previous is None or previous == 0.0:
                        print(f"    ({low:6d},{high:6d}]  {value:+.4e}"
                              f"          --                    --")
                    else:
                        ratio = abs(value) / abs(previous)
                        slope = -np.log2(ratio) if ratio > 0 else np.nan
                        print(f"    ({low:6d},{high:6d}]  {value:+.4e}"
                              f"     {ratio:9.4f}            {slope:+8.3f}")
                    previous = value

        print(f"\n  Reading: band ~ L^(-p). p > 0 means the tail dies and the"
              f" sum converges;\n  p <= 0 means it does not.")


if __name__ == "__main__":
    main()
