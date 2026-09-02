"""Coordinate-map figure for the cosmology section: chi(lambda) and z(lambda).

Compact two-panel figure that visualises the one-to-one
``(lambda, n_hat) <-> (z, r_hat)`` map of the manuscript. Both panels share
the affine-parameter lambda axis on x; the left panel renders chi(lambda),
the right panel z(lambda). The cosmology is the project-wide fiducial
(Omega_m = 0.3160919980475834, h = 0.6711), so the curves coincide with the
radial grid the analysis pipeline consumes.

Provenance. Promoted on 2026-09-02 from
``_archive/canoes_pipeline/scripts/plot_coord_maps.py`` (2026-05-18), with the
two background helpers of the archived ``canoes_pipeline/cosmo.py`` and the
canoes default constants vendored below, so the figure needs neither the
canoes package nor the pre-2026-06-17 package layout. The rendering code is
unchanged; the output reproduces the deployed ``figures/cosmology_coord_maps.pdf``
byte for byte apart from the PDF creation date (checked with the PyCCL env,
matplotlib 3.10.8).

Output (default): ``figures/outputs/cosmology_coord_maps.{pdf,png}`` next to
this script. The paper copy under the repository ``figures/`` is deployed by
the reproduction driver, not by this script; pass ``--out-stem`` to redirect.

Run:
    /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python \
        SFT-lensing-paper-analyses/figures/make_cosmology_coord_maps.py
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import cumulative_trapezoid  # noqa: E402

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from _plot_style import PALETTE, apply_rcparams  # noqa: E402

apply_rcparams()

# canoes.cosmo.background defaults (canoes d7b20ce), vendored verbatim so the
# arithmetic, and therefore the double-precision values, match the archived
# generator exactly.
C_KM_S: float = 299792.458
H_DEFAULT: float = 0.6711
OMEGA_C_H2_DEFAULT: float = 0.12029
OMEGA_B_H2_DEFAULT: float = 0.02207
OMEGA_M_DEFAULT: float = (OMEGA_C_H2_DEFAULT + OMEGA_B_H2_DEFAULT) / H_DEFAULT**2
HUBBLE_RADIUS_MPC_H: float = C_KM_S / 100.0  # c / (100 km/s/Mpc) in Mpc/h

Z_MAX: float = 5.0
N_GRID: int = 401
COLOR_CHI = PALETTE[0]  # blue
COLOR_Z = PALETTE[1]    # vermilion


@dataclass(frozen=True)
class FiducialCosmology:
    """Minimal fiducial cosmology (Omega_m, h), as in the archived pipeline."""

    Omega_m: float = OMEGA_M_DEFAULT
    h: float = H_DEFAULT
    n_s: float = 0.97


def _chi_h_from_z(
    z: np.ndarray, *, Omega_m: float, n_grid: int = 5001, z_max: float = 6.0
) -> np.ndarray:
    """Cumulative-trapezoid comoving distance in Mpc/h (flat LCDM, no radiation)."""
    z_grid = np.linspace(0.0, max(z_max, float(np.max(z)) + 0.1), n_grid)
    E = np.sqrt(Omega_m * (1.0 + z_grid) ** 3 + (1.0 - Omega_m))
    dchi_dz = HUBBLE_RADIUS_MPC_H / E
    chi_grid = np.zeros_like(z_grid)
    chi_grid[1:] = cumulative_trapezoid(dchi_dz, z_grid)
    return np.interp(z, z_grid, chi_grid)


def project_lambda_h_from_z(
    z: np.ndarray, *, Omega_m: float, n_grid: int = 5001, z_max: float = 6.0
) -> np.ndarray:
    """Photon affine parameter (E_0 = 1 convention) in Mpc/h.

    Project convention: dlambda/dchi = 1/(1+z)^2 = a^2(chi), so
    dlambda/dz = HUBBLE_RADIUS / ((1+z)^2 E(z)).
    """
    z_grid = np.linspace(0.0, max(z_max, float(np.max(z)) + 0.1), n_grid)
    E = np.sqrt(Omega_m * (1.0 + z_grid) ** 3 + (1.0 - Omega_m))
    dlam_dz = HUBBLE_RADIUS_MPC_H / ((1.0 + z_grid) ** 2 * E)
    lam_grid = np.zeros_like(z_grid)
    lam_grid[1:] = cumulative_trapezoid(dlam_dz, z_grid)
    return np.interp(z, z_grid, lam_grid)


def _build_curves(
    cosmo: FiducialCosmology, z_max: float, n: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(lambda_Mpc, chi_Mpc, z)`` on a fine grid in physical Mpc."""
    z_grid = np.linspace(0.0, z_max, n)
    chi_h = _chi_h_from_z(z_grid, Omega_m=cosmo.Omega_m)
    lam_h = project_lambda_h_from_z(z_grid, Omega_m=cosmo.Omega_m)
    chi_Mpc = chi_h / cosmo.h
    lam_Mpc = lam_h / cosmo.h
    return lam_Mpc, chi_Mpc, z_grid


def render(out_stem: Path) -> None:
    cosmo = FiducialCosmology()
    lam, chi, z = _build_curves(cosmo, Z_MAX, N_GRID)

    # Bump font sizes substantially for this paper-figure insertion
    # (overrides the _plot_style defaults locally; restored after savefig).
    panel_rc = {
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
    }
    saved = {k: plt.rcParams[k] for k in panel_rc}
    plt.rcParams.update(panel_rc)

    # Wider + slightly taller frame so the bigger labels have room to breathe.
    fig, (ax_chi, ax_z) = plt.subplots(
        1, 2, figsize=(11.6, 4.6), sharex=True, constrained_layout=True
    )

    # --- left: chi(lambda) ---
    ax_chi.plot(lam, chi, color=COLOR_CHI, lw=2.0)
    ax_chi.set_xlabel(r"$\lambda\;[\mathrm{Mpc}]$")
    ax_chi.set_ylabel(r"$\chi(\lambda)\;[\mathrm{Mpc}]$")
    ax_chi.set_xlim(0.0, lam.max())
    ax_chi.set_ylim(0.0, chi.max() * 1.04)

    # --- right: z(lambda) ---
    ax_z.plot(lam, z, color=COLOR_Z, lw=2.0)
    ax_z.set_xlabel(r"$\lambda\;[\mathrm{Mpc}]$")
    ax_z.set_ylabel(r"$z(\lambda)$")
    ax_z.set_xlim(0.0, lam.max())
    ax_z.set_ylim(0.0, z.max() * 1.04)

    # Reference markers at z in {0.5, 1, 2, 3, 5} for orientation.
    # Labels are placed to the LEFT of each marker on the chi panel only
    # (the z panel reads its own value off the y-axis).
    for z_ref in (0.5, 1.0, 2.0, 3.0, 5.0):
        if z_ref > z.max():
            continue
        i = int(np.argmin(np.abs(z - z_ref)))
        ax_chi.plot(lam[i], chi[i], marker="o", ms=6.5,
                    color=COLOR_CHI, mfc="white", mew=1.5, zorder=4)
        ax_z.plot(lam[i], z[i], marker="o", ms=6.5,
                  color=COLOR_Z, mfc="white", mew=1.5, zorder=4)
        ax_chi.annotate(
            f"$z={z_ref:g}$",
            xy=(lam[i], chi[i]),
            xytext=(-9, 0), textcoords="offset points",
            ha="right", va="center",
            fontsize=16, color="0.20",
        )

    pdf = out_stem.with_suffix(".pdf")
    png = out_stem.with_suffix(".png")
    fig.savefig(pdf)
    fig.savefig(png)
    plt.close(fig)
    plt.rcParams.update(saved)
    print(f"wrote {pdf}\n      {png}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_out = HERE / "outputs" / "cosmology_coord_maps"
    parser.add_argument("--out-stem", type=Path, default=default_out)
    args = parser.parse_args()
    args.out_stem.parent.mkdir(parents=True, exist_ok=True)
    render(args.out_stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
