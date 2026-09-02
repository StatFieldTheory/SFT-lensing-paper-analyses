"""Fiducial cosmology, linear P(k), growth and sigma8 for the B_model chain.

Everything here mirrors what the deployed equal_time_limber kappa3 build
actually uses, so that a tree-level B_model call reproduces the canoes
internal tree assembly bit-for-bit:

* cosmology  : canoes defaults (Omega_m, h) plus n_s = 0.97
* P(k)       : the CAMB matter table the 3-point build reads, in h-units
               (k in h/Mpc, P in (Mpc/h)^3), interpolated by the same
               ``canoes.cosmo.pk._CambTablePk`` the build uses
* growth     : ``canoes.cosmo.growth.growth_D``, D(0) = 1

Note the documented split of P(k) tables between the two- and three-point
sides of the pipeline: the 3-point vertex is built from ``PCAMBz0.txt``
(sigma8 = 0.808988), the corr_op 2-point side from
``PCAMB_pyccl_stf_fid_z0.txt`` (sigma8 = 0.810000). This module defaults to
the 3-point table because that is what the vertex rebuild must match.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# canoes is the source of truth for growth and P(k) interpolation.
from canoes.cosmo.background import H_DEFAULT, OMEGA_M_DEFAULT
from canoes.cosmo.growth import growth_D
from canoes.cosmo.pk import _CambTablePk

CANOES_DATA = Path(__file__).resolve().parents[4]
_CANOES_ROOT_ENV = "CANOES_ROOT"


def _default_data_dir() -> Path:
    """Locate the canoes examples/data directory holding the CAMB tables."""
    import os

    root = os.environ.get(_CANOES_ROOT_ENV)
    if root:
        return Path(root) / "examples" / "data"
    import canoes

    # src layout: .../canoes/src/canoes/__init__.py -> repo root is parents[2]
    return Path(canoes.__file__).resolve().parents[2] / "examples" / "data"


#: P(k) table used by the deployed equal_time_limber kappa3 build.
PK_TABLE_3PT = "PCAMBz0.txt"
#: P(k) table used by the corr_op two-point side (for cross-checks only).
PK_TABLE_2PT = "PCAMB_pyccl_stf_fid_z0.txt"

#: sigma8 of each table, from direct top-hat integration (see tests).
SIGMA8_3PT = 0.808988
SIGMA8_2PT = 0.810000


@dataclass(frozen=True)
class Fiducial:
    """Fiducial cosmology of the STF lensing pipeline (flat LCDM, c = 1)."""

    Omega_m: float = OMEGA_M_DEFAULT
    h: float = H_DEFAULT
    n_s: float = 0.97
    #: Not pinned by canoes; needed only by baryon models. Planck-like value.
    Omega_b: float = 0.049004
    #: sigma8 of the 3-point P(k) table (the one the vertex is built from).
    sigma8: float = SIGMA8_3PT

    @property
    def Omega_de(self) -> float:
        return 1.0 - self.Omega_m


FIDUCIAL = Fiducial()


def load_pk_table(
    name: str = PK_TABLE_3PT, data_dir: Path | None = None
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return the raw ``(k [h/Mpc], P [(Mpc/h)^3])`` columns of a CAMB table."""
    directory = data_dir if data_dir is not None else _default_data_dir()
    arr = np.loadtxt(Path(directory) / name)
    return arr[:, 0].copy(), arr[:, 1].copy()


def load_pk_lin(
    name: str = PK_TABLE_3PT,
    n_s: float = FIDUCIAL.n_s,
    data_dir: Path | None = None,
) -> _CambTablePk:
    """Linear matter P(k, z=0) callable, interpolated exactly as canoes does."""
    k_h, p_m = load_pk_table(name, data_dir)
    return _CambTablePk(k_h, p_m, n_s=float(n_s), extrapolation_low_k=True)


def sigma8_of_table(
    name: str = PK_TABLE_3PT, data_dir: Path | None = None
) -> float:
    """sigma8 by direct top-hat integration of a tabulated linear P(k)."""
    k, p = load_pk_table(name, data_dir)
    x = 8.0 * k
    window = 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
    integrand = k**3 * p * window**2 / (2.0 * np.pi**2)
    return float(np.sqrt(np.trapezoid(integrand, np.log(k))))


def growth(z: float | NDArray[np.float64], cosmo: Fiducial = FIDUCIAL):
    """Linear growth factor D(z) normalised to D(0) = 1 (canoes convention)."""
    return growth_D(z, Omega_m=cosmo.Omega_m)
