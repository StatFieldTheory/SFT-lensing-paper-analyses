"""Vendored minimal cosmology + matter-P(k) loader for the equal_time_limber κ³.

Self-contained copy of the two helpers the build needs (trimmed from the
archived ``canoes_pipeline/cosmo.py`` + ``pk.py``).  Kept HERE so this callable
folder has no cross-folder import dependency.

The cosmology is the canoes default == the STF fiducial (verified):
``OMEGA_M_DEFAULT = 0.3160919980475834``, ``H_DEFAULT = 0.6711`` — identical to
the corr_op / R_contracted callables, so all four sachs_sft callables share one
cosmology.

P(k) note: the equal_time_limber build feeds the MATTER ``P_δ(k, z=0)`` (NOT
``P_φ``); ``spt_kind="tree_phi"`` applies the per-leg Poisson factor A(a)
internally.  The +3.08e-5 FK baseline was measured with ``PCAMBz0.txt``, so this
loader uses that table (documented difference vs the corr_op/R_contracted
``PCAMB_pyccl_stf_fid_z0.txt``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from canoes.cosmo.background import H_DEFAULT, OMEGA_M_DEFAULT
from canoes.cosmo.pk import _CambTablePk

# The +3.08e-5 baseline P(k) table (canoes examples).
def _canoes_examples_data() -> Path:
    """canoes examples/data directory, following the package rather than a
    hardcoded checkout path (the repository has been moved twice)."""
    import os

    env = os.environ.get("CANOES_ROOT")
    if env:
        return Path(env) / "examples" / "data"
    import canoes

    return Path(canoes.__file__).resolve().parents[2] / "examples" / "data"


_DEFAULT_CAMB_TABLE = _canoes_examples_data() / "PCAMBz0.txt"


@dataclass(frozen=True)
class FiducialCosmology:
    """CosmoLike (Omega_m, h) + n_s — canoes default == STF fiducial."""

    Omega_m: float = OMEGA_M_DEFAULT
    h: float = H_DEFAULT
    n_s: float = 0.97


def load_pk_delta(
    n_s: float = 0.97, table_path: Path = _DEFAULT_CAMB_TABLE,
) -> _CambTablePk:
    """Load raw CAMB matter ``P_δ(k, z=0)`` [h-units] for the tree_phi κ³ build.

    No Poisson conversion here: ``compute_kappa3_zeta_table(spt_kind="tree_phi")``
    applies the per-leg Poisson amplitude A(a) internally.
    """
    arr = np.loadtxt(table_path)
    k_h, p_m = arr[:, 0], arr[:, 1]
    return _CambTablePk(k_h.copy(), p_m.copy(), n_s=float(n_s),
                        extrapolation_low_k=True)
