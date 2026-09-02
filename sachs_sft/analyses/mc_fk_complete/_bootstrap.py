"""Path and vertex-callable wiring shared by every script in this folder.

The Monte-Carlo engine lives in the paper's analysis tree
(``SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt``) and is imported
read-only; nothing here writes to it.  ``driver_stats`` reaches the three-point
vertex through the module-level handle ``_k3``, so binding the corrected,
permutation-aware callable to the converged table is a single assignment --
the pattern established in ``analyses/revision_2026-08/regenerate_mc.py``.
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
REPO = _HERE.parents[4]
MC_DIR = REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "analyses" / "mc_sachs_2pt"
FIX_DIR = _HERE.parents[2] / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"
PRODUCTS = _HERE.parents[2] / "sftwick_outputs" / "2PCF" / "C_corr_op_K_limber_FK_cut15360_permfix"

TABLE = FIX_DIR / "table_permclosed.npz"
FOLD = PRODUCTS / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"


def wire(table: Path | None = None):
    """Put the MC engine on the path and bind the corrected vertex callable.

    Returns ``(driver_stats, background, sachs_mc_core)``.
    """
    for p in (str(MC_DIR), str(FIX_DIR)):
        if p not in sys.path:
            sys.path.insert(0, p)
    import driver_stats as ds
    import background as bg
    import perm_aware_kappa3_callable as pa

    pa.TABLE_PATH = Path(table or TABLE).resolve()
    pa._CACHE = None
    ds._k3 = pa
    import sachs_mc_core as core
    return ds, bg, core
