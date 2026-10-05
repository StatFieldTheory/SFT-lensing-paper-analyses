"""Selected C23 input with finite-amplitude and sampled-sign acceptance.

The package-local table uses ridge interpolation, output_h_power=0 and
ell_lo_threshold=21. This binding grants no global covariance accuracy.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

def _resolve_canoes_root() -> Path:
    """Use an explicit checkout or the installed canoes package."""
    import os
    env = os.environ.get("CANOES_ROOT")
    if env:
        return Path(env).resolve()
    try:
        import canoes
        return Path(canoes.__file__).resolve().parents[2]
    except ImportError as error:
        raise RuntimeError("Install canoes or set CANOES_ROOT") from error


_CANOES_ROOT = _resolve_canoes_root()
for _p in (_CANOES_ROOT / "src", _CANOES_ROOT):
    if _p.is_dir() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    CorrOpTable2D,
    load_table_2d_from_npz,
    make_C_fn_lambda_project,
)

TABLE_PATH = Path(__file__).resolve().parent / 'inputs/corr_op_table_zs5_23pt.npz'
TABLE_SHA256 = '5b371f8632fd4c24c829a23d5ae0d946223c2fbddcba63b8c7c42a109037aba2'
INTERP_METHOD = "ridge"

# --- sft-wick contract (the only values the config / sft-wick needs) ----------
QUERY_COORDINATE = "lambda_project_mpc"   # t1, t2 are λ_project in physical Mpc
QUERY_UNITS = "physical Mpc"
COMPONENT_ORDER = ("Phi_00", "Re Psi_0", "Im Psi_0")
OUTPUT_H_POWER = 0.0                       # no h-rescale (value is self-consistent)


def _load_table() -> CorrOpTable2D:
    digest = hashlib.sha256(TABLE_PATH.read_bytes()).hexdigest()
    if digest != TABLE_SHA256:
        raise ValueError(f"Candidate table SHA256 mismatch: {digest}, expected {TABLE_SHA256}")
    table = load_table_2d_from_npz(str(TABLE_PATH))
    # Fail-loud: corr_op's natural coordinate is λ_project; reject a χ-indexed
    # table (the chi-vs-lambda trap). Read from the TABLE config, never hand-typed.
    src = str(table.config.source_coord)
    if src != "lambda_project_mpc":
        raise ValueError(
            f"corr_op table source_coord is {src!r}, expected 'lambda_project_mpc'. "
            f"This callable serves the λ_project (affine-source) corr_op propagator; "
            f"a χ-indexed table is the wrong product here."
        )
    if table.config.cosmology is None:
        raise ValueError(
            "corr_op table is missing config.cosmology (needed for the "
            "λ_project → χ conversion in make_C_fn_lambda_project)."
        )
    if table.config.ell_lo_threshold != 21:
        raise ValueError("Candidate requires ell_lo_threshold=21")
    return table


_TABLE = _load_table()
_C_fn_impl = make_C_fn_lambda_project(
    _TABLE, interp_method=INTERP_METHOD, output_h_power=OUTPUT_H_POWER
)
if not hasattr(_C_fn_impl, "batch"):
    raise RuntimeError(
        "make_C_fn_lambda_project did not expose a `.batch` path; the vectorized "
        "C_fn_batch (sft-wick c_closed_form_vectorized=True) cannot be served."
    )


def C_fn(
    n1: NDArray[np.float64], t1: float,
    n2: NDArray[np.float64], t2: float,
) -> NDArray[np.float64]:
    """Correlation propagator ⟨Φ Φ⟩ as a (3, 3) tensor (scalar path).

    Parameters
    ----------
    n1, n2 : unit direction vectors (3,).
    t1, t2 : source positions in λ_project [physical Mpc].
    """
    return np.asarray(_C_fn_impl(n1, t1, n2, t2), dtype=np.float64)


def C_fn_batch(
    n1: NDArray[np.float64], t1: NDArray[np.float64],
    n2: NDArray[np.float64], t2: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Vectorized ⟨Φ Φ⟩ over a batch of samples → (n, 3, 3).

    This is the sft-wick **vectorized-C contract** (``c_closed_form_vectorized=
    True``): ``C_fn_batch(x1, t1, x2, t2) -> (n, 3, 3)`` where ``t1, t2`` are
    ``(n,)`` arrays of λ_project [physical Mpc] (a scalar is accepted and yields
    ``(1, 3, 3)``, which sft-wick squeezes), and ``x1, x2`` are direction vectors
    ``(n, 3)`` or a single ``(3,)``.

    Delegates to the candidate table grouped-batch implementation.
    """
    t1a = np.atleast_1d(np.asarray(t1, dtype=np.float64))
    t2a = np.atleast_1d(np.asarray(t2, dtype=np.float64))
    return np.asarray(_C_fn_impl.batch(n1, t1a, n2, t2a), dtype=np.float64)


# Convention attributes (read from the table; for the config/manifest to record).
C_fn.query_coordinate = QUERY_COORDINATE       # type: ignore[attr-defined]
C_fn.query_units = QUERY_UNITS                 # type: ignore[attr-defined]
C_fn.component_order = COMPONENT_ORDER          # type: ignore[attr-defined]
C_fn.source_coord = str(_TABLE.config.source_coord)  # type: ignore[attr-defined]
C_fn.cosmology = _TABLE.config.cosmology        # type: ignore[attr-defined]
C_fn.ell_hi_threshold = _TABLE.config.ell_hi_threshold  # type: ignore[attr-defined]
C_fn.output_h_power = OUTPUT_H_POWER            # type: ignore[attr-defined]
C_fn.table_path = str(TABLE_PATH)              # type: ignore[attr-defined]

C_fn.interp_method = INTERP_METHOD
C_fn.table_hash = TABLE_SHA256
C_fn.table_sha256 = TABLE_SHA256
C_fn.ell_lo_threshold = _TABLE.config.ell_lo_threshold

# Mirror the convention attributes onto the vectorized variant + flag it.
C_fn_batch.vectorized = True                    # type: ignore[attr-defined]
for _a in ("query_coordinate", "query_units", "component_order", "source_coord",
           "cosmology", "ell_hi_threshold", "output_h_power", "table_path",
           "interp_method", "table_hash", "table_sha256", "ell_lo_threshold"):
    setattr(C_fn_batch, _a, getattr(C_fn, _a))

__all__ = ["C_fn", "C_fn_batch", "TABLE_PATH", "QUERY_COORDINATE", "COMPONENT_ORDER"]
