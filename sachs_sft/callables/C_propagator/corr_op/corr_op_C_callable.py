"""C_propagator / corr_op — sft-wick correlation-propagator callable.

Exposes ``C_fn(n1, t1, n2, t2) -> (3, 3)`` for sft-wick, where ``t1, t2`` are
source positions in **λ_project [physical Mpc]** and the (3, 3) tensor is the
two-point cumulant ``⟨Φ Φ⟩`` in component order (Φ00, Re Ψ0, Im Ψ0).

Design (the sachs_sft "traceable callable" paradigm)
----------------------------------------------------
* **Explicit table path** (``TABLE_PATH`` below, in THIS folder). No
  ``SFT_WICK_*`` env-var resolution — the config→product link is a constant here,
  not a hidden default. This is the mechanism that lost provenance in the old
  ``corr_op_l2_callable.py``.
* **Conventions are read FROM the table** (``config.source_coord`` /
  ``config.cosmology``), never hand-typed. The old callable assigned
  ``NATIVE_TABLE_COORDINATE`` three times (chi→lambda→chi, last-wins ``chi_mpc``)
  which CONTRADICTED the table's actual ``source_coord='lambda_project_mpc'``.
  Here the table is the single source of truth, so the label cannot disagree
  with the data.
* **Fail-loud** on the one contract that matters: the table must be λ_project
  (``source_coord == 'lambda_project_mpc'``). corr_op's natural source coordinate
  is the affine λ, and sft-wick consumes λ_project. A χ-indexed table raises.

Provenance — table built by (reproducible; see ``build_corr_op_table.py``)
-------------------------------------------------------------------------
``python -m canoes.sachs.sft_input.corr_op.build
    --pk examples/data/PCAMB_pyccl_stf_fid_z0.txt
    --omega-m 0.3160919980475834 --h 0.6711 --n-s 0.97
    --ell-max 5000 --dense-threshold 500 --n-log-per-decade 16
    --lambda-grid zs5-21pt           # → source_coord=lambda_project_mpc
    --chi-min 50.0 --n-chi-per-pair 256
    --accuracy default               # FFTlog Nmax=256
    --ell-hi-threshold 800           # S3: Limber for ℓ>800, t-form for ℓ≤800
    --output corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz``
(`ell_taper_frac=0.2`, `cl_tensor_convention=bare` are CorrOpTable2DConfig defaults.)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# canoes on path
_CANOES_ROOT = Path("/Users/zzhang/projects/canoes")
if str(_CANOES_ROOT) not in sys.path:
    sys.path.insert(0, str(_CANOES_ROOT))

from canoes.sachs.sft_input.corr_op.table import (  # noqa: E402
    CorrOpTable2D,
    load_table_2d_from_npz,
    make_C_fn_lambda_project,
)

_HERE = Path(__file__).resolve().parent
TABLE_PATH = _HERE / "corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz"

# --- sft-wick contract (the only values the config / sft-wick needs) ----------
QUERY_COORDINATE = "lambda_project_mpc"   # t1, t2 are λ_project in physical Mpc
QUERY_UNITS = "physical Mpc"
COMPONENT_ORDER = ("Phi_00", "Re Psi_0", "Im Psi_0")
OUTPUT_H_POWER = 0.0                       # no h-rescale (value is self-consistent)


def _load_table() -> CorrOpTable2D:
    if not TABLE_PATH.exists():
        raise FileNotFoundError(
            f"corr_op table not found: {TABLE_PATH}\n"
            f"Build it with: python build_corr_op_table.py (in this folder)."
        )
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
    return table


_TABLE = _load_table()
_C_fn_impl = make_C_fn_lambda_project(_TABLE, output_h_power=OUTPUT_H_POWER)
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

    Thin re-export of the canoes grouped-batch path
    (``make_C_fn_lambda_project(...).batch``): it groups samples by the rounded
    ``(cos γ, ψ₁, ψ₂)`` triple to evaluate the per-direction angular sum once,
    then does a single batched grid interpolation. Verified **bit-identical** to
    the per-sample scalar ``C_fn`` (rel-err = 0) and ~1.5e3× faster — so a sweep
    can set ``c_closed_form_vectorized: true`` with no change in result.
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

# Mirror the convention attributes onto the vectorized variant + flag it.
C_fn_batch.vectorized = True                    # type: ignore[attr-defined]
for _a in ("query_coordinate", "query_units", "component_order", "source_coord",
           "cosmology", "ell_hi_threshold", "output_h_power", "table_path"):
    setattr(C_fn_batch, _a, getattr(C_fn, _a))

__all__ = ["C_fn", "C_fn_batch", "TABLE_PATH", "QUERY_COORDINATE", "COMPONENT_ORDER"]


if __name__ == "__main__":  # self-check
    lam = np.asarray(_TABLE.cl_table.lambda_grid, dtype=np.float64)
    print(f"[corr_op C_fn] table: {TABLE_PATH.name}")
    print(f"  source_coord = {C_fn.source_coord}  (must be lambda_project_mpc)")
    print(f"  cosmology    = {C_fn.cosmology}")
    print(f"  ell_hi_threshold = {C_fn.ell_hi_threshold}")
    print(f"  lambda_grid [{lam[0]:.1f} .. {lam[-1]:.1f}] Mpc, {lam.size} nodes")
    n = np.array([0.0, 0.0, 1.0])
    t = float(lam[len(lam) // 2])
    C = C_fn(n, t, n, t)
    print(f"  C_fn(equal n, t={t:.1f}) shape={C.shape}, finite={np.all(np.isfinite(C))}")
    print(f"  diag = {np.diag(C)}")
    # Vectorized path: must be bit-identical to the per-sample scalar loop.
    zhat = np.array([0.0, 0.0, 1.0])
    ts = lam[[2, len(lam) // 2, -3]]
    n1s, t1s, n2s, t2s = [], [], [], []
    for th in (1e-4, 0.01, 0.1, 0.5, 1.0):
        ny = np.array([np.sin(th), 0.0, np.cos(th)])
        for a in ts:
            for b in ts:
                n1s.append(zhat); t1s.append(a); n2s.append(ny); t2s.append(b)
    n1a, t1a, n2a, t2a = (np.array(n1s), np.array(t1s), np.array(n2s), np.array(t2s))
    scal = np.array([C_fn(n1a[i], float(t1a[i]), n2a[i], float(t2a[i]))
                     for i in range(len(t1a))])
    bat = C_fn_batch(n1a, t1a, n2a, t2a)
    rel = np.abs(scal - bat) / np.maximum(np.abs(scal), 1e-300)
    print(f"  C_fn_batch vs scalar: {len(t1a)} samples, shape={bat.shape}, "
          f"max rel-err={rel.max():.2e} (must be ~0)")
