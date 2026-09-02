"""A drop-in replacement for the production corr_op C callable that RECORDS
every (t1, t2) pair sft-wick asks for, then forwards to the real callable.

Writes one .npy per process to $FOLD_BIAS_REC_DIR at interpreter exit, so it
works under joblib/loky too.
"""
from __future__ import annotations

import atexit
import os
import sys
from pathlib import Path

import numpy as np

_PROD = Path(
    "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/"
    "sachs_sft/callables/C_propagator/corr_op"
)
if str(_PROD) not in sys.path:
    sys.path.insert(0, str(_PROD))

import corr_op_C_callable as _prod  # noqa: E402

TABLE_PATH = _prod.TABLE_PATH
QUERY_COORDINATE = _prod.QUERY_COORDINATE
QUERY_UNITS = _prod.QUERY_UNITS
COMPONENT_ORDER = _prod.COMPONENT_ORDER
OUTPUT_H_POWER = _prod.OUTPUT_H_POWER

_REC: list[np.ndarray] = []
_OUT = Path(os.environ.get("FOLD_BIAS_REC_DIR", "/tmp/fold_bias_rec"))


def _flush() -> None:
    if not _REC:
        return
    _OUT.mkdir(parents=True, exist_ok=True)
    arr = np.concatenate(_REC, axis=0)
    np.save(_OUT / f"rec_{os.getpid()}.npy", arr)


atexit.register(_flush)


def _log(n1, t1, n2, t2) -> None:
    t1a = np.atleast_1d(np.asarray(t1, dtype=float)).ravel()
    t2a = np.atleast_1d(np.asarray(t2, dtype=float)).ravel()
    n = max(t1a.size, t2a.size)
    t1a = np.broadcast_to(t1a, (n,))
    t2a = np.broadcast_to(t2a, (n,))
    a1 = np.atleast_2d(np.asarray(n1, dtype=float))
    a2 = np.atleast_2d(np.asarray(n2, dtype=float))
    cg = np.broadcast_to(np.sum(a1 * a2, axis=-1).ravel(), (n,))
    _REC.append(np.stack([t1a, t2a, cg], axis=-1))


def C_fn(n1, t1, n2, t2):
    _log(n1, t1, n2, t2)
    return _prod.C_fn(n1, t1, n2, t2)


def C_fn_batch(n1, t1, n2, t2):
    _log(n1, t1, n2, t2)
    return _prod.C_fn_batch(n1, t1, n2, t2)


for _f in (C_fn, C_fn_batch):
    _f.query_coordinate = QUERY_COORDINATE
    _f.query_units = QUERY_UNITS
    _f.component_order = COMPONENT_ORDER
