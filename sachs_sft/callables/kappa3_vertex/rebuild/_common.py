"""Shared loaders for the factored kappa3 vertex rebuild.

`build_equal_time_limber_table.py` builds LOW and HIGH in one process, so
every (cutoff, bispectrum model, grid density) combination pays for the
expensive LOW branch again. Three properties of the calculation make a
factored build exact rather than approximate:

* rows are independent. zeta[i] depends only on cosine_triples[i]; verified
  bit-exactly for both branches by `tests/test_row_independence.py`. So one
  build on the densest triple set yields every coarser subset by selection.
* LOW is independent of the HIGH cutoff and of the bispectrum model. It is
  always tree-level, because its FFTlog machinery needs a separable
  bispectrum. So one LOW build serves every variant.
* the HIGH octave bands are disjoint in max-leg, so their sum is the same
  integral. One band set therefore yields every cutoff by partial summation.

These three together turn a build matrix into one LOW build plus one band
set per bispectrum model.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
#: sachs_sft/callables/kappa3_vertex, which holds the b_model package.
_CODE_ROOT = _HERE.parent
#: repository root, STF_lensing.
_REPO = _CODE_ROOT.parents[3]

_ETL = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
        / "kappa3_vertex" / "equal_time_limber")

#: The L2 NPZ that fixes the deployed lambda-shell grid.
ARCHIVED_L2 = _ETL / "inputs" / "l2_lambda_grid_z5covgrid16.npz"

CHANNELS = ("TTT", "TTP", "TPP", "PPP")

#: Linear P(k) table of every piece built before 2026-10-03 (sigma8 0.808988,
#: no header). A piece whose build metadata has no "pk_table" key used it.
PK_TABLE_DEFAULT = "PCAMBz0.txt"
#: Table of the two-point side (C table, FF, Order 0, PyCCL reference; sigma8
#: 0.810000, cosmology in its header). The T-001 rebuild uses it (FIX_LIST N11).
PK_TABLE_2PT = "PCAMB_pyccl_stf_fid_z0.txt"
PK_TABLES = (PK_TABLE_DEFAULT, PK_TABLE_2PT)

#: leg_order of the canoes HIGH builders (canoes 30f2e3b and later; default "auto").
LEG_ORDERS = ("auto", "given")


def _resolve_canoes_root() -> Path:
    env = os.environ.get("CANOES_ROOT")
    if env:
        return Path(env)
    for candidate in (Path("/Users/zzhang/projects/angular_statistics/canoes"),
                      Path("/Users/zzhang/projects/canoes")):
        if (candidate / "src" / "canoes").is_dir():
            return candidate
    raise RuntimeError("canoes checkout not found; set CANOES_ROOT")


def bootstrap_paths() -> None:
    """Put canoes, the production callable folder and b_model on sys.path."""
    root = _resolve_canoes_root()
    for p in (root / "src", root, _ETL, _CODE_ROOT):
        if p.is_dir() and str(p) not in sys.path:
            sys.path.insert(0, str(p))


def pk_table_path(pk_table: str) -> Path:
    """Path of one of the known CAMB tables in canoes examples/data."""
    if pk_table not in PK_TABLES:
        raise ValueError(f"unknown P(k) table {pk_table!r}; known: {PK_TABLES}")
    return _resolve_canoes_root() / "examples" / "data" / pk_table


def load_cosmology(pk_table: str = PK_TABLE_DEFAULT):
    """The fiducial cosmology and the linear P(k) of one CAMB table.

    The default is the table the deployed (pre-2026-10-03) pieces were built on.
    """
    from _local_cosmo_pk import FiducialCosmology, load_pk_delta

    cosmo = FiducialCosmology()
    return cosmo, load_pk_delta(n_s=cosmo.n_s, table_path=pk_table_path(pk_table))


def pk_table_record(pk_table: str) -> dict:
    """Name and sha256 of the CAMB table, for a piece's build metadata."""
    digest = hashlib.sha256(pk_table_path(pk_table).read_bytes()).hexdigest()
    return {"pk_table": pk_table, "pk_table_sha256": digest}


def build_pk_table(build: dict) -> str:
    """The P(k) table a piece was built on, from its build metadata."""
    return str(build.get("pk_table", PK_TABLE_DEFAULT))


def build_pk_identity(build: dict) -> tuple[str, str | None]:
    """Validated basename and hash, with ``None`` identifying legacy provenance.

    A known hash and a missing hash are distinct identities, even when the
    basenames agree. This does not require the LOW and HIGH identities to match.
    """
    table = build_pk_table(build)
    if (not isinstance(build.get("pk_table", PK_TABLE_DEFAULT), str)
            or not table or table in {".", ".."}
            or any(character in table for character in ("/", "\\", ":", "\0"))):
        raise ValueError("P(k) table metadata must contain a nonempty basename")
    digest = build.get("pk_table_sha256")
    if digest is not None and (
        not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    ):
        raise ValueError("P(k) table SHA-256 must be 64 lowercase hexadecimal characters")
    return table, digest


def build_leg_order(build: dict) -> str:
    """The leg order a piece was built with.

    Pieces without the key predate canoes 30f2e3b, which evaluated every row in
    the given order.
    """
    return str(build.get("leg_order", "given"))


def canoes_record() -> dict:
    """canoes HEAD and whether its src/ differs from HEAD, for a piece's build metadata."""
    import subprocess

    root = _resolve_canoes_root()
    head = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--", "src"],
                           capture_output=True, text=True, check=True).stdout.strip()
    return {"canoes_commit": head, "canoes_src_clean": not dirty}


def load_lambda_grid(l2_npz: Path, h: float) -> tuple[np.ndarray, np.ndarray]:
    """(lambda physical Mpc, lambda in Mpc/h) from an L2 NPZ."""
    d = np.load(l2_npz, allow_pickle=False)
    lam = np.asarray(d["lambda_grid"], dtype=np.float64)
    if lam.ndim != 1 or not np.all(np.diff(lam) > 0):
        raise ValueError(f"{l2_npz}: lambda_grid must be strictly increasing 1-D")
    return lam, lam * float(h)


def load_triples(path: Path) -> np.ndarray:
    return np.asarray(np.load(path)["cosine_triples"], dtype=np.float64)


def octave_bands(low: int, high: int) -> list[tuple[int, int]]:
    """Disjoint octave windows covering (low, high] in max-leg."""
    bands, edge = [], int(low)
    while edge < high:
        upper = min(2 * edge, int(high))
        bands.append((edge, upper))
        edge = upper
    return bands


def dump_meta(meta: dict) -> np.ndarray:
    """Serialise a metadata dict to the uint8 byte array canoes uses."""
    return np.frombuffer(json.dumps(meta, default=str).encode(), dtype=np.uint8)


def load_meta(arr: np.ndarray) -> dict:
    return json.loads(np.asarray(arr, dtype=np.uint8).tobytes().decode())


def grid_fingerprint(triples: np.ndarray, lam: np.ndarray) -> dict:
    """Identity of the grid a piece was built on, checked at assembly time.

    Pieces built on different grids must never be summed. The fingerprint is
    compared exactly, so a piece from a different triple set or shell grid
    fails loudly instead of producing a silently wrong table.
    """
    return dict(
        n_triples=int(triples.shape[0]),
        n_shells=int(lam.size),
        triples_hash=int(
            np.frombuffer(
                np.ascontiguousarray(np.round(triples, 12)).tobytes(),
                dtype=np.uint64,
            ).sum(dtype=np.uint64)
        ),
        lam_first=float(lam[0]),
        lam_last=float(lam[-1]),
    )
