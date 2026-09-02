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

import json
import os
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
#: driver_field_emulators/code, which holds the b_model package.
_CODE_ROOT = _HERE.parent
#: repository root, STF_lensing.
_REPO = _CODE_ROOT.parents[1]

_ETL = (_REPO / "SFT-lensing-paper-analyses" / "sachs_sft" / "callables"
        / "kappa3_vertex" / "equal_time_limber")

#: The L2 NPZ that fixes the deployed lambda-shell grid.
ARCHIVED_L2 = (
    _REPO / "SFT-lensing-paper-analyses" / "_archive" / "canoes_pipeline"
    / "inputs" / "l2_callables_2026_05_26"
    / "windowed_squeezed_kappa3_lmax1000_z5_covgrid16_omega0316_h06711_physMpc.npz"
)

CHANNELS = ("TTT", "TTP", "TPP", "PPP")


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


def load_cosmology():
    """The fiducial cosmology and linear P(k) the deployed table was built on."""
    from _local_cosmo_pk import FiducialCosmology, load_pk_delta

    cosmo = FiducialCosmology()
    return cosmo, load_pk_delta(n_s=cosmo.n_s)


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
