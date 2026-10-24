"""Guarded z_s = 3.0 replay of the T-001 cutoff-ladder rungs.

Candidate z_s=3 finite products, no scientific acceptance.

Subcommands
  source-map  write the z_s=3.0 endpoint record (new file, published complete or not at all)
  prepare     materialise config, loader and contract per target; no fold
  run         one serial fold under the shared lock, with started/completed records;
              a z_s=3 fold first needs a PASS z_s=5 gate report
  gate        compare z_s=5 replays with the accepted T-001 products (1e-12 relative)
  assemble    write the z_s=3 ladder members and assembly_manifest.json

The reviewed modules (main_guard.py, continuation_guard.py, run_fk_variant.py,
replay_suite.py, aligned_response.py and the review batch driver) are imported
by path and never edited. Every generated file is a new file in a fresh
directory. Every record binds the sha256 of this file, so an edit to it
invalidates the source map, the prepared batches and the gate report written
before the edit. Run with the sft-wick interpreter, -B, PYTHONDONTWRITEBYTECODE=1,
PYTHONPATH=<frozen framework export>:$CANOES_ROOT/src and CANOES_ROOT set.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import ctypes
import dataclasses
import difflib
import errno
import hashlib
import importlib.util
import json
import math
import os
import re
import resource
import secrets
import subprocess
import sys
import time
import zipfile
from collections.abc import Callable, Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = Path(__file__).resolve().parent
DRIVER = Path(__file__).resolve()
PKG = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses")
SACHS = PKG / "sachs_sft"
R1 = SACHS / "analyses/r1_sft061"
PLANNED_HOME = R1 / "zs3p0_ladder_t001_20261023"
CAND = R1 / "t001_candidate_inputs"
V6 = CAND / "main_fk_gl36_optimized_v6"
FW = V6 / "framework_export/src"
MAIN_EXEC = V6 / "execution/main_fk_gl36"
REVIEW = CAND / "point6_preflight/review_batch"
CONT = CAND / "point6_preflight/continuation_batch"
LADDER_V8 = CAND / "final_ladder_v8_20261005/assembly_manifest.json"
ACCEPTED = R1 / "t001_accepted"
KV = SACHS / "callables/kappa3_vertex"
PIECES = KV / "rebuild/products/pieces_t001"
PCHIP = KV / "equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable_pchip.py"
REFINED = PIECES / "table_permclosed_np512_t001_refined.npz"
CANOES = Path("/Users/zzhang/projects/angular_statistics/canoes")
C23_TABLE = CANOES / "docs/handoff/2026_10_03_T-001_delivery/corr_op_table_zs5_23pt_stf_fid_t001.npz"
PYTHON = Path("/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python")
LOCK = R1 / ".replay_suite.lock"
REVIEWED_SOURCE_MAP = R1 / "source_map.json"
EXPECTED_PYTHONPATH = f"{FW}{os.pathsep}{CANOES / 'src'}"
FRAMEWORK_COMMIT = "5198fb9845334982d103059ac7b84f7afb2228ca"
CANOES_COMMIT = "53850f3423ee25b9c1813cff8c2a2b282c017c70"
EMPTY_DIFF_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
MAIN_GEOMETRY = "6102de8ee3fef96fcafda2773046e71f26394b55c1b15d1f56513e9e567996d3"
LADDER_GEOMETRY = "bad69362cf45bab24e366e1dadd3d0d1a49f3a57b5b24bc819a69acd7cb66188"
#: Source planes this driver serves: the z_s=5 gates on the reviewed source map and
#: the z_s=3 production endpoint on the map this driver writes.
GATE_Z = 5.0
PRODUCTION_Z = 3.0
REVIEWED_Z5_LAMBDA = 2317.9695958203943
EXPECTED_LAMBDA = {PRODUCTION_Z: 2255.3607197066854}
MAX_AFFINE_MPC = 2353.929874188644
ROUNDTRIP_ATOL = 5e-13
MAIN_SUBSET = [0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
N_LADDER_ROWS = len(MAIN_SUBSET)
GL_ACCEPTED = 36
TEMPLATE_GL = 24
REL_TOL = 1e-12
MULTIZ_GAMMA_ATOL_ARCMIN = 1e-3
ORDERS = "0,2"
SCOPE = "candidate z_s=3 finite products, no scientific acceptance"
MARKER = 'TABLE_PATH = _HERE / "table_permclosed.npz"'
TAG_PATTERN = re.compile(r"[A-Za-z0-9_]+")
GATE_WORK_NAME = "gates_zs5_gl36"
#: renamex_np flag (macOS <stdio.h>): fail with EEXIST instead of replacing the target.
RENAME_EXCL = 0x00000004
#: Measured peak sampled group RSS per sweep worker count (T-001 GL36 records:
#: main v6 resource_completion.json for 4; final_ladder_v8 admissions for 8 and 12).
MEASURED_PEAK_BYTES = {4: 27544076288, 8: 50694569984, 12: 67749150720}
MEMORY_FLOOR_BYTES = 12 * 2**30
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")

#: Reviewed or accepted files whose bytes this driver depends on.
PINNED = {
    R1 / "aligned_response.py": "888f1593e55342f088ac486e7709688853b6827506d8d43f90ec35ee10e1c5f4",
    REVIEWED_SOURCE_MAP: "83a926ae4d52dc5435dbe9ab7d016dddbbbf25b149c5a640e419afcd842bdc2b",
    R1 / "replay_suite.py": "9a829970732753dc9496ecd1fcfe864f18a9718cc8ddcfd447af086040f41117",
    PKG / "reproduce/archived_inputs.py": "1f9f0297ed8014f09d8c2350259af58fe627300a70fed970792b5b19f4509dde",
    REVIEW / "main_guard.py": "0e55eef3dbdaae77f00dc135f305a090a01aea7950cd83d99a10f44231eb9e63",
    REVIEW / "batch_driver.py": "4fcb3b4034c19a3a37f007c67c1a12e7e41eae7f8f39c4bad6a066ed11b0c20d",
    REVIEW / "review_plan.json": "62f8b6f019d374e28e226101e45f9b0ebba8ed052006077d160e11f366a99883",
    CONT / "continuation_guard.py": "49a56a5b99b53175258925fcf65ae461516a0d7a5691b35513b822c416d0c15a",
    PCHIP: "b98d25c5f629ca79ebbae25bde05fd94c630bdfeff1b3efcaf158c9c9e42512d",
    REFINED: "5ba7497b76b3fdf1c9ccc13efb06263be2204dc7054b7867699975a406e62f8d",
    PIECES / "ladder_manifest.json": "9b6ee5403adbd3f56a5027be1c752ff139239b43836068c15f8bee793996dd5f",
    PIECES / "ladder/table_tree_auto_cut7680.npz": "4711eb7122f641ba1c27200f98fb093d3edd428d7ed7868ffe9d6dff999bd0ac",
    PIECES / "ladder/table_tree_auto_cut15360.npz": "dd8bb15e62cca0b72973b753d5081deb675c85178b8885f483c4a92f79b65deb",
    PIECES / "ladder/table_bihalofit_auto_cut7680.npz":
        "4c354769cd4ece2b80fe0eb20a95105cb4687ffcac11ec4d75444e590403b1bc",
    PIECES / "ladder/table_bihalofit_auto_cut15360.npz":
        "97fb0d7ab4d75f5c57057285ddfb4d2a5cf39b235faa1c360306faa681fe75a7",
    PIECES / "ladder/table_bihalofit_given_cut7680.npz":
        "e3a378393e76bb3815fcd1f9ed2ddf83cf64bdd96e98b92680b6bda6eb748d13",
    PIECES / "ladder/table_bihalofit_given_cut15360.npz":
        "dcf8fd5937dc96dc587de10a0d58ee9cdb117490c1fff8d90e6c3af040fb6a9d",
    CAND / "corr_op_C_callable.py": "1d0c80839d7a3f3866b5fff36d07b9e08d945d7005cc5b5ec5183b2c21906dd2",
    C23_TABLE: "d487c51b9dcf1d87520d82eac7d37aee17b872b3c83356f8007ff82d472562ee",
    SACHS / "scripts/inputs/F_tensor.npy": "63b682be4315a5783d1554e7ea48364229e8bcf17181c6482d14269101e1d0db",
    SACHS / "scripts/kappa2_callable.py": "49372414bab52dce803f786c301c07219b638b39f5038f60e0bb9858ce18cb29",
    SACHS / "scripts/inputs/kappa2_grid_lcut0_low_empty.npz":
        "1e03b492d46ed93db3a4ad6d953323325a20322470371412f1ac3587b1139e63",
    SACHS / "scripts/run_FF_single.py": "5ca112bbcad3e38ff8ee6f6c6544661dc962eeeff7aea3423a0f5dd9abc7e443",
    SACHS / "scripts/spin_rotation.py": "0322515d53b0dcc8832e585fb8bed3b242d465873e79f71440adf7e0f71abc8c",
    KV / "rebuild/run_fk_variant.py": "0ea02c4ce96cd77f1c146a5cc38214a059530b8833f8a61b005b490ec84407b2",
    CAND / "configs/tree_7680.yaml": "a1e4a5e96884b8459a97dc22a4db603588c1f16bfba6737aacd247bf34bfeb02",
    CAND / "configs/tree_15360.yaml": "3dc8972bd90ad084e4b6d0c56a0aa5698d4db87e8f157b881119c3bc6b3e79b9",
    CAND / "configs/bihalofit_7680.yaml": "a10f20375fc87deee3132d8a7dfc1e66037be9469f6bef65111e7d560dc17c94",
    CAND / "configs/bihalofit_15360.yaml": "754a7e287fe999867a69c49100fd11baf8d71ab100f5176276fe92c932a396a9",
    CANOES / "src/canoes/cosmo/lightcone.py": "69f57f5e95545c094a4ad95e70a603635f5a886178df5e82703191a849fb3d02",
    CANOES / "src/canoes/cosmo/growth.py": "c931cca5f6d66b4b3082a6a0a11924604d12583219bbe31ca8df804cd371e08c",
    PYTHON.resolve(): "17f40a0285f0eddfc7eb226d625390ecd8fbd3ed907ae37bb98e81881ee95784",
    V6 / "source_snapshot.json": "19d794073859d9addcd153039f9cb3b6198a5871fe4a179935b1260143105d16",
    MAIN_EXEC / "prepared.json": "e319cc4cc695bbdfde15259b1fc85225921195775d8205cf460aff7ca9004dd6",
    MAIN_EXEC / "config_main.yaml": "60e3268fe80170cc594f3d186b92f2454b5ad1bcf277b6918324ab44ad75ef20",
    MAIN_EXEC / "callable_main.py": "65e136196e80df18a9f1aee47d7d9742f6069b4a9723fc0872bd3b8c548f6cde",
    MAIN_EXEC / "main_contract.json": "f39dd012e356968328074981180d40824d23a97be8ba886cd823ba3206f07e86",
    MAIN_EXEC / "completed.json": "e88e2d2c37227dbfca82a86a4fd963720576277f9285ff5a0268d701ea73ed50",
    MAIN_EXEC / "xi_main.npz": "d4bbc17d26db1ca50c27c8316238edcca02c998e8882ab88116be371b1d10e3a",
    LADDER_V8: "088bb80d110ee183ca53e3b99b8510d7ea2069ee448e26fe90fd8a3501d66a8a",
    ACCEPTED / "cutoff_ladder/bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz":
        "eb504b9fa7bd66fdda1b64cff3d9c7fd1aeef4b99372c7fcb37762a45794f792",
    ACCEPTED / "cutoff_ladder/tree_auto/table_tree_cut15360_r4_xi.npz":
        "e89b13a65e0608d6f65df5588a1c503f97891dd0920fd6458a572b1bd911d5eb",
    ACCEPTED / "products/main_fk.npz": "d2c6f5c434c2ab6970cea7e42091e3e7ab9f0424a9c63e7efb801124809137e3",
    ACCEPTED / "products/multiz.npz": "4bc94805b5c57939a55454a1c706ec8f718d848c8872a5ecffaabed8b22f71cb",
}
#: Fold inputs pinned in every contract, in addition to tables, templates,
#: the source map, this driver, the framework export and canoes sources.
FOLD_INPUTS = [R1 / "aligned_response.py", R1 / "replay_suite.py", PKG / "reproduce/archived_inputs.py",
               REVIEW / "main_guard.py", REVIEW / "batch_driver.py", CONT / "continuation_guard.py",
               PCHIP, PIECES / "ladder_manifest.json", CAND / "corr_op_C_callable.py", C23_TABLE,
               SACHS / "scripts/inputs/F_tensor.npy", SACHS / "scripts/kappa2_callable.py",
               SACHS / "scripts/inputs/kappa2_grid_lcut0_low_empty.npz", SACHS / "scripts/run_FF_single.py",
               SACHS / "scripts/spin_rotation.py", KV / "rebuild/run_fk_variant.py", PYTHON.resolve()]
REFERENCE_FILES = sorted(set(PINNED) - set(FOLD_INPUTS), key=str)


@dataclasses.dataclass(frozen=True)
class Spec:
    """One ladder rung: model, high-ell leg order, cutoff, table guard scope and accepted analogue."""

    model: str
    leg: str
    cutoff: int
    scope: str
    table: Path
    template: str
    accepted: str | None
    member: str | None


TARGETS = {
    "tree_main_15360": Spec("tree", "auto", 15360, "main", REFINED, "tree_15360", "main_fk_gl36",
                            "tree_auto/table_tree_cut15360_r4_xi.npz"),
    "tree_auto_7680": Spec("tree", "auto", 7680, "ladder", PIECES / "ladder/table_tree_auto_cut7680.npz",
                           "tree_7680", "tree_auto_7680_gl36", "tree_auto/table_tree_cut7680_r4_xi.npz"),
    "bihalofit_auto_7680": Spec("bihalofit", "auto", 7680, "ladder",
                                PIECES / "ladder/table_bihalofit_auto_cut7680.npz", "bihalofit_7680",
                                "bihalofit_auto_7680_gl36", "bihalofit_auto/table_bihalofit_cut7680_r4_xi.npz"),
    "bihalofit_auto_15360": Spec("bihalofit", "auto", 15360, "ladder",
                                 PIECES / "ladder/table_bihalofit_auto_cut15360.npz", "bihalofit_15360",
                                 "bihalofit_auto_15360_gl36", "bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz"),
    # Optional controls (plan step 8). They are never assembled.
    "tree_auto_15360": Spec("tree", "auto", 15360, "ladder", PIECES / "ladder/table_tree_auto_cut15360.npz",
                            "tree_15360", None, None),
    "bihalofit_given_7680": Spec("bihalofit", "given", 7680, "ladder",
                                 PIECES / "ladder/table_bihalofit_given_cut7680.npz", "bihalofit_7680",
                                 "bihalofit_given_7680_gl36", None),
    "bihalofit_given_15360": Spec("bihalofit", "given", 15360, "ladder",
                                  PIECES / "ladder/table_bihalofit_given_cut15360.npz", "bihalofit_15360",
                                  "bihalofit_given_15360_gl36", None),
}
PRODUCTION = ("tree_main_15360", "tree_auto_7680", "bihalofit_auto_7680", "bihalofit_auto_15360")
GATES = ("bihalofit_auto_15360", "tree_main_15360")
GATE_MATCH_KEYS = ("table_sha256", "template_config", "n_jobs", "gl")


class Refusal(RuntimeError):
    """A guard refused the requested action."""


# --------------------------------------------------------------------- basics
def digest(path: Path) -> str:
    """sha256 of a file's bytes."""
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def record(path: Path) -> dict:
    """Resolved path and sha256 of a file."""
    return {"path": str(Path(path).resolve()), "sha256": digest(path)}


def canonical(data: Any) -> bytes:
    """Canonical JSON bytes (sorted keys, no whitespace); the hash basis of every record."""
    return json.dumps(data, sort_keys=True, separators=(",", ":")).encode()


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def check_records(records: list[dict]) -> None:
    """Refuse if a recorded file is missing, unreadable or changed."""
    for item in records:
        try:
            current = digest(Path(item["path"]))
        except OSError as exc:
            raise Refusal(f"Pinned file unreadable: {item['path']} ({exc})") from exc
        if current != item["sha256"]:
            raise Refusal(f"Pinned file changed: {item['path']}")


def verify_pinned() -> None:
    """Every reviewed or accepted file the driver depends on still has its pinned sha256."""
    check_records([{"path": str(path), "sha256": expected} for path, expected in PINNED.items()])


def read_json(path: Path) -> Any:
    """Parse a JSON record; a missing or truncated file refuses."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Refusal(f"Cannot read JSON record {path}: {exc}") from exc


def exclusive_rename(src: Path, dst: Path) -> None:
    """Atomic rename that fails if ``dst`` exists, file or directory (macOS renamex_np RENAME_EXCL)."""
    function = getattr(ctypes.CDLL(None, use_errno=True), "renamex_np", None)
    if function is None:
        raise Refusal("renamex_np(RENAME_EXCL) is unavailable; this driver runs on macOS only")
    function.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
    function.restype = ctypes.c_int
    if function(os.fsencode(src), os.fsencode(dst), RENAME_EXCL) != 0:
        code = ctypes.get_errno()
        if code == errno.EEXIST:
            raise FileExistsError(code, os.strerror(code), str(dst))
        raise OSError(code, os.strerror(code), str(dst))


def publish_new(path: Path, data: bytes, validate: Callable[[Path], object] | None = None) -> None:
    """Create ``path`` holding ``data``: complete or absent, never replacing anything.

    The bytes go to a hidden temporary file in the same directory, are flushed to
    disk and optionally validated there, and are then renamed into place with
    RENAME_EXCL, which fails if ``path`` exists. The temporary name is removed in
    every case.
    """
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise Refusal(f"Refusing to replace an existing file: {path}")
    tmp = path.parent / f".{path.name}.{os.getpid()}.{secrets.token_hex(4)}.partial"
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        if validate is not None:
            validate(tmp)
        exclusive_rename(tmp, path)
    except FileExistsError as exc:
        raise Refusal(f"Refusing to replace an existing file: {path}") from exc
    finally:
        tmp.unlink(missing_ok=True)


def write_new_text(path: Path, text: str) -> None:
    publish_new(path, text.encode("utf-8"))


def write_new_bytes(path: Path, data: bytes) -> None:
    publish_new(path, data)


def write_new_json(path: Path, value: dict, validate: Callable[[Path], object] | None = None) -> None:
    """Same bytes as replay_suite.dump_new (indent 2, trailing newline), published atomically."""
    publish_new(path, (json.dumps(value, indent=2) + "\n").encode("utf-8"), validate)


def inside(path: Path, parent: Path) -> bool:
    """``path`` lies strictly below ``parent`` after resolving symlinks."""
    return path.resolve().is_relative_to(parent.resolve()) and path.resolve() != parent.resolve()


# ------------------------------------------------------ reviewed modules by path
_MODULES: dict[str, ModuleType] = {}


def load_module(name: str, path: Path, extra_path: Path | None = None) -> ModuleType:
    """Import a reviewed module by explicit path, executing the same bytes whose sha256 was checked."""
    if name in _MODULES:
        return _MODULES[name]
    try:
        source = path.read_bytes()
    except OSError as exc:
        raise Refusal(f"Cannot read reviewed module {path}: {exc}") from exc
    if sha_bytes(source) != PINNED[path]:
        raise Refusal(f"Reviewed module changed: {path}")
    if extra_path is not None and str(extra_path) not in sys.path:
        sys.path.insert(0, str(extra_path))
    spec = importlib.util.spec_from_loader(name, loader=None, origin=str(path))
    if spec is None:
        raise Refusal(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(path)
    sys.modules[name] = module
    try:
        # dont_inherit: the module must not pick up this driver's __future__ flags.
        exec(compile(source, str(path), "exec", dont_inherit=True), module.__dict__)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    _MODULES[name] = module
    return module


def shared() -> ModuleType:
    """Review batch driver: check_active_folds, fold_lock, validate_result, read_geometry."""
    return load_module("zs3_t001_review_batch_driver", REVIEW / "batch_driver.py", REVIEW)


def continuation_guard() -> ModuleType:
    return load_module("zs3_continuation_guard", CONT / "continuation_guard.py")


def replay() -> ModuleType:
    return load_module("zs3_replay_suite", R1 / "replay_suite.py")


def variant() -> ModuleType:
    return load_module("zs3_run_fk_variant", KV / "rebuild/run_fk_variant.py")


def aligned() -> ModuleType:
    if Path(os.environ.get("CANOES_ROOT", "")) != CANOES:
        raise Refusal("CANOES_ROOT must equal the recorded canoes root before importing aligned_response")
    return load_module("zs3_aligned_response", R1 / "aligned_response.py")


# ---------------------------------------------------------- runtime binding
def require_invocation() -> None:
    """The caller is the sft-wick interpreter with the pinned PYTHONPATH, CANOES_ROOT and no bytecode."""
    problems = []
    if Path(sys.executable).resolve() != PYTHON.resolve():
        problems.append(f"interpreter {sys.executable} is not {PYTHON}")
    if os.environ.get("PYTHONDONTWRITEBYTECODE") != "1" or not sys.flags.dont_write_bytecode:
        problems.append("set PYTHONDONTWRITEBYTECODE=1 (the documented commands also pass -B)")
    if os.environ.get("PYTHONPATH") != EXPECTED_PYTHONPATH:
        problems.append(f"PYTHONPATH must be exactly {EXPECTED_PYTHONPATH}")
    if os.environ.get("CANOES_ROOT") != str(CANOES):
        problems.append(f"CANOES_ROOT must be {CANOES}")
    if problems:
        raise Refusal("; ".join(problems))


def framework_state() -> dict:
    """The frozen 5198fb9 export must be what sft_wick resolves to, byte for byte."""
    found = importlib.util.find_spec("sft_wick")
    origin = Path(found.origin).resolve() if found is not None and found.origin else None
    if origin != FW / "sft_wick/__init__.py":
        raise Refusal(f"sft_wick resolves to {origin}, not the frozen export {FW}")
    snapshot = read_json(V6 / "source_snapshot.json")
    if snapshot["commit"] != FRAMEWORK_COMMIT or Path(snapshot["framework_src"]) != FW:
        raise Refusal("Framework source snapshot does not record the frozen 5198fb9 export")
    roster = {item["export"]["path"]: item["export"]["sha256"] for item in snapshot["source_roster"]}
    entries = sorted(FW.rglob("*"))
    files = {str(p): digest(p) for p in entries if p.is_file()}
    stray = [str(p) for p in entries if not p.is_file() and not p.is_dir() or p.name == "__pycache__"]
    if stray or files != roster:
        raise Refusal("Framework export differs from its recorded roster")
    return {"commit": FRAMEWORK_COMMIT, "framework_src": str(FW), "files": len(files),
            "source_snapshot": record(V6 / "source_snapshot.json")}


def runtime_state() -> dict:
    """Dependency versions of the fold interpreter equal the T-001 actual_runtime record."""
    actual = read_json(MAIN_EXEC / "prepared.json")["actual_runtime"]
    identity = shared().runtime_identity()
    if (identity["dependencies"] != actual["dependencies"]
            or identity["python_version"] != actual["python_version"]
            or identity["executable"] != actual["executable"]
            or actual["framework_commit"] != FRAMEWORK_COMMIT or Path(actual["framework_src"]) != FW):
        raise Refusal("Interpreter or dependencies differ from the T-001 actual_runtime record")
    return {"python": str(PYTHON), "python_version": identity["python_version"],
            "dependencies": identity["dependencies"], "actual_runtime_record": record(MAIN_EXEC / "prepared.json")}


def canoes_state() -> dict:
    """canoes HEAD and tracked diff equal the T-001 review record (53850f3, clean)."""
    state = shared().git_state(CANOES)
    recorded = read_json(REVIEW / "review_plan.json")["canoes"]
    if (state["commit"] != CANOES_COMMIT or state["tracked_diff_sha256"] != EMPTY_DIFF_SHA256
            or state != recorded):
        raise Refusal(f"canoes state {state} differs from the T-001 record {recorded}")
    return state


def fold_environment(root: Path) -> dict:
    """Caller environment without SFT_WICK_* or PYTHON* variables, plus the pinned fold variables."""
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SFT_WICK_", "PYTHON"))}
    env.update(PYTHONPATH=EXPECTED_PYTHONPATH, CANOES_ROOT=str(CANOES), PYTHONDONTWRITEBYTECODE="1",
               PYTHONNOUSERSITE="1", MPLCONFIGDIR=str(root / "mpl_cache"))
    env.update({name: "1" for name in THREAD_VARS})
    return env


def subprocess_binding(env: dict, cwd: Path) -> None:
    """Import resolution inside the fold environment, with run_FF_single's path order.

    find_spec resolves without executing the packages, so no matplotlib cache is created.
    """
    code = ("import sys, importlib.util as u\nfrom pathlib import Path\n"
            "for d in sys.argv[3:5]:\n"
            "    sys.path.remove(d) if d in sys.path else None\n    sys.path.insert(0, d)\n"
            "for name, expected in (('sft_wick', sys.argv[1]), ('canoes', sys.argv[2])):\n"
            "    spec = u.find_spec(name)\n"
            "    assert spec is not None and Path(spec.origin).resolve() == Path(expected), (name, spec)\n"
            "print('BINDING_OK')\n")
    argv = [str(PYTHON), "-B", "-c", code, str(FW / "sft_wick/__init__.py"),
            str(CANOES / "src/canoes/__init__.py"), str(PKG.parent), str(SACHS / "scripts")]
    done = subprocess.run(argv, env=env, cwd=cwd, capture_output=True, text=True, check=False)
    if done.returncode != 0 or done.stdout.strip() != "BINDING_OK":
        raise Refusal(f"Fold-environment import binding failed: {done.stderr.strip()[-2000:]}")


# ------------------------------------------------------------ source endpoint
def tag_for(z: float) -> str:
    return "zs" + repr(float(z)).replace(".", "p").replace("-", "m")


def require_tag(tag: Any) -> str:
    if not isinstance(tag, str) or not TAG_PATTERN.fullmatch(tag):
        raise Refusal(f"Invalid target tag {tag!r}")
    return tag


def read_source_map_checked(path: Path) -> dict:
    """replay_suite.read_source_map (the reviewed rules); a failure refuses."""
    try:
        return replay().read_source_map(Path(path))
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise Refusal(f"Source map {path} fails replay_suite.read_source_map: {exc}") from exc


def validate_source_map_file(path: Path) -> None:
    read_source_map_checked(path)


def source_map_command(z: float, out: Path) -> int:
    """Write the z_s=3.0 endpoint record; it is validated before it appears under its name."""
    require_invocation()
    z = float(z)
    if z not in EXPECTED_LAMBDA:
        raise Refusal(f"This driver writes source maps only for z_s in {sorted(EXPECTED_LAMBDA)}")
    out = Path(os.path.abspath(out)).resolve()
    if not inside(out, HERE) or out.exists() or not out.parent.is_dir():
        raise Refusal(f"--out must be a new file inside {HERE} with an existing parent: {out}")
    verify_pinned()
    canoes = canoes_state()
    endpoint, checks = derive_endpoint(z, validate=True)
    reviewed = read_json(REVIEWED_SOURCE_MAP)
    sources = [R1 / "aligned_response.py", CANOES / "src/canoes/cosmo/lightcone.py",
               CANOES / "src/canoes/cosmo/growth.py", CANOES / "src/canoes/sachs/sft_input/corr_op/table.py",
               CANOES / "src/canoes/sachs/sft_input/corr_op/_limber.py"]
    provenance = dict(reviewed["mapping_provenance"])
    provenance.update(canoes_commit=canoes["commit"], canoes_tracked_diff_sha256=canoes["tracked_diff_sha256"],
                      source_sha256={str(p): digest(p) for p in sources},
                      C_table_path=str(C23_TABLE), C_table_sha256=digest(C23_TABLE),
                      C_wrapper=record(CAND / "corr_op_C_callable.py"),
                      source_mapping=("lambda_s = np.interp(z, _Z, _LAMBDA) on the cached canoes knots, the rule "
                                      "of aligned_response.write_source_map"),
                      reviewed_source_map=record(REVIEWED_SOURCE_MAP))
    payload = {"z": [endpoint["z"]], "lambda_mpc": [endpoint["lambda_mpc"]], "chi_mpc": [endpoint["chi_mpc"]],
               "D_mpc": [endpoint["D_mpc"]], "target_tag": tag_for(z), "mapping_provenance": provenance,
               "response_module": str(R1 / "aligned_response.py"), "response_attr": "response_R",
               "validation": checks, "producer": {"driver": record(DRIVER), "command": "source-map"},
               "scope": ("Candidate source endpoint for the T-001 ladder replay. The reviewed source_map.json "
                         "is unchanged and remains pinned by its own contracts.")}
    write_new_json(out, payload, validate=validate_source_map_file)
    print(json.dumps({"path": str(out), "sha256": digest(out), "z": payload["z"],
                      "lambda_mpc": payload["lambda_mpc"], "chi_mpc": payload["chi_mpc"],
                      "D_mpc": payload["D_mpc"], "checks": {k: v for k, v in checks.items()
                                                            if k != "aligned_response_validate"}}, indent=2))
    return 0


def derive_endpoint(z: float, validate: bool = False) -> tuple[dict, dict]:
    """Recompute lambda_s(z) and reproduce the reviewed six endpoints bitwise."""
    ar = aligned()
    reviewed = read_json(REVIEWED_SOURCE_MAP)
    six = np.interp(ar.SOURCE_REDSHIFTS, ar._Z, ar._LAMBDA)
    exact = {"z": ar.SOURCE_REDSHIFTS.tolist() == reviewed["z"],
             "lambda_mpc": six.tolist() == reviewed["lambda_mpc"],
             "chi_mpc": np.asarray(ar.chi_of_lambda(six)).tolist() == reviewed["chi_mpc"],
             "D_mpc": np.asarray(ar.D_at(six)).tolist() == reviewed["D_mpc"]}
    max_affine = float(ar.MAX_AFFINE_MPC)
    if (not all(exact.values()) or max_affine != MAX_AFFINE_MPC
            or reviewed["mapping_provenance"]["max_affine_mpc"] != MAX_AFFINE_MPC):
        raise Refusal(f"aligned_response no longer reproduces the reviewed source map: {exact}, {max_affine}")
    z = float(z)
    lam = float(np.interp(z, ar._Z, ar._LAMBDA))
    chi, dist, back = float(ar.chi_of_lambda(lam)), float(ar.D_at(lam)), float(ar.z_of_lambda(lam))
    if not (0.0 < lam < max_affine) or abs(back - z) > ROUNDTRIP_ATOL:
        raise Refusal(f"Endpoint for z={z} fails the round trip: lambda={lam}, z(lambda)={back}")
    if z in EXPECTED_LAMBDA and lam != EXPECTED_LAMBDA[z]:
        raise Refusal(f"lambda_s({z}) = {lam!r} differs from the recorded {EXPECTED_LAMBDA[z]!r}")
    lc = ar._lightcone
    public_chi = lc.chi_of_lambda(ar.COSMOLOGY, np.array([lam]), convention="project",
                                  z_max=ar.BACKGROUND_Z_MAX, distance_unit="Mpc")
    public_z = lc.z_of_chi(ar.COSMOLOGY, public_chi, z_max=ar.BACKGROUND_Z_MAX, distance_unit="Mpc")
    if float(public_chi[0]) != chi or float(public_z[0]) != back:
        raise Refusal("Cached affine maps differ from the public canoes maps at the endpoint")
    with np.load(REFINED, allow_pickle=False) as table:
        shells = np.asarray(table["lambda_shells_Mpc"], float)
    if not shells[0] <= lam <= shells[-1]:
        raise Refusal(f"lambda_s={lam} lies outside the vertex shells")
    upper = int(np.searchsorted(shells, lam))
    checks = {"reviewed_six_bitwise": exact, "max_affine_mpc": max_affine, "z_roundtrip": back,
              "z_roundtrip_abs_error": abs(back - z), "public_map_parity": "exact",
              "expected_lambda_mpc": EXPECTED_LAMBDA.get(z),
              "vertex_shell_bracket_mpc": [float(shells[upper - 1]), float(shells[upper])]}
    if validate:
        checks["aligned_response_validate"] = ar.validate()
    return {"z": z, "lambda_mpc": lam, "chi_mpc": chi, "D_mpc": dist}, checks


def load_endpoint(source_map: Path, z: float) -> dict:
    """Read a source map with the reviewed rules and re-derive the requested endpoint.

    Only two endpoints are served: z_s=5 on the reviewed map (gates) and z_s=3 on a
    map this driver wrote (production). The target tag follows from z_s.
    """
    source_map, z = source_map.resolve(), float(z)
    mapping = read_source_map_checked(source_map)
    if source_map == REVIEWED_SOURCE_MAP:
        if z != GATE_Z:
            raise Refusal(f"The reviewed source map serves only the z_s={GATE_Z:g} gate batch, not z={z}")
        tag = f"zs{z:g}".replace(".", "p")
    else:
        producer = mapping.get("producer", {}).get("driver", {})
        if producer.get("sha256") != digest(DRIVER) or not inside(source_map, HERE):
            raise Refusal("A new source map must be written by this driver, inside its directory")
        if z != PRODUCTION_Z:
            raise Refusal(f"A driver-written source map serves only z_s={PRODUCTION_Z:g}, not z={z}")
        tag = mapping.get("target_tag")
        if tag != tag_for(z):
            raise Refusal(f"Source-map target_tag {tag!r} differs from {tag_for(z)!r}")
    require_tag(tag)
    if Path(mapping["response_module"]) != R1 / "aligned_response.py" or mapping["response_attr"] != "response_R":
        raise Refusal("Source map response differs from aligned_response.response_R")
    hits = np.flatnonzero(np.isclose(mapping["z"], z, rtol=0, atol=1e-12))
    if len(hits) != 1:
        raise Refusal(f"Source map has no unique endpoint for z={z}")
    i = int(hits[0])
    entry = {k: float(mapping[k][i]) for k in ("z", "lambda_mpc", "chi_mpc", "D_mpc")}
    derived, _ = derive_endpoint(entry["z"])
    if derived != entry or entry["z"] != z:
        raise Refusal(f"Source-map entry {entry} differs from the re-derived endpoint {derived}")
    return {**entry, "t_max": float(math.ceil(entry["lambda_mpc"])), "tag": tag,
            "source_map": record(source_map)}


# -------------------------------------------------------------- accepted records
def accepted_records() -> dict:
    """Accepted T-001 admissions: main v6 execution files and final_ladder_v8 continuation targets."""
    manifest = read_json(LADDER_V8)
    main = {"id": "main_fk_gl36", "root": str(MAIN_EXEC),
            "config": yaml.safe_load((MAIN_EXEC / "config_main.yaml").read_text(encoding="utf-8")),
            "loader_text": (MAIN_EXEC / "callable_main.py").read_text(encoding="utf-8"),
            "contract": read_json(MAIN_EXEC / "main_contract.json")}
    first = manifest["admissions"][0]["target"]
    if (first["id"] != "main_fk_gl36" or first["config"] != main["config"]
            or first["loader_text"] != main["loader_text"]):
        raise Refusal("final_ladder_v8 main admission differs from the main execution files")
    ladder = {a["target"]["id"]: a["target"] for a in manifest["admissions"][1:]}
    return {"main_fk_gl36": main, **ladder}


def accepted_n_jobs(spec: Spec, gl: int, accepted: dict) -> int:
    if spec.accepted is None or gl != GL_ACCEPTED:
        raise Refusal(f"No accepted GL{gl} record for this target; pass an explicit --n-jobs")
    return int(accepted[spec.accepted]["config"]["sweep"]["n_jobs"])


# ------------------------------------------------------------------ pins
def pin_records(paths: set[Path]) -> list[dict]:
    unique = {str(p.resolve()): p for p in paths}
    return [{"path": key, "sha256": digest(unique[key])} for key in sorted(unique)]


def shared_pins() -> list[dict]:
    """Inputs common to every batch: reviewed modules, runner, driver, framework export, canoes sources."""
    paths = set(FOLD_INPUTS) | {DRIVER}
    paths |= {p for p in FW.rglob("*") if p.is_file()}
    paths |= set((CANOES / "src/canoes").rglob("*.py"))
    return pin_records(paths)


def batch_pins(names: list[str], source_map: Path, common: list[dict]) -> list[dict]:
    """Contract pins: shared inputs plus this batch's tables, templates and source map."""
    paths = {source_map}
    for name in names:
        spec = TARGETS[name]
        paths |= {spec.table, CAND / "configs" / f"{spec.template}.yaml"}
    merged = {p["path"]: p for p in common + pin_records(paths)}
    return [merged[key] for key in sorted(merged)]


def check_canoes_pins(pins: list[dict], accepted: dict) -> None:
    """Canoes Python sources (tracked and untracked) equal the accepted contract pins exactly."""
    prefix = str(CANOES / "src/canoes") + "/"
    reference = {p["path"]: p["sha256"] for p in accepted["bihalofit_auto_15360_gl36"]["contract"]["pins"]
                 if p["path"].startswith(prefix)}
    mine = {p["path"]: p["sha256"] for p in pins if p["path"].startswith(prefix)}
    if mine != reference:
        raise Refusal("canoes source files differ from the accepted T-001 contract pins")


# ---------------------------------------------------------- materialisation
def check_template(cfg: dict) -> None:
    """The T-001 template is the GL24 z_s=5 config with 11 directions and one K vertex."""
    expected = {("system", "linear", "R_time_module"): str(R1 / "aligned_response.py"),
                ("propagators", "c_closed_form_module"): str(CAND / "corr_op_C_callable.py"),
                ("propagators", "t_max"): 2318.0, ("sweep", "t_final_grid"): [REVIEWED_Z5_LAMBDA],
                ("sweep", "component_pairs"): [[0, 0]], ("sweep", "n_gauss"): TEMPLATE_GL, ("sweep", "n_jobs"): 4,
                ("expand", "n_jobs"): 1, ("propagators", "n_jobs"): 4, ("sweep", "vertex_types"): ["FK"]}
    for keys, value in expected.items():
        node: Any = cfg
        for key in keys:
            node = node[key]
        if node != value:
            raise Refusal(f"Template {'.'.join(keys)}={node!r}, expected {value!r}")
    if len(cfg["sweep"]["positions_grid"]["y"]) != N_LADDER_ROWS or len(cfg["output"]) != 1:
        raise Refusal(f"Template must have {N_LADDER_ROWS} ladder directions and one output")
    if [v["name"] for v in cfg["system"]["nonlocal_vertices"]] != ["K"]:
        raise Refusal("Template must have exactly one nonlocal K vertex")


def materialise_config(template: dict, root: Path, result: Path, source: dict, gl: int, n_jobs: int) -> dict:
    """Deep copy of the template with paths, endpoint, GL, workers and the T-001 R_time_vectorized change."""
    cfg = copy.deepcopy(template)
    cfg["system"]["nonlocal_vertices"][0]["coupling_module"] = str(root / "callable_main.py")
    cfg["system"]["linear"]["R_time_vectorized"] = True
    cfg["expand"]["cache_path"] = str(root / "cache_expand")
    cfg["propagators"]["cache_path"] = str(root / "cache_propagators")
    cfg["propagators"]["t_max"] = source["t_max"]
    cfg["sweep"]["t_final_grid"] = [source["lambda_mpc"]]
    cfg["sweep"]["n_gauss"] = gl
    cfg["sweep"]["n_jobs"] = n_jobs
    cfg["sweep"]["component_pairs"] = [[0, 0]]
    cfg["output"][0]["path"] = str(result)
    return cfg


def table_info(spec: Spec) -> dict:
    """continuation_guard.validate_table plus geometry, refinement overlay and manifest checks."""
    scope = "multiz" if spec.scope == "main" else "ladder"
    try:
        info = continuation_guard().validate_table(spec.table, scope, spec.model, spec.cutoff, spec.leg)
    except (ValueError, KeyError, OSError) as exc:
        raise Refusal(f"Table {spec.table} fails continuation_guard.validate_table: {exc}") from exc
    # The refined main table carries its overlay; the 34-row ladder tables are unrefined.
    geometry = MAIN_GEOMETRY if spec.scope == "main" else LADDER_GEOMETRY
    if info["geometry_sha256"] != geometry or info["refinement_overlay_present"] != (spec.scope == "main"):
        raise Refusal(f"Unexpected geometry or refinement overlay in {spec.table}")
    if spec.scope != "main":
        manifest = {t["path"]: t["sha256"] for t in read_json(PIECES / "ladder_manifest.json")["tables"]}
        if manifest.get(str(spec.table)) != digest(spec.table):
            raise Refusal(f"Ladder table is missing from or stale in ladder_manifest.json: {spec.table}")
    return info


def make_contract(spec: Spec, info: dict, directions: list, t_max: float, pins: list[dict]) -> dict:
    contract = {"pins": pins, "table": str(spec.table), "table_sha256": digest(spec.table),
                "shape": info["shape"], "shells": info["shells"], "geometry_sha256": info["geometry_sha256"],
                "directions": directions, "query_interval": [0.0, t_max]}
    if spec.scope != "main":
        contract.update(scope="ladder", model=spec.model, cutoff=spec.cutoff, leg_order=spec.leg)
    return contract


def make_loader(spec: Spec, contract_path: Path, contract_sha: str) -> str:
    """Same construction as the accepted T-001 loaders: PCHIP source, TABLE_PATH, guard install."""
    source = PCHIP.read_text(encoding="utf-8")
    if source.count(MARKER) != 1:
        raise Refusal("Ambiguous TABLE_PATH marker in the PCHIP source")
    text = source.replace(MARKER, f"TABLE_PATH = Path({str(spec.table)!r})")
    if spec.scope == "main":
        text += "\nimport importlib.util as _guard_import\n"
        text += (f'_guard_spec = _guard_import.spec_from_file_location("candidate_main_guard", '
                 f"{str(REVIEW / 'main_guard.py')!r})\n")
        text += "_guard_mod = _guard_import.module_from_spec(_guard_spec)\n_guard_spec.loader.exec_module(_guard_mod)\n"
        text += f"_guard_mod.install(globals(), {str(contract_path)!r}, {contract_sha!r})\n"
    else:
        text += "\nimport importlib.util as _gi\n"
        text += (f'_gs = _gi.spec_from_file_location("candidate_continuation_guard", '
                 f"{str(CONT / 'continuation_guard.py')!r})\n")
        text += "_gm = _gi.module_from_spec(_gs)\n_gs.loader.exec_module(_gm)\n"
        text += f"_gm.install(globals(), {str(contract_path)!r}, {contract_sha!r})\n"
    return text


# ----------------------------------------------------- diff against accepted
CONFIG_CATEGORIES = {
    "expand.cache_path": "path", "propagators.cache_path": "path", "output.0.path": "path",
    "system.nonlocal_vertices.0.coupling_module": "path", "sweep.t_final_grid": "endpoint",
    "propagators.t_max": "endpoint", "sweep.positions_grid.y": "main_subset",
    "sweep.component_pairs": "main_subset", "sweep.n_jobs": "workers", "sweep.n_gauss": "quadrature",
}
CONTRACT_CATEGORIES = {"query_interval": "endpoint", "directions": "main_subset"}


def flatten(value: Any, prefix: tuple = ()) -> Iterator[tuple[str, Any]]:
    """Leaf fields of nested dicts (and lists of dicts) as dotted paths."""
    if isinstance(value, dict):
        for key in sorted(value):
            yield from flatten(value[key], prefix + (str(key),))
    elif isinstance(value, list) and value and all(isinstance(v, dict) for v in value):
        for index, item in enumerate(value):
            yield from flatten(item, prefix + (str(index),))
    else:
        yield ".".join(prefix), value


def leaf_diff(mine: dict, accepted: dict, categories: dict, allowed: set) -> list[dict]:
    left, right = dict(flatten(mine)), dict(flatten(accepted))
    rows = []
    for field in sorted(set(left) | set(right)):
        a, b = left.get(field, "<absent>"), right.get(field, "<absent>")
        if a != b:
            category = categories.get(field, "unexpected")
            rows.append({"field": field, "materialised": a, "accepted": b, "category": category,
                         "allowed": category in allowed})
    return rows


def pin_diff(mine: list[dict], accepted: list[dict]) -> dict:
    left = {p["path"]: p["sha256"] for p in mine}
    right = {p["path"]: p["sha256"] for p in accepted}
    changed = sorted(p for p in set(left) & set(right) if left[p] != right[p])
    return {"materialised_count": len(left), "accepted_count": len(right), "common": len(set(left) & set(right)),
            "only_materialised": sorted(set(left) - set(right)), "only_accepted": sorted(set(right) - set(left)),
            "changed_common": changed}


def loader_diff(mine: str, mine_contract: tuple[str, str], accepted: str, accepted_contract: tuple[str, str]) -> dict:
    def normal(text: str, path: str, sha: str) -> str:
        return text.replace(repr(path), "'<CONTRACT_PATH>'").replace(repr(sha), "'<CONTRACT_SHA256>'")
    lines = list(difflib.unified_diff(accepted.splitlines(), mine.splitlines(), "accepted", "materialised",
                                      n=0, lineterm=""))
    return {"equal_after_contract_binding": normal(mine, *mine_contract) == normal(accepted, *accepted_contract),
            "raw_diff": lines}


def diff_against_accepted(target: dict, spec: Spec, accepted: dict, source: dict, gl: int) -> dict | None:
    """Field-by-field diff with the accepted rung; any difference outside the allowed categories refuses."""
    if spec.accepted is None:
        return None
    ref = accepted[spec.accepted]
    allowed = {"path"}
    if source["lambda_mpc"] != REVIEWED_Z5_LAMBDA:
        allowed.add("endpoint")
    if spec.scope == "main":
        allowed.add("main_subset")
    if target["n_jobs_source"] == "cli":
        allowed.add("workers")
    if gl != GL_ACCEPTED:
        allowed.add("quadrature")
    ref_contract = ref["contract"]
    ref_path = str(Path(ref["root"]) / "main_contract.json")
    contract = {k: v for k, v in target["contract"].items() if k != "pins"}
    report = {
        "accepted_id": spec.accepted, "allowed_categories": sorted(allowed),
        "config": leaf_diff(target["config"], ref["config"], CONFIG_CATEGORIES, allowed),
        "contract": leaf_diff(contract, {k: v for k, v in ref_contract.items() if k != "pins"},
                              CONTRACT_CATEGORIES, allowed),
        "pins": pin_diff(target["contract"]["pins"], ref_contract["pins"]),
        "loader": loader_diff(target["loader_text"], (str(Path(target["root"]) / "main_contract.json"),
                                                      target["contract_sha256"]),
                              ref["loader_text"], (ref_path, sha_bytes(canonical(ref_contract)))),
    }
    unexpected = [row["field"] for row in report["config"] + report["contract"] if not row["allowed"]]
    if report["pins"]["changed_common"]:
        unexpected.append("pins changed: " + ", ".join(report["pins"]["changed_common"]))
    if not report["loader"]["equal_after_contract_binding"]:
        unexpected.append("loader text")
    if unexpected:
        raise Refusal(f"{target['id']} differs from accepted {spec.accepted} beyond the allowed fields: {unexpected}")
    return report


# ------------------------------------------------------------------ batch plan
def build_batch(args: dict) -> dict:
    """Deterministic in-memory batch plan; reads and hashes, never writes."""
    verify_pinned()
    runtime, framework, canoes = runtime_state(), framework_state(), canoes_state()
    source = load_endpoint(Path(args["source_map"]), args["z"])
    accepted = accepted_records()
    common = shared_pins()
    check_canoes_pins(common, accepted)
    pins = batch_pins(args["targets"], Path(args["source_map"]), common)
    work = Path(args["work_dir"])
    main_y = accepted["main_fk_gl36"]["config"]["sweep"]["positions_grid"]["y"]
    targets = [build_target(name, work, source, args, accepted, pins, main_y) for name in args["targets"]]
    return {"schema": 1, "kind": "zs3p0 T-001 ladder replay batch", "scope": SCOPE,
            "status": "CANDIDATE_NOT_SCIENTIFICALLY_ACCEPTED", "args": args, "driver": record(DRIVER),
            "stage_test": bool(args["stage_test"]), "runtime": runtime, "framework": framework, "canoes": canoes,
            "source": source, "reference_pins": [record(p) for p in REFERENCE_FILES], "shared_pins": common,
            "pins": pins, "targets": targets}


def build_target(name: str, work: Path, source: dict, args: dict, accepted: dict, pins: list[dict],
                 main_y: list) -> dict:
    """In-memory config, contract and loader text of one target, diffed against its accepted rung."""
    spec = TARGETS[name]
    ident = f"{name}_{require_tag(source['tag'])}_gl{args['gl']}"
    root, template_path = work / ident, CAND / "configs" / f"{spec.template}.yaml"
    result = root / f"xi_{ident}.npz"
    template = yaml.safe_load(template_path.read_text(encoding="utf-8"))
    check_template(template)
    if spec.scope == "main" and template["sweep"]["positions_grid"]["y"] != [main_y[i] for i in MAIN_SUBSET]:
        raise Refusal("Ladder directions are not the accepted main subset [0, 5..14]")
    if args["n_jobs"] == "accepted":
        n_jobs, n_source = accepted_n_jobs(spec, args["gl"], accepted), f"accepted:{spec.accepted}"
    else:
        n_jobs, n_source = int(args["n_jobs"]), "cli"
    cfg = materialise_config(template, root, result, source, args["gl"], n_jobs)
    config_text = canonical(cfg).decode() + "\n"
    if yaml.safe_load(config_text) != cfg:
        raise Refusal("Canonical JSON config does not round-trip through YAML")
    info = table_info(spec)
    directions = [cfg["sweep"]["positions_grid"]["x"][0], *cfg["sweep"]["positions_grid"]["y"]]
    contract = make_contract(spec, info, directions, source["t_max"], pins)
    contract_sha = sha_bytes(canonical(contract))
    loader_text = make_loader(spec, root / "main_contract.json", contract_sha)
    target = {"id": ident, "name": name, **{k: v for k, v in dataclasses.asdict(spec).items()
                                            if k not in ("table", "accepted", "member")},
              "guard": "main_guard" if spec.scope == "main" else "continuation_guard",
              "table": str(spec.table), "table_sha256": digest(spec.table),
              "table_info": {k: info[k] for k in ("shape", "geometry_sha256", "band_quadrature",
                                                  "refinement_overlay_present", "native_metadata")},
              "template_config": record(template_path), "root": str(root), "result": str(result),
              "config": cfg, "config_sha256": sha_bytes(config_text.encode()), "loader_text": loader_text,
              "loader_sha256": sha_bytes(loader_text.encode()), "contract": contract,
              "contract_sha256": contract_sha, "n_jobs": n_jobs, "n_jobs_source": n_source, "gl": args["gl"],
              "expected_rows": N_LADDER_ROWS, "member": spec.member, "accepted_analog": spec.accepted}
    target["accepted_diff"] = diff_against_accepted(target, spec, accepted, source, args["gl"])
    return target


def generated_records(target: dict) -> list[dict]:
    root = Path(target["root"])
    return [{"path": str(root / "main_contract.json"), "sha256": target["contract_sha256"]},
            {"path": str(root / "callable_main.py"), "sha256": target["loader_sha256"]},
            {"path": str(root / "config_main.yaml"), "sha256": target["config_sha256"]}]


def loader_self_check(target: dict) -> str:
    """Import the loader in the fold environment: its guard must install (no vertex query)."""
    root = Path(target["root"])
    code = ("import importlib.util as u, sys\n"
            "s = u.spec_from_file_location('zs3_loader_self_check', sys.argv[1])\n"
            "m = u.module_from_spec(s)\ns.loader.exec_module(m)\n"
            "assert m.ALLOW_INTERPOLATION is False and m._CACHE is None\n"
            "print('GUARD_INSTALLED', m.coupling_fn_batch.__name__, str(m.TABLE_PATH))\n")
    done = subprocess.run([str(PYTHON), "-B", "-c", code, str(root / "callable_main.py")],
                          env=fold_environment(root), cwd=root, capture_output=True, text=True, check=False)
    expected = f"GUARD_INSTALLED batch {target['table']}"
    if done.returncode != 0 or done.stdout.strip() != expected:
        raise Refusal(f"Loader guard self-check failed for {target['id']}: {done.stdout}{done.stderr[-2000:]}")
    return done.stdout.strip()


def refuse_production_location(target: dict, work: Path) -> None:
    """No output beside a template config or under a PRODUCTION-marked folder; checked before mkdir."""
    result, template = Path(target["result"]), Path(target["template_config"]["path"])
    if not inside(result, work) or result.parent.parent != work:
        raise Refusal(f"Output {result} is not in its own target root inside {work}")
    markers = [folder / "PRODUCTION" for folder in (result.parent, work, HERE)]
    if result.parent.resolve() == template.parent.resolve() or any(m.exists() for m in markers):
        raise Refusal(f"Output would sit in a production or template directory: {result}")


# ---------------------------------------------------------------- prepare
def n_jobs_arg(text: str) -> int | str:
    """argparse type: a positive ASCII integer, or the literal 'accepted'."""
    if text == "accepted":
        return text
    if not re.fullmatch(r"[1-9][0-9]*", text):
        raise argparse.ArgumentTypeError("must be a positive integer or 'accepted'")
    return int(text)


def gl_arg(text: str) -> int:
    if not re.fullmatch(r"[1-9][0-9]*", text):
        raise argparse.ArgumentTypeError("must be a positive integer")
    return int(text)


def materialise_batch(work: Path, plan: dict) -> dict:
    """Write each target's contract, loader and config; check them and the loader guard in the fold env."""
    checks = {}
    for target in plan["targets"]:
        root = Path(target["root"])
        root.mkdir(exist_ok=False)
        write_new_bytes(root / "main_contract.json", canonical(target["contract"]))
        write_new_text(root / "callable_main.py", target["loader_text"])
        write_new_text(root / "config_main.yaml", canonical(target["config"]).decode() + "\n")
        check_records(generated_records(target))
        written = yaml.safe_load((root / "config_main.yaml").read_text(encoding="utf-8"))
        try:
            replay().validate_config_paths(written, work)
        except (RuntimeError, OSError) as exc:
            raise Refusal(f"{target['id']} config fails replay_suite.validate_config_paths: {exc}") from exc
        subprocess_binding(fold_environment(root), root)
        checks[target["id"]] = loader_self_check(target)
    return checks


def record_prepare_failure(work: Path, exc: BaseException) -> None:
    """Mark a batch directory whose preparation failed; it never becomes runnable."""
    note = {"schema": 1, "status": "PREPARE_FAILED", "error": repr(exc),
            "note": "Not a batch: there is no prepared.json. Prepare again with a new --work-dir."}
    try:
        write_new_json(work / "prepare_failed.json", note)
    except (OSError, Refusal) as marker_error:
        print(f"WARNING: could not write {work / 'prepare_failed.json'}: {marker_error}", file=sys.stderr)


def prepare_command(ns: argparse.Namespace) -> int:
    """Build the plan in memory, then materialise it into a new --work-dir and write prepared.json."""
    require_invocation()
    work = Path(os.path.abspath(ns.work_dir)).resolve()
    if HERE != PLANNED_HOME and not ns.allow_stage_test:
        raise Refusal(f"The driver must sit at {PLANNED_HOME} (use --allow-stage-test only for test prepares)")
    if not inside(work, HERE) or work.exists() or not work.parent.is_dir():
        raise Refusal(f"--work-dir must be a new directory inside {HERE} with an existing parent: {work}")
    args = {"source_map": str(Path(ns.source_map).resolve()), "z": float(ns.z), "gl": int(ns.gl),
            "n_jobs": ns.n_jobs, "targets": list(dict.fromkeys(ns.targets)), "work_dir": str(work),
            "stage_test": bool(ns.allow_stage_test)}
    plan = build_batch(args)
    for target in plan["targets"]:
        refuse_production_location(target, work)
    work.mkdir(exist_ok=False)
    try:
        checks = materialise_batch(work, plan)
        prepared = {"schema": 1, "status": "PREPARED_NOT_STARTED", "scope": SCOPE, "prepared_unix": time.time(),
                    "plan_sha256": sha_bytes(canonical(plan)), "loader_self_check": checks,
                    "generated": {t["id"]: generated_records(t) for t in plan["targets"]}, "plan": plan}
        write_new_json(work / "prepared.json", prepared)
    except BaseException as exc:
        record_prepare_failure(work, exc)
        raise
    print(json.dumps({"prepared": str(work / "prepared.json"), "plan_sha256": prepared["plan_sha256"],
                      "stage_test": plan["stage_test"], "source": {k: plan["source"][k] for k in
                                                                   ("z", "lambda_mpc", "t_max", "tag")},
                      "targets": [{"id": t["id"], "n_jobs": t["n_jobs"], "n_jobs_source": t["n_jobs_source"],
                                   "contract_sha256": t["contract_sha256"], "pins": len(t["contract"]["pins"])}
                                  for t in plan["targets"]]}, indent=2))
    return 0


def load_prepared(work: Path) -> tuple[dict, dict]:
    """prepared.json of a runnable batch: intact plan hash, this driver, its own --work-dir, not a stage test."""
    work = Path(os.path.abspath(work)).resolve()
    path = work / "prepared.json"
    if not path.is_file():
        raise Refusal(f"No prepared.json in {work}; a failed or missing prepare needs a new --work-dir")
    prepared = read_json(path)
    plan = prepared["plan"]
    if sha_bytes(canonical(plan)) != prepared["plan_sha256"]:
        raise Refusal("prepared.json plan hash mismatch")
    if plan["stage_test"] or HERE != PLANNED_HOME:
        raise Refusal("Stage-test batches, or a driver outside its planned directory, cannot run, gate or assemble")
    if plan["driver"]["sha256"] != digest(DRIVER):
        raise Refusal("The driver changed after preparation")
    if Path(plan["args"]["work_dir"]) != work:
        raise Refusal(f"prepared.json was prepared for {plan['args']['work_dir']}, not {work}")
    return prepared, plan


def require_current(plan: dict, prepared: dict) -> None:
    """Rebuild the plan from current inputs; any difference refuses."""
    if canonical(build_batch(plan["args"])) != canonical(plan):
        raise Refusal("Inputs, runtime or records changed since preparation")
    for records in prepared["generated"].values():
        check_records(records)


def is_gate_plan(plan: dict) -> bool:
    """A z_s=5 batch on the reviewed source map."""
    source = plan["source"]
    return (source["source_map"]["path"] == str(REVIEWED_SOURCE_MAP) and source["z"] == GATE_Z
            and source["lambda_mpc"] == REVIEWED_Z5_LAMBDA)


def is_production_plan(plan: dict) -> bool:
    """A z_s=3.0 batch on a source map this driver wrote."""
    source = plan["source"]
    return (source["source_map"]["path"] != str(REVIEWED_SOURCE_MAP) and source["z"] == PRODUCTION_Z
            and source["lambda_mpc"] == EXPECTED_LAMBDA[PRODUCTION_Z] and source["tag"] == tag_for(PRODUCTION_Z))


def gate_report_path(ns: argparse.Namespace) -> Path:
    return Path(ns.gate_report) if ns.gate_report is not None else HERE / GATE_WORK_NAME / "gate_report.json"


# -------------------------------------------------------------------- run
def parse_vm_stat(text: str) -> int:
    """Free plus inactive bytes from vm_stat output."""
    fields = {}
    for key, pattern in (("page", r"page size of (\d+) bytes"), ("free", r"Pages free:\s+(\d+)"),
                         ("inactive", r"Pages inactive:\s+(\d+)")):
        match = re.search(pattern, text)
        if match is None:
            raise Refusal(f"Cannot parse vm_stat output ({key})")
        fields[key] = int(match.group(1))
    return (fields["free"] + fields["inactive"]) * fields["page"]


def available_memory_bytes() -> int:
    try:
        done = subprocess.run(["vm_stat"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise Refusal(f"vm_stat failed: {exc}") from exc
    return parse_vm_stat(done.stdout)


def memory_admission(n_jobs: int, min_gib: float | None) -> dict:
    """Free + inactive memory must cover the measured T-001 peak for this worker count plus 12 GiB."""
    if min_gib is None:
        if n_jobs not in MEASURED_PEAK_BYTES:
            raise Refusal(f"No measured peak for {n_jobs} workers; pass --min-available-gib")
        required = MEASURED_PEAK_BYTES[n_jobs] + MEMORY_FLOOR_BYTES
    else:
        required = int(min_gib * 2**30)
    available = available_memory_bytes()
    if available < required:
        raise Refusal(f"Available memory {available / 2**30:.1f} GiB < required {required / 2**30:.1f} GiB")
    return {"method": "vm_stat free+inactive, observation not reservation", "available_bytes": available,
            "required_bytes": required, "rule": "explicit" if min_gib is not None else "measured peak + 12 GiB"}


def refuse_existing(target: dict) -> None:
    root, result = Path(target["root"]), Path(target["result"])
    names = ("started.json", "completed.json", "failed.json", "cache_expand", "cache_propagators")
    if any((root / n).exists() for n in names) or result.exists():
        raise Refusal(f"Refusing an existing or interrupted run in {root}; prepare a fresh batch")


@contextlib.contextmanager
def shared_fold_lock() -> Iterator[None]:
    """The reviewed fold_lock on the shared r1_sft061 lock; a held lock or an I/O error refuses."""
    with contextlib.ExitStack() as stack:
        try:
            stack.enter_context(shared().fold_lock(LOCK))
        except (RuntimeError, OSError) as exc:
            raise Refusal(f"Cannot take the shared fold lock {LOCK}: {exc}") from exc
        yield


def require_no_active_folds() -> None:
    """The reviewed ps check; an active fold, or ps being unavailable (sandbox), refuses."""
    try:
        shared().check_active_folds()
    except (RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        raise Refusal(f"Active-fold check refused: {exc}") from exc


def run_command(ns: argparse.Namespace) -> int:
    """One fold. A z_s=3 batch first needs a PASS gate report bound to the same inputs."""
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    target = next((t for t in plan["targets"] if t["id"] == ns.target), None)
    if target is None:
        raise Refusal(f"Unknown target {ns.target}")
    refuse_existing(target)
    gate = None if is_gate_plan(plan) else gate_evidence(gate_report_path(ns), plan)
    with shared_fold_lock():
        require_no_active_folds()
        require_current(plan, prepared)
        refuse_existing(target)
        memory = memory_admission(target["n_jobs"], ns.min_available_gib)
        root = Path(target["root"])
        env = fold_environment(root)
        subprocess_binding(env, root)
        loader_self_check(target)
        execute_fold(work, prepared, plan, target, env, memory, gate)
    return 0


def execute_fold(work: Path, prepared: dict, plan: dict, target: dict, env: dict, memory: dict,
                 gate: dict | None = None) -> None:
    """started.json, run_fk_variant.run in the fold environment, result validation, completed.json."""
    root, result, config = Path(target["root"]), Path(target["result"]), Path(target["root"]) / "config_main.yaml"
    fold = variant().run
    write_new_json(root / "started.json", {
        "schema": 1, "status": "STARTED_NOT_COMPLETED", "target_id": target["id"], "pid": os.getpid(),
        "start_unix": time.time(), "plan_sha256": prepared["plan_sha256"], "prepared": record(work / "prepared.json"),
        "lock": str(LOCK), "n_jobs": target["n_jobs"], "gl": target["gl"], "orders": ORDERS, "memory": memory,
        "gate": gate, "environment": {k: env[k] for k in ("PYTHONPATH", "CANOES_ROOT", "PYTHONNOUSERSITE",
                                                           "MPLCONFIGDIR", *THREAD_VARS)}})
    caller_env = dict(os.environ)
    started = time.monotonic()
    try:
        # run_fk_variant.run passes os.environ to run_FF_single.py; it takes no env argument.
        os.environ.clear()
        os.environ.update(env)
        fold(config, result, int(target["gl"]), ORDERS)
        check_records(plan["pins"] + plan["reference_pins"] + prepared["generated"][target["id"]])
        if not result.is_file():
            raise Refusal("Fold finished without a result file")
        validation = validate_fold_result(result, target)
    except BaseException as exc:
        write_new_json(root / "failed.json", {"schema": 1, "status": "failed", "target_id": target["id"],
                                              "error": repr(exc), "elapsed_seconds": time.monotonic() - started})
        raise
    finally:
        os.environ.clear()
        os.environ.update(caller_env)
    write_new_json(root / "completed.json", {
        "schema": 1, "status": "completed", "target_id": target["id"], "plan_sha256": prepared["plan_sha256"],
        "prepared": record(work / "prepared.json"), "started": record(root / "started.json"),
        "elapsed_seconds": time.monotonic() - started, "returncode": 0, "result": record(result),
        "output_validation": validation, "scientific_acceptance": False,
        "max_single_process_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        "rss_note": ("getrusage(RUSAGE_CHILDREN).ru_maxrss: the largest single reaped descendant, in bytes "
                     "on macOS; not the summed process-group RSS quoted in the T-001 memory records")})
    print(json.dumps({"completed": target["id"], "result": record(result)}, indent=2))


# --------------------------------------------------------- result handling
def read_rows(path: Path) -> dict:
    """Numeric members without pickle; raw object geometry through the reviewed restricted reader."""
    with np.load(path, allow_pickle=False) as data:
        rows = {k: np.asarray(data[k]) for k in ("t_final", "a", "b", "order", "value")}
    count = rows["value"].shape[0]
    with zipfile.ZipFile(path) as archive:
        for key in ("x", "y"):
            header = archive.read(key + ".npy")[:128]
            if b"'O'" in header or b"|O" in header:
                rows[key] = shared().read_geometry(archive, key, count)
            else:
                with np.load(path, allow_pickle=False) as data:
                    rows[key] = np.asarray(data[key], dtype=float)
    return rows


def validate_fold_result(result: Path, target: dict) -> dict:
    """Reviewed validate_result plus the 11-row (0,0) order-2 endpoint and direction check."""
    try:
        schema = shared().validate_result(result, target["config"])
        rows = read_rows(result)
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        raise Refusal(f"Result {result} fails the reviewed result validation: {exc}") from exc
    lam = target["config"]["sweep"]["t_final_grid"][0]
    directions = np.asarray(target["contract"]["directions"], float)
    ok = (rows["value"].shape == (N_LADDER_ROWS,) and np.all(rows["a"] == 0) and np.all(rows["b"] == 0)
          and np.all(rows["order"] == 2) and np.all(rows["t_final"] == lam)
          and np.all(rows["x"] == directions[0])
          and sorted(map(tuple, rows["y"])) == sorted(map(tuple, directions[1:])))
    if not ok:
        raise Refusal(f"Result {result} fails the {N_LADDER_ROWS}-row (0,0) order-2 endpoint check")
    return {**schema, "rows": N_LADDER_ROWS, "t_final": lam, "pair": [0, 0],
            "value_sha256": sha_bytes(np.ascontiguousarray(rows["value"]).tobytes())}


def completed_rows(target: dict, plan_sha: str) -> tuple[dict, dict]:
    """Rows of a completed fold whose completion, started record and result are all intact."""
    root = Path(target["root"])
    if (root / "failed.json").exists():
        raise Refusal(f"{target['id']} has failed.json")
    if not (root / "completed.json").is_file():
        raise Refusal(f"{target['id']} has not completed (no completed.json in {root})")
    completed = read_json(root / "completed.json")
    if completed["target_id"] != target["id"] or completed["plan_sha256"] != plan_sha:
        raise Refusal(f"{target['id']} completion belongs to another plan")
    check_records([completed["result"], completed["started"]])
    if completed["result"]["path"] != str(Path(target["result"]).resolve()):
        raise Refusal(f"{target['id']} completion points at another result")
    if validate_fold_result(Path(target["result"]), target) != completed["output_validation"]:
        raise Refusal(f"{target['id']} result no longer matches its completion record")
    return read_rows(Path(target["result"])), {"completed": record(root / "completed.json"),
                                               "result": completed["result"],
                                               "elapsed_seconds": completed["elapsed_seconds"]}


def row_keys(rows: dict) -> list[tuple]:
    return [(tuple(rows["x"][i]), tuple(rows["y"][i]), float(rows["t_final"][i]), int(rows["a"][i]),
             int(rows["b"][i]), int(rows["order"][i])) for i in range(rows["value"].shape[0])]


def compare_rows(mine: dict, ref: dict, label: str, reference: Path) -> dict:
    """Row-matched comparison: max relative difference, bitwise flag, pass at REL_TOL."""
    keys, ref_keys = row_keys(mine), row_keys(ref)
    if len(set(keys)) != len(keys) or sorted(keys) != sorted(ref_keys):
        return {"label": label, "reference": record(reference), "pass": False, "reason": "row sets differ"}
    order = [ref_keys.index(k) for k in keys]
    value, ref_value = np.asarray(mine["value"], float), np.asarray(ref["value"], float)[order]
    abs_diff = np.abs(value - ref_value)
    scale = np.abs(ref_value)
    rel = np.where(scale > 0, abs_diff / np.where(scale > 0, scale, 1.0), np.where(abs_diff > 0, np.inf, 0.0))
    max_rel = float(np.max(rel))
    return {"label": label, "reference": record(reference), "rows": len(keys),
            "same_row_order": order == list(range(len(keys))), "max_abs_diff": float(np.max(abs_diff)),
            "max_rel_diff": max_rel, "bitwise": value.tobytes() == ref_value.tobytes(),
            "tolerance_rel": REL_TOL, "pass": bool(max_rel <= REL_TOL)}


def main_subset_rows() -> dict:
    """Accepted main FK pair (0,0) at main direction indices [0, 5..14]."""
    path = ACCEPTED / "products/main_fk.npz"
    with np.load(path, allow_pickle=False) as data:
        sel = (data["a"] == 0) & (data["b"] == 0) & (data["order"] == 2)
        rows = {k: np.asarray(data[k])[sel] for k in ("x", "y", "t_final", "a", "b", "order", "value")}
    grid = read_json(MAIN_EXEC / "config_main.yaml")["sweep"]["positions_grid"]["y"]
    if rows["value"].shape != (40,) or not np.array_equal(rows["y"], np.asarray(grid, float)):
        raise Refusal("Accepted main FK does not hold the 40 (0,0) rows in main-grid order")
    return {k: v[MAIN_SUBSET] for k, v in rows.items()}


def reference_rows(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {k: np.asarray(data[k]) for k in ("x", "y", "t_final", "a", "b", "order", "value")}


# -------------------------------------------------------------------- gate
def evaluate_gates(work: Path, prepared: dict, plan: dict) -> dict:
    """G1/G2 report of a z_s=5 gate batch from its completed folds; reads only."""
    if not is_gate_plan(plan):
        raise Refusal("Gates must be a z_s=5 batch on the reviewed source map")
    by_name = {t["name"]: t for t in plan["targets"]}
    if any(name not in by_name or by_name[name]["gl"] != GL_ACCEPTED for name in GATES):
        raise Refusal(f"Gate batch must contain GL{GL_ACCEPTED} {GATES}")
    bihalofit, bh_record = completed_rows(by_name["bihalofit_auto_15360"], prepared["plan_sha256"])
    tree, tree_record = completed_rows(by_name["tree_main_15360"], prepared["plan_sha256"])
    ladder_bh = ACCEPTED / "cutoff_ladder/bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz"
    ladder_tree = ACCEPTED / "cutoff_ladder/tree_auto/table_tree_cut15360_r4_xi.npz"
    checks = {"G1_bihalofit_auto_15360": [compare_rows(bihalofit, reference_rows(ladder_bh),
                                                       "accepted bihalofit_auto l_max=15360 rung", ladder_bh)],
              "G2_tree_main_15360": [compare_rows(tree, main_subset_rows(),
                                                  "accepted main FK (0,0) at [0,5..14]",
                                                  ACCEPTED / "products/main_fk.npz"),
                                     compare_rows(tree, reference_rows(ladder_tree),
                                                  "accepted tree_auto l_max=15360 rung", ladder_tree)]}
    passed = all(c["pass"] for group in checks.values() for c in group)
    workers = {name: {"used": by_name[name]["n_jobs"], "source": by_name[name]["n_jobs_source"]} for name in GATES}
    return {"schema": 1, "status": "PASS" if passed else "FAIL", "criterion": "max |value/reference - 1| <= 1e-12",
            "bitwise_all": all(c.get("bitwise", False) for group in checks.values() for c in group),
            "prepared": record(work / "prepared.json"), "plan_sha256": prepared["plan_sha256"],
            "completions": {"bihalofit_auto_15360": bh_record, "tree_main_15360": tree_record},
            "workers": workers, "checks": checks, "scientific_acceptance": False}


def gate_command(ns: argparse.Namespace) -> int:
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    out = work / "gate_report.json"
    if out.exists() or out.is_symlink():
        raise Refusal(f"{out} exists; a gate report is never replaced (a failed gate needs a new gate batch)")
    require_current(plan, prepared)
    report = evaluate_gates(work, prepared, plan)
    write_new_json(out, report)
    print(json.dumps({k: report[k] for k in ("status", "bitwise_all", "workers")}
                     | {"checks": {k: [{x: c.get(x) for x in ("label", "max_rel_diff", "bitwise", "pass")}
                                       for c in v] for k, v in report["checks"].items()}}, indent=2))
    return 0 if report["status"] == "PASS" else 1


def gate_evidence(gate_report: Path, plan: dict) -> dict:
    """Re-verify a PASS gate report from its own batch and bind it to ``plan``'s inputs.

    The gate batch must be current, its report must equal a recomputation from its
    completed folds, and it must share the driver, runtime, framework, canoes state
    and shared pins with ``plan``. Every GL36 target of ``plan`` that has a gate
    must use the gate's table, template, GL and worker count.
    """
    gate_report = Path(os.path.abspath(gate_report)).resolve()
    if not gate_report.is_file():
        raise Refusal(f"No gate report at {gate_report}; run and gate the z_s=5 batch first")
    report = read_json(gate_report)
    if not isinstance(report, dict) or report.get("status") != "PASS":
        raise Refusal(f"Gate report {gate_report} is not PASS")
    prepared_record = report.get("prepared")
    if not isinstance(prepared_record, dict) or not isinstance(prepared_record.get("path"), str):
        raise Refusal(f"Gate report {gate_report} has no prepared-batch record")
    gate_work = Path(prepared_record["path"]).parent
    if gate_report.parent != gate_work:
        raise Refusal(f"Gate report {gate_report} does not sit in the batch it reports on ({gate_work})")
    gate_prepared, gate_plan = load_prepared(gate_work)
    if (report["prepared"] != record(gate_work / "prepared.json")
            or report["plan_sha256"] != gate_prepared["plan_sha256"]):
        raise Refusal("Gate report belongs to another prepared batch")
    if not is_gate_plan(gate_plan):
        raise Refusal("Gate batch is not a z_s=5 batch on the reviewed source map")
    require_current(gate_plan, gate_prepared)
    keys = ("driver", "runtime", "framework", "canoes", "reference_pins", "shared_pins")
    if any(gate_plan[k] != plan[k] for k in keys):
        raise Refusal("Gate batch ran with another driver, runtime, framework, canoes state or shared inputs")
    if json.loads(json.dumps(evaluate_gates(gate_work, gate_prepared, gate_plan))) != report:
        raise Refusal(f"Gate report {gate_report} differs from a recomputation from its batch")
    gates = {t["name"]: t for t in gate_plan["targets"]}
    for target in plan["targets"]:
        if target["name"] in GATES and target["gl"] == GL_ACCEPTED:
            ref = gates[target["name"]]
            if any(ref[k] != target[k] for k in GATE_MATCH_KEYS):
                raise Refusal(f"{target['id']} differs from gate target {ref['id']} in {GATE_MATCH_KEYS}")
    return {"gate_report": record(gate_report), "status": report["status"], "bitwise_all": report["bitwise_all"],
            "workers": report["workers"], "gate_plan_sha256": gate_prepared["plan_sha256"]}


# ---------------------------------------------------------------- assemble
def separation_arcmin(x: np.ndarray, y: np.ndarray) -> float:
    cosang = np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y))
    return float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))) * 60.0)


def plausibility(tree: dict, lam: float) -> dict:
    """Tree FK at z_s next to the accepted multiz z=2.5 and z=3.2 planes; no acceptance threshold."""
    with np.load(ACCEPTED / "products/multiz.npz", allow_pickle=False) as data:
        multiz = {k: np.asarray(data[k]) for k in ("z", "lam", "fk", "gamma")}
    lo, hi = int(np.flatnonzero(multiz["z"] == 2.5)[0]), int(np.flatnonzero(multiz["z"] == 3.2)[0])
    main = main_subset_rows()
    columns = [MAIN_SUBSET[[tuple(v) for v in main["y"]].index(tuple(y))] for y in tree["y"]]
    rows = []
    for i, column in enumerate(columns):
        theta = separation_arcmin(tree["x"][i], tree["y"][i])
        if abs(theta - float(multiz["gamma"][column])) > MULTIZ_GAMMA_ATOL_ARCMIN:
            raise Refusal(f"multiz column {column} does not match direction {i} ({theta} arcmin)")
        fk, low, high = float(tree["value"][i]), float(multiz["fk"][lo, column]), float(multiz["fk"][hi, column])
        rows.append({"theta_arcmin": theta,
                     "gamma_nominal_arcmin": float(multiz["gamma"][column]), "main_column": column,
                     "fk_multiz_z2p5": low, "fk_candidate": fk, "fk_multiz_z3p2": high,
                     "ratio_to_z2p5": fk / low, "ratio_to_z3p2": fk / high,
                     "fraction_between_planes": (fk - low) / (high - low)})
    rows.sort(key=lambda r: r["theta_arcmin"])
    lam_lo, lam_hi = float(multiz["lam"][lo]), float(multiz["lam"][hi])
    return {"acceptance_threshold": None, "multiz": record(ACCEPTED / "products/multiz.npz"),
            "lambda_mpc": {"z2p5": lam_lo, "candidate": lam, "z3p2": lam_hi},
            "affine_fraction_between_planes": (lam - lam_lo) / (lam_hi - lam_lo),
            "note": "Plausibility listing only. Multiz uses the same refined table and C23 covariance.",
            "rows": rows}


def member_path(name: str) -> str:
    member = TARGETS[name].member
    if member is None:
        raise Refusal(f"{name} has no assembled member")
    return member


def write_member(path: Path, rows: dict, meta: dict) -> dict:
    """Write one member npz (float64 x/y, metadata_json) and check it reloads exactly."""
    arrays = {"x": np.asarray(rows["x"], float), "y": np.asarray(rows["y"], float), "t_final": rows["t_final"],
              "a": rows["a"], "b": rows["b"], "order": rows["order"], "value": rows["value"],
              "metadata_json": np.array(canonical(meta).decode())}
    with path.open("xb") as handle:
        np.savez(handle, **arrays)
    with zipfile.ZipFile(path) as archive:
        members = {name: sha_bytes(archive.read(name)) for name in archive.namelist()}
    with np.load(path, allow_pickle=False) as check:
        for key, value in arrays.items():
            if not np.array_equal(check[key], value) or check[key].dtype != value.dtype:
                raise Refusal(f"Assembled member {path}:{key} does not reload exactly")
    return {"sha256": digest(path), "member_sha256": members}


def assembly_members(plan: dict, prepared: dict, loaded: dict) -> dict:
    """Member path -> (rows, metadata, manifest extras) for every production target; writes nothing."""
    by_name = {t["name"]: t for t in plan["targets"]}
    members = {}
    for name in PRODUCTION:
        target, (rows, completion) = by_name[name], loaded[name]
        meta = {"scope": SCOPE, "target_id": target["id"], "source_redshift": plan["source"]["z"],
                "lambda_mpc": plan["source"]["lambda_mpc"], "ell_max": target["cutoff"], "model": target["model"],
                "leg_order": target["leg"], "table_sha256": target["table_sha256"], "gl": target["gl"],
                "n_jobs": target["n_jobs"], "raw_result": completion["result"], "plan_sha256": prepared["plan_sha256"],
                "export": "raw numeric members copied; object x/y converted to float64 (n,3)"}
        members[member_path(name)] = (rows, meta, {"target_id": target["id"], "completion": completion["completed"]})
    return members


def publish_assembly(out: Path, manifest: dict, members: dict) -> dict:
    """Write members and manifest into a sibling staging directory, then rename it to ``out``.

    The rename uses RENAME_EXCL, so an existing ``out`` (even an empty directory)
    refuses. A failure leaves ``out`` absent; the staging directory keeps the
    partial files and never blocks a retry under the same --out.
    """
    staging = out.parent / f".{out.name}.{os.getpid()}.{secrets.token_hex(4)}.staging"
    staging.mkdir(exist_ok=False)
    try:
        outputs = {}
        for member, (rows, meta, extra) in members.items():
            path = staging / member
            path.parent.mkdir(exist_ok=True)
            outputs[member] = {"path": str(out / member), **write_member(path, rows, meta), **extra}
        manifest = {**manifest, "outputs": outputs}
        write_new_json(staging / "assembly_manifest.json", manifest)
        try:
            exclusive_rename(staging, out)
        except FileExistsError as exc:
            raise Refusal(f"--out appeared during assembly: {out}") from exc
    except BaseException:
        print(f"Assembly failed; partial files remain in {staging}; --out {out} was not created", file=sys.stderr)
        raise
    check_records([{"path": item["path"], "sha256": item["sha256"]} for item in outputs.values()])
    return manifest


def assemble_command(ns: argparse.Namespace) -> int:
    """Validate the z_s=3 production folds and the gate, then publish the ladder members and manifest."""
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    out = Path(os.path.abspath(ns.out)).resolve()
    if not inside(out, HERE) or out.exists() or out.is_symlink() or not out.parent.is_dir():
        raise Refusal(f"--out must be a new directory inside {HERE}: {out}")
    if not is_production_plan(plan):
        raise Refusal(f"Assembly takes only the z_s={PRODUCTION_Z:g} batch on this driver's source map")
    by_name = {t["name"]: t for t in plan["targets"]}
    if sorted(by_name) != sorted(PRODUCTION) or any(t["gl"] != GL_ACCEPTED for t in plan["targets"]):
        raise Refusal(f"Assembly needs exactly the GL{GL_ACCEPTED} production targets {PRODUCTION}")
    require_current(plan, prepared)
    gate = gate_evidence(gate_report_path(ns), plan)
    loaded = {name: completed_rows(by_name[name], prepared["plan_sha256"]) for name in PRODUCTION}
    members = assembly_members(plan, prepared, loaded)
    manifest = {"schema": 1, "kind": "zs3p0 cutoff ladder candidate assembly", "scope": SCOPE,
                "status": "CANDIDATE_NO_SCIENTIFIC_ACCEPTANCE", "source": plan["source"],
                "driver": record(DRIVER), "prepared": record(work / "prepared.json"),
                "plan_sha256": prepared["plan_sha256"], "gate": gate,
                "highest_tree_origin": "tree_main_15360 fold (refined main table, 11 ladder directions)",
                "plausibility": plausibility(loaded["tree_main_15360"][0], plan["source"]["lambda_mpc"]),
                "caveats": ["T-001 acceptance covers finite amplitude and sampled sign at z_s=5 only.",
                            "The 16-shell radial grid is concentrated near z=5; radial accuracy at z_s=3 is not "
                            "established.", "C23 source-grid accuracy is not granted.",
                            "No GL24/GL36 comparison at z_s=3 unless the optional controls are run.",
                            "Quote every FK value with its l_max."],
                "active_manifest_modified": False}
    manifest = publish_assembly(out, manifest, members)
    print(json.dumps({"assembly_manifest": record(out / "assembly_manifest.json"),
                      "outputs": {k: v["sha256"] for k, v in manifest["outputs"].items()}}, indent=2))
    return 0


# -------------------------------------------------------------------- CLI
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("source-map", help="write the z_s=3.0 source endpoint record")
    p.add_argument("--z", type=float, required=True)
    p.add_argument("--out", type=Path, default=None, help="default: <driver dir>/source_map_<tag>.json")
    p = sub.add_parser("prepare", help="materialise configs, loaders and contracts; no fold")
    p.add_argument("--source-map", type=Path, required=True)
    p.add_argument("--z", type=float, required=True)
    p.add_argument("--gl", type=gl_arg, required=True)
    p.add_argument("--n-jobs", type=n_jobs_arg, required=True,
                   help="sweep workers: a positive integer, or 'accepted' (per-target T-001 record)")
    p.add_argument("--targets", nargs="+", choices=sorted(TARGETS), required=True)
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--allow-stage-test", action="store_true",
                   help="TEST ONLY: allow a driver outside its planned r1_sft061 directory; such batches cannot run")
    gate_help = f"PASS gate report (default: <driver dir>/{GATE_WORK_NAME}/gate_report.json)"
    p = sub.add_parser("run", help="run one prepared fold under the shared lock")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--min-available-gib", type=float, default=None)
    p.add_argument("--gate-report", type=Path, default=None, help=gate_help + "; not used by a z_s=5 batch")
    p = sub.add_parser("gate", help="compare z_s=5 gate folds with the accepted products")
    p.add_argument("--work-dir", type=Path, required=True)
    p = sub.add_parser("assemble", help="write the candidate ladder members and manifest")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--gate-report", type=Path, default=None, help=gate_help)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    ns = parse_args(argv)
    try:
        if ns.command == "source-map":
            return source_map_command(ns.z, ns.out or HERE / f"source_map_{tag_for(ns.z)}.json")
        return {"prepare": prepare_command, "run": run_command, "gate": gate_command,
                "assemble": assemble_command}[ns.command](ns)
    except Refusal as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    except (OSError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"REFUSED ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
