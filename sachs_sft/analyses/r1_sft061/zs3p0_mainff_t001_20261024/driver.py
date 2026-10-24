"""Guarded z_s = 3.0 Order-0, FF and main-FK replay for the NHFP figure (version 2).

Candidate z_s=3 finite products, no scientific acceptance.

Version 2 of zs3p0_ladder_t001_20261023/driver.py. It keeps that driver's guards
and reuses its z_s = 3.0 source map, v1 products and the T-001 inputs. Pipelines:

  order0  SFT source-folded scalar convergence (source_folded_reference.py, PyCCL
          interpreter): parity at z=3.2 (recorded argv) and at z=5 (full channels),
          then the z=3.0 production, five-multipole fine probe and full-channel source
  ff      C23 FF, Gaussian driving fields, pair (0,0), GL24: canary at z_s=5 on main
          directions 0 and 6, production at z_s=3.0 on all 40 directions
  fk      tree-level main FK, refined table, l_max=15360, GL36, main_guard: regression
          at z=3.2 on main directions 0,5,10,15,25,35, production at z_s=3.0 (40 rows)

Subcommands
  plan         print one stage's plan; reads and hashes, writes nothing
  prepare      materialise one stage in a new directory; no fold and no Order-0 job
  run          one fold or one Order-0 job under the shared lock; a production
               target needs its pipeline's PASS gate report
  gate         compare a gate stage with the accepted products (writes gate_report.json)
  consistency  compare the z_s=3.0 main FK with the v1 tree_main_15360 fold
  export       write the NHFP package-style copies and export_manifest.json

Reviewed modules (main_guard.py, continuation_guard.py, batch_driver.py,
run_fk_variant.py, replay_suite.py, aligned_response.py) are imported by path and
never edited. Every generated file is new, in a fresh directory. Every record binds
the sha256 of this file. Run with the sft-wick interpreter, -B,
PYTHONDONTWRITEBYTECODE=1, PYTHONPATH=<frozen framework export>:$CANOES_ROOT/src and
CANOES_ROOT set.
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
PLANNED_HOME = R1 / "zs3p0_mainff_t001_20261024"
V1 = R1 / "zs3p0_ladder_t001_20261023"
V1_DRIVER = V1 / "driver.py"
V1_SOURCE_MAP = V1 / "source_map_zs3p0.json"
V1_TREE = V1 / "folds_zs3p0_gl36/tree_main_15360_zs3p0_gl36"
V1_TREE_RESULT = V1_TREE / "xi_tree_main_15360_zs3p0_gl36.npz"
V1_LADDER = V1 / "cutoff_ladder_zs3p0"
LADDER_MEMBERS = ("tree_auto/table_tree_cut7680_r4_xi.npz", "tree_auto/table_tree_cut15360_r4_xi.npz",
                  "bihalofit_auto/table_bihalofit_cut7680_r4_xi.npz",
                  "bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz")
LADDER_TOP = "tree_auto/table_tree_cut15360_r4_xi.npz"
CAND = R1 / "t001_candidate_inputs"
V6 = CAND / "main_fk_gl36_optimized_v6"
FW = V6 / "framework_export/src"
MAIN_EXEC = V6 / "execution/main_fk_gl36"
REVIEW = CAND / "point6_preflight/review_batch"
CONT = CAND / "point6_preflight/continuation_batch"
MULTIZ_V10 = CAND / "final_multiz_mixed_v10_20261005/assembly_manifest.json"
MAIN_V11 = CAND / "final_main_v11_20261005/assembly_manifest.json"
ACCEPTED = R1 / "t001_accepted"
KV = SACHS / "callables/kappa3_vertex"
PIECES = KV / "rebuild/products/pieces_t001"
ACCEPTANCE = PIECES / "acceptance"
FF_MAIN_RELEASE = ACCEPTANCE / "ff_maincap66_full_main_preparation_20261005/release.json"
FF_Z3P2_RELEASE = ACCEPTANCE / "ff_singlepair24_selected_releases_20261005/full_ff_z3p2_gl24/release.json"
PCHIP = KV / "equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable_pchip.py"
REFINED = PIECES / "table_permclosed_np512_t001_refined.npz"
CANOES = Path("/Users/zzhang/projects/angular_statistics/canoes")
C23_TABLE = CANOES / "docs/handoff/2026_10_03_T-001_delivery/corr_op_table_zs5_23pt_stf_fid_t001.npz"
C_WRAPPER = CAND / "corr_op_C_callable.py"
PYTHON = Path("/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python")
PYCCL = Path("/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python")
LOCK = R1 / ".replay_suite.lock"
REVIEWED_SOURCE_MAP = R1 / "source_map.json"
R1A = SACHS / "analyses/analysis1/r1_aligned"
O0_SCRIPT = R1A / "source_folded_reference.py"
O0_CCL = R1A / "ccl_reference.py"
O0_ANGLES = R1A / "outputs/order0_table_integral/order0_observables.npz"
O0_PLAN = R1A / "source_folded_multiz_plan.json"
O0_Z3P2 = R1A / "source_folded_multiz_z3p2.npz"
O0_Z3P2_JSON = R1A / "source_folded_multiz_z3p2.json"
O0_Z3P2_PROBE_JSON = R1A / "source_folded_multiz_z3p2_probe.json"
O0_Z5 = R1A / "source_folded_full_boundary_n4096_k2048.npz"
O0_Z5_JSON = R1A / "source_folded_full_boundary_n4096_k2048.json"
O0_FIG = R1 / "figure_products/order0.npz"
O0_FIG_JSON = R1 / "figure_products/order0.json"
O0_PLOT_CONFIG = SACHS / "sftwick_outputs/2PCF/C_corr_op_O0/config_L2.yaml"
PREPARE_PLOT = R1 / "prepare_plot_inputs.py"
PCAMB = CANOES / "examples/data/PCAMB_pyccl_stf_fid_z0.txt"
#: NHFP generator files whose reading conventions the export reproduces. Recorded, not
#: pinned: the generator may change when the figure is switched to z_s=3.
NHFP_REFERENCES = (Path("/Users/zzhang/MyVault/Career/job2026/NHFP/figures/src/nhfp_fig_evidence.py"),
                   Path("/Users/zzhang/MyVault/Career/job2026/NHFP/figures/src/ccl_born_limber.py"))
EXPECTED_PYTHONPATH = f"{FW}{os.pathsep}{CANOES / 'src'}"
FRAMEWORK_COMMIT = "5198fb9845334982d103059ac7b84f7afb2228ca"
CANOES_COMMIT = "53850f3423ee25b9c1813cff8c2a2b282c017c70"
EMPTY_DIFF_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
MAIN_GEOMETRY = "6102de8ee3fef96fcafda2773046e71f26394b55c1b15d1f56513e9e567996d3"
#: The three source planes this driver serves (v1 map for z_s=3.0, reviewed map otherwise).
Z_PRODUCTION, Z_REGRESSION, Z_CANARY = 3.0, 3.2, 5.0
EXPECTED_LAMBDA = {Z_PRODUCTION: 2255.3607197066854, Z_REGRESSION: 2266.590169124188,
                   Z_CANARY: 2317.9695958203943}
TAGS = {Z_PRODUCTION: "zs3p0", Z_REGRESSION: "zs3p2", Z_CANARY: "zs5"}
MAX_AFFINE_MPC = 2353.929874188644
ROUNDTRIP_ATOL = 5e-13
N_MAIN = 40
LADDER_INDICES = (0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
REGRESSION_INDICES = (0, 5, 10, 15, 25, 35)
CANARY_INDICES = (0, 6)
REL_TOL = 1e-12
ORDER0_PEAK_TOL = 1e-5
FINE_PROBE_FLAG = 5e-3
RESOLUTION_FLAG = 1e-3
NODE_ATOL_ARCMIN = 1e-9
ORDERS = "0,2"
SCOPE = "candidate z_s=3 finite products, no scientific acceptance"
MARKER = 'TABLE_PATH = _HERE / "table_permclosed.npz"'
GUARD_MARKERS = ("\nimport importlib.util as _guard_import\n", "\nimport importlib.util as _gi\n")
BOOLEAN_FLAGS = frozenset({"--cached-direct", "--observer-boundary"})
#: renamex_np flag (macOS <stdio.h>): fail with EEXIST instead of replacing the target.
RENAME_EXCL = 0x00000004
#: Peak memory per (pipeline, sweep workers). FK: T-001 GL36 sampled group RSS (main v6
#: resource_completion.json for 4, final_ladder_v8 admissions for 8 and 12). FF: the reviewed
#: C23 GL24 four-worker cap (16 GiB), a bound rather than a measurement. Order-0: the
#: recorded z=5 full-channel run (one process).
MEASURED_PEAK_BYTES = {("fk", 4): 27544076288, ("fk", 8): 50694569984, ("fk", 12): 67749150720,
                       ("ff", 4): 17179869184, ("order0", 1): 1128103936}
MEMORY_FLOOR_BYTES = 12 * 2**30
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
#: The recorded multiz Order-0 environment (source_folded_multiz_plan.json "environment").
ORDER0_ENV = {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "CANOES_SUPPRESS_METAL_WARNING": "1"}
ORDER0_TIMEOUT_S = 7200
#: Process names that must not run during a fold or an Order-0 job, beyond the reviewed list.
EXTRA_ACTIVE = ("source_folded_reference.py", "ccl_reference.py", "replay_suite.py", "companion_suite.py",
                "permaware_ladder/pieces.py", "permaware_ladder/suite.py", "zs3p0_ladder_t001_20261023/driver.py")

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
    C_WRAPPER: "1d0c80839d7a3f3866b5fff36d07b9e08d945d7005cc5b5ec5183b2c21906dd2",
    C23_TABLE: "d487c51b9dcf1d87520d82eac7d37aee17b872b3c83356f8007ff82d472562ee",
    SACHS / "scripts/inputs/F_tensor.npy": "63b682be4315a5783d1554e7ea48364229e8bcf17181c6482d14269101e1d0db",
    SACHS / "scripts/kappa2_callable.py": "49372414bab52dce803f786c301c07219b638b39f5038f60e0bb9858ce18cb29",
    SACHS / "scripts/inputs/kappa2_grid_lcut0_low_empty.npz":
        "1e03b492d46ed93db3a4ad6d953323325a20322470371412f1ac3587b1139e63",
    SACHS / "scripts/run_FF_single.py": "5ca112bbcad3e38ff8ee6f6c6544661dc962eeeff7aea3423a0f5dd9abc7e443",
    SACHS / "scripts/spin_rotation.py": "0322515d53b0dcc8832e585fb8bed3b242d465873e79f71440adf7e0f71abc8c",
    KV / "rebuild/run_fk_variant.py": "0ea02c4ce96cd77f1c146a5cc38214a059530b8833f8a61b005b490ec84407b2",
    CANOES / "src/canoes/cosmo/lightcone.py": "69f57f5e95545c094a4ad95e70a603635f5a886178df5e82703191a849fb3d02",
    CANOES / "src/canoes/cosmo/growth.py": "c931cca5f6d66b4b3082a6a0a11924604d12583219bbe31ca8df804cd371e08c",
    PYTHON.resolve(): "17f40a0285f0eddfc7eb226d625390ecd8fbd3ed907ae37bb98e81881ee95784",
    PYCCL.resolve(): "6c93640f2dd0fc2b45ff0fca2b233da318921099893b897f37a029e9c29d1c25",
    V6 / "source_snapshot.json": "19d794073859d9addcd153039f9cb3b6198a5871fe4a179935b1260143105d16",
    MAIN_EXEC / "prepared.json": "e319cc4cc695bbdfde15259b1fc85225921195775d8205cf460aff7ca9004dd6",
    MAIN_EXEC / "config_main.yaml": "60e3268fe80170cc594f3d186b92f2454b5ad1bcf277b6918324ab44ad75ef20",
    MAIN_EXEC / "callable_main.py": "65e136196e80df18a9f1aee47d7d9742f6069b4a9723fc0872bd3b8c548f6cde",
    MAIN_EXEC / "main_contract.json": "f39dd012e356968328074981180d40824d23a97be8ba886cd823ba3206f07e86",
    MAIN_EXEC / "completed.json": "e88e2d2c37227dbfca82a86a4fd963720576277f9285ff5a0268d701ea73ed50",
    MAIN_EXEC / "xi_main.npz": "d4bbc17d26db1ca50c27c8316238edcca02c998e8882ab88116be371b1d10e3a",
    MULTIZ_V10: "0c1390bce3ea73683e00edbf50eb0efe2c22cd7f2dd842a36359111512e98ec2",
    MAIN_V11: "f34f19d850585c792adb2bf57031d81124fb70a2bd367bfad690682bbe35e5c3",
    FF_MAIN_RELEASE: "6e2b0a5587392da380385004e0740af929265b614b0080bc44955a3ada42e1ce",
    FF_Z3P2_RELEASE: "f124e90bd989fb0dd3360a72b54e62cd9e3102eeca4fc7e4826e7ca218ab938c",
    ACCEPTED / "products/main_fk.npz": "d2c6f5c434c2ab6970cea7e42091e3e7ab9f0424a9c63e7efb801124809137e3",
    ACCEPTED / "products/main_ff.npz": "340aff7ae5688b5e135af7769c1809f97ef03de6701f22f84a2ffb8ee3a44fae",
    ACCEPTED / "products/multiz.npz": "4bc94805b5c57939a55454a1c706ec8f718d848c8872a5ecffaabed8b22f71cb",
    ACCEPTED / "cutoff_ladder" / LADDER_MEMBERS[0]:
        "cddd936c711f0de6cef7535b884885d54df6966ce78c8473d151b7ae38c77d1c",
    ACCEPTED / "cutoff_ladder" / LADDER_MEMBERS[1]:
        "e89b13a65e0608d6f65df5588a1c503f97891dd0920fd6458a572b1bd911d5eb",
    ACCEPTED / "cutoff_ladder" / LADDER_MEMBERS[2]:
        "c0033c497ca95843f8563f627fb63b36f7a2f05d009555dfa258875313474fce",
    ACCEPTED / "cutoff_ladder" / LADDER_MEMBERS[3]:
        "eb504b9fa7bd66fdda1b64cff3d9c7fd1aeef4b99372c7fcb37762a45794f792",
    V1_DRIVER: "dd6ae0af0b598cb374d6a5e4e131ada1cb8a7da234208e17ed020ac81f2c642b",
    V1_SOURCE_MAP: "bca33d13f261cfe2772d83b40d5a6bebafe2a7f9056a2229ba3ac85b172e1008",
    V1 / "folds_zs3p0_gl36/prepared.json": "4ccd967dcbc349fbe7eb044506c923514cc87de4b8f73f3ab90c06a205769954",
    V1_TREE / "started.json": "0d56f3a5285d0da966ef15ffa7a32b2e0e3e14f758f43a5096c3bf65aedf6607",
    V1_TREE / "completed.json": "048c767558080d96fe235a33a26d0ee093b054586a52a286216bb9cc481c6d3f",
    V1_TREE_RESULT: "5cc27a56bdd87258a4742ae857f22b3cf790341380bfabb77ddab1d28729a681",
    V1 / "gates_zs5_gl36/gate_report.json": "d7c711739ab348fa55ea7c00dab2566f90736d0f284e62e81a87fbd4de9094ed",
    V1_LADDER / "assembly_manifest.json": "aefd3085cc0a11327e03fa564938a9719a6ee97cc1d01177c595ea1b52a97dde",
    V1_LADDER / LADDER_MEMBERS[0]: "7db412fb71315b05bceb4888aa9b6ebc95b745135632c072c76c348cb2c77253",
    V1_LADDER / LADDER_MEMBERS[1]: "799b7c9ec2211443eb120c267b6646b8594fa9081380dd4c8898abfbddf34775",
    V1_LADDER / LADDER_MEMBERS[2]: "2d0cd15cc1006fc070bcb4df344a993870de623d5cf50e10884e265c8dc07c16",
    V1_LADDER / LADDER_MEMBERS[3]: "078ba1c069225911f9079e2975e8aa2698f047800085adbc0a330740246bbaa5",
    O0_SCRIPT: "47b855b1b438ca323f85e3bcf9fe9aee466d38dc9c5830423624eb7b70913424",
    O0_CCL: "520f5bc26ba0465b7588d63bb881ead33ca35004f24e23f952390428538ed57a",
    O0_ANGLES: "5c4f13a7042e83275d62ad1e8c1559ec731bd09d5a37dce7eaeb865c46701b5e",
    O0_PLAN: "a4202e8ad4c7f987b1ed03b7a3d6f6f88dc5e3f75d8918091fb6880033334912",
    O0_Z3P2: "9b13642eae19b1a091540ee77e59b7d4e53d547f5195124c9b7f8791aa654f6f",
    O0_Z3P2_JSON: "17b6669cc181a67fd7e36ab4b3d0192016e0921cc2a01d4aded569c4faabd2b6",
    O0_Z3P2_PROBE_JSON: "f39ef68c214e93beb97844ed60a24ca6c70c43c84342da262095a50a5e602d11",
    O0_Z5: "6dbc612cddd2e9e29591e6b3465ac71bddfe2dbf8cc2b3a45d4b4756b8945a14",
    O0_Z5_JSON: "f723d917316f36d995747c633d2ffa3681f1d779baa8aeaf5dd184f136cbe77c",
    O0_FIG: "975bd9564e1d503d096fd9d57d291d0987f2e35df51f861a39ffb04e35115c87",
    O0_FIG_JSON: "38e93fc3268e51dfb1c61b792857eec52cb9f03410a3870d01a93b8e676629b1",
    O0_PLOT_CONFIG: "94f82602970b3a20ae5df41c8e03bacb2fba40812d17b627e6c7a65651c8958a",
    PREPARE_PLOT: "377619f153512036a74ffc2bffc0d46db618f406ab5f54d68624c2fb564788c9",
    PCAMB: "2e5600748706ff4b67c99df7f9331b2a6ee1fa8623324e25dde836681f8f9776",
}
#: Inputs pinned in every batch (fold and Order-0 alike), in addition to this driver,
#: the framework export and the canoes sources.
SHARED_INPUTS = [R1 / "aligned_response.py", R1 / "replay_suite.py", PKG / "reproduce/archived_inputs.py",
                 REVIEW / "main_guard.py", REVIEW / "batch_driver.py", CONT / "continuation_guard.py", PCHIP,
                 C_WRAPPER, C23_TABLE, SACHS / "scripts/inputs/F_tensor.npy", SACHS / "scripts/kappa2_callable.py",
                 SACHS / "scripts/inputs/kappa2_grid_lcut0_low_empty.npz", SACHS / "scripts/run_FF_single.py",
                 SACHS / "scripts/spin_rotation.py", KV / "rebuild/run_fk_variant.py", PYTHON.resolve(),
                 O0_SCRIPT, O0_CCL, O0_ANGLES, PCAMB, PYCCL.resolve()]
REFERENCE_FILES = sorted(set(PINNED) - set(SHARED_INPUTS), key=str)


# ------------------------------------------------------------------- stages and targets
@dataclasses.dataclass(frozen=True)
class Stage:
    """One preparable batch: its pipeline, gate or production role, targets and default directory."""

    pipeline: str
    role: str
    targets: tuple[str, ...]
    default_dir: str


@dataclasses.dataclass(frozen=True)
class FoldSpec:
    """One sft-wick fold through run_fk_variant.run -> run_FF_single.py."""

    pipeline: str
    z: float
    source_map: Path
    base: str
    y_indices: tuple[int, ...] | None
    gl: int
    accepted: tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class Order0Spec:
    """One source_folded_reference.py job: the recorded z=3.2 argv with explicit flag changes.

    ``source_map`` None keeps the recorded argv, whose default is the reviewed map.
    ``reference`` is the accepted product JSON whose recorded settings the job reproduces
    apart from out, z_source and source_map.
    """

    z: float
    source_map: Path | None
    reference: Path
    changes: tuple[tuple[str, str], ...]
    arrays: tuple[str, ...]


STAGES = {
    "order0_parity": Stage("order0", "gate", ("o0_parity_zs3p2", "o0_parity_full_zs5"), "order0_parity"),
    "order0_production": Stage("order0", "production", ("o0_prod_zs3p0", "o0_probe5_zs3p0", "o0_full_zs3p0"),
                               "order0_production_zs3p0"),
    "ff_canary": Stage("ff", "gate", ("ff_canary_zs5",), "ff_canary_zs5_gl24"),
    "ff_production": Stage("ff", "production", ("ff_zs3p0",), "ff_production_zs3p0_gl24"),
    "fk_regression": Stage("fk", "gate", ("fk_regress_zs3p2",), "fk_regression_zs3p2_gl36"),
    "fk_production": Stage("fk", "production", ("fk_main_zs3p0",), "fk_production_zs3p0_gl36"),
}
GATE_STAGE = {"order0": "order0_parity", "ff": "ff_canary", "fk": "fk_regression"}
PRODUCTION_STAGE = {"order0": "order0_production", "ff": "ff_production", "fk": "fk_production"}
FOLDS = {
    "ff_canary_zs5": FoldSpec("ff", Z_CANARY, REVIEWED_SOURCE_MAP, "ff_main_release", CANARY_INDICES, 24,
                              ("ff_main_release", "ff_main_v11")),
    "ff_zs3p0": FoldSpec("ff", Z_PRODUCTION, V1_SOURCE_MAP, "ff_z3p2_release", None, 24,
                         ("ff_z3p2_release", "ff_main_release")),
    "fk_regress_zs3p2": FoldSpec("fk", Z_REGRESSION, REVIEWED_SOURCE_MAP, "main_fk_gl36", REGRESSION_INDICES, 36,
                                 ("main_fk_gl36", "multiz_production_gl36")),
    "fk_main_zs3p0": FoldSpec("fk", Z_PRODUCTION, V1_SOURCE_MAP, "main_fk_gl36", None, 36, ("main_fk_gl36",)),
}
SPARSE_ARRAYS = ("ell", "cl_pp", "sampled_ell", "sampled_cl_pp", "theta_arcmin", "kk", "ell_taper")
FULL_ARRAYS = (*SPARSE_ARRAYS, "cl_pw", "cl_ww", "xip", "xim", "kgt")
PROBE_ARRAYS = ("ell", "cl_pp")
PROBE_CHANGES = (("--ells", "2,10,100,500,800"), ("--n-chi", "4096"), ("--low-ell-n-k", "2048"), ("--ell-block", "4"))
FULL_CHANGES = (("--n-chi", "4096"), ("--low-ell-n-k", "2048"), ("--ell-block", "4"),
                ("--channels", "cl_pp,cl_pw,cl_ww"))
ORDER0 = {
    "o0_parity_zs3p2": Order0Spec(Z_REGRESSION, None, O0_Z3P2_JSON, (), SPARSE_ARRAYS),
    "o0_parity_full_zs5": Order0Spec(Z_CANARY, REVIEWED_SOURCE_MAP, O0_Z5_JSON, FULL_CHANGES, FULL_ARRAYS),
    "o0_prod_zs3p0": Order0Spec(Z_PRODUCTION, V1_SOURCE_MAP, O0_Z3P2_JSON, (), SPARSE_ARRAYS),
    "o0_probe5_zs3p0": Order0Spec(Z_PRODUCTION, V1_SOURCE_MAP, O0_Z3P2_PROBE_JSON, PROBE_CHANGES, PROBE_ARRAYS),
    "o0_full_zs3p0": Order0Spec(Z_PRODUCTION, V1_SOURCE_MAP, O0_Z5_JSON, FULL_CHANGES, FULL_ARRAYS),
}
#: Order-0 parity: reference product, arrays that must be bitwise equal, arrays held to
#: ORDER0_PEAK_TOL of their peak, arrays only reported.
ORDER0_GATES = {
    "o0_parity_zs3p2": (O0_Z3P2, ("theta_arcmin", "ell_taper", "ell", "sampled_ell"), ("kk",),
                        ("cl_pp", "sampled_cl_pp")),
    "o0_parity_full_zs5": (O0_Z5, ("theta_arcmin", "ell_taper", "ell", "sampled_ell"), ("kk", "xip", "xim", "kgt"),
                           ("cl_pp", "cl_pw", "cl_ww")),
}
GATE_MATCH_KEYS = ("gl", "n_jobs", "guard", "c_module", "table_sha256")
CONFIG_CATEGORIES = {
    "expand.cache_path": "path", "propagators.cache_path": "path", "output.0.path": "path",
    "system.nonlocal_vertices.0.coupling_module": "path", "propagators.c_closed_form_module": "c_module",
    "sweep.t_final_grid": "endpoint", "propagators.t_max": "endpoint",
    "sweep.positions_grid.y": "direction_subset", "sweep.component_pairs": "pair_subset",
    "sweep.n_jobs": "workers", "sweep.n_gauss": "quadrature",
}
CONTRACT_CATEGORIES = {"query_interval": "endpoint", "directions": "direction_subset", "scope": "guard_scope",
                       "model": "guard_scope", "cutoff": "guard_scope", "leg_order": "guard_scope"}
SETTINGS_CATEGORIES = {"out": "path", "z_source": "endpoint", "source_map": "endpoint"}
#: Allowed difference categories of each materialised fold against each accepted record.
ALLOWED = {
    ("fk_regress_zs3p2", "main_fk_gl36"): {"path", "endpoint", "direction_subset", "pair_subset"},
    ("fk_regress_zs3p2", "multiz_production_gl36"): {"path", "endpoint", "direction_subset", "workers",
                                                     "guard_scope"},
    ("fk_main_zs3p0", "main_fk_gl36"): {"path", "endpoint", "pair_subset"},
    ("ff_canary_zs5", "ff_main_release"): {"path", "c_module", "direction_subset", "pair_subset"},
    ("ff_canary_zs5", "ff_main_v11"): {"path", "c_module", "direction_subset", "pair_subset"},
    ("ff_zs3p0", "ff_z3p2_release"): {"path", "c_module", "endpoint"},
    ("ff_zs3p0", "ff_main_release"): {"path", "c_module", "endpoint", "pair_subset"},
}
#: Fields in which the two reviewed C23 FF release configs may differ.
FF_RELEASE_FREE_FIELDS = {"expand.cache_path", "propagators.cache_path", "output.0.path",
                          "propagators.c_closed_form_module", "propagators.t_max", "sweep.t_final_grid",
                          "sweep.component_pairs"}


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
    """Create ``path`` holding ``data``: complete or absent, never replacing anything."""
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


def write_new_json(path: Path, value: dict) -> None:
    """Indented JSON with a trailing newline, published atomically."""
    publish_new(path, (json.dumps(value, indent=2) + "\n").encode("utf-8"))


def inside(path: Path, parent: Path) -> bool:
    """``path`` lies strictly below ``parent`` after resolving symlinks."""
    return path.resolve().is_relative_to(parent.resolve()) and path.resolve() != parent.resolve()


def new_dir_arg(path: Path | None, default: str) -> Path:
    """Absolute, resolved directory argument with a default under this driver's directory."""
    chosen = Path(path) if path is not None else HERE / default
    return Path(os.path.abspath(chosen)).resolve()


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
    return load_module("zs3v2_review_batch_driver", REVIEW / "batch_driver.py", REVIEW)


def continuation_guard() -> ModuleType:
    return load_module("zs3v2_continuation_guard", CONT / "continuation_guard.py")


def replay() -> ModuleType:
    return load_module("zs3v2_replay_suite", R1 / "replay_suite.py")


def variant() -> ModuleType:
    return load_module("zs3v2_run_fk_variant", KV / "rebuild/run_fk_variant.py")


def aligned() -> ModuleType:
    if Path(os.environ.get("CANOES_ROOT", "")) != CANOES:
        raise Refusal("CANOES_ROOT must equal the recorded canoes root before importing aligned_response")
    return load_module("zs3v2_aligned_response", R1 / "aligned_response.py")


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


def order0_environment(root: Path | None) -> dict:
    """Caller environment without SFT_WICK_* or PYTHON* variables, plus the recorded Order-0 variables."""
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SFT_WICK_", "PYTHON"))}
    env.update(ORDER0_ENV, PYTHONDONTWRITEBYTECODE="1")
    if root is not None:
        env["MPLCONFIGDIR"] = str(root / "mpl_cache")
    return env


def order0_runtime(cwd: Path | None = None) -> dict:
    """Identity of the PyCCL interpreter that runs source_folded_reference.py (metadata only)."""
    code = ("import sys, json, importlib.metadata as m\n"
            "print(json.dumps({'python_version': sys.version, 'executable': sys.executable,"
            " 'dependencies': {n: m.version(n) for n in ('numpy', 'scipy', 'pyccl')}}))\n")
    done = subprocess.run([str(PYCCL), "-B", "-c", code], env=order0_environment(None), cwd=cwd or HERE,
                          capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise Refusal(f"PyCCL interpreter query failed: {done.stderr.strip()[-2000:]}")
    identity = json.loads(done.stdout)
    if Path(identity["executable"]).resolve() != PYCCL.resolve():
        raise Refusal(f"PyCCL interpreter resolves to {identity['executable']}")
    return {"python": str(PYCCL), **identity, "binary": record(PYCCL.resolve())}


def fold_environment(root: Path) -> dict:
    """Caller environment without SFT_WICK_* or PYTHON* variables, plus the pinned fold variables."""
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SFT_WICK_", "PYTHON"))}
    env.update(PYTHONPATH=EXPECTED_PYTHONPATH, CANOES_ROOT=str(CANOES), PYTHONDONTWRITEBYTECODE="1",
               PYTHONNOUSERSITE="1", MPLCONFIGDIR=str(root / "mpl_cache"))
    env.update({name: "1" for name in THREAD_VARS})
    return env


def subprocess_binding(env: dict, cwd: Path) -> None:
    """Import resolution inside the fold environment, with run_FF_single's path order (no execution)."""
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


# ------------------------------------------------------------ source endpoints
def read_source_map_checked(path: Path) -> dict:
    """replay_suite.read_source_map (the reviewed rules); a failure refuses."""
    try:
        return replay().read_source_map(Path(path))
    except (ValueError, KeyError, TypeError, OSError) as exc:
        raise Refusal(f"Source map {path} fails replay_suite.read_source_map: {exc}") from exc


def derive_endpoint(z: float) -> dict:
    """Recompute lambda_s(z) with aligned_response and reproduce the reviewed six endpoints bitwise."""
    ar = aligned()
    reviewed = read_json(REVIEWED_SOURCE_MAP)
    six = np.interp(ar.SOURCE_REDSHIFTS, ar._Z, ar._LAMBDA)
    exact = {"z": ar.SOURCE_REDSHIFTS.tolist() == reviewed["z"],
             "lambda_mpc": six.tolist() == reviewed["lambda_mpc"],
             "chi_mpc": np.asarray(ar.chi_of_lambda(six)).tolist() == reviewed["chi_mpc"],
             "D_mpc": np.asarray(ar.D_at(six)).tolist() == reviewed["D_mpc"]}
    if not all(exact.values()) or float(ar.MAX_AFFINE_MPC) != MAX_AFFINE_MPC:
        raise Refusal(f"aligned_response no longer reproduces the reviewed source map: {exact}")
    z = float(z)
    lam = float(np.interp(z, ar._Z, ar._LAMBDA))
    chi, dist, back = float(ar.chi_of_lambda(lam)), float(ar.D_at(lam)), float(ar.z_of_lambda(lam))
    if not (0.0 < lam < MAX_AFFINE_MPC) or abs(back - z) > ROUNDTRIP_ATOL:
        raise Refusal(f"Endpoint for z={z} fails the round trip: lambda={lam}, z(lambda)={back}")
    if lam != EXPECTED_LAMBDA[z]:
        raise Refusal(f"lambda_s({z}) = {lam!r} differs from the recorded {EXPECTED_LAMBDA[z]!r}")
    return {"z": z, "lambda_mpc": lam, "chi_mpc": chi, "D_mpc": dist}


def load_endpoint(source_map: Path, z: float) -> dict:
    """Serve exactly three endpoints: z=3.0 on the v1 map, z=3.2 and z=5 on the reviewed map."""
    source_map, z = Path(source_map).resolve(), float(z)
    if z not in EXPECTED_LAMBDA:
        raise Refusal(f"Only z in {sorted(EXPECTED_LAMBDA)} is served, not z={z}")
    mapping = read_source_map_checked(source_map)
    if source_map == REVIEWED_SOURCE_MAP:
        if z == Z_PRODUCTION:
            raise Refusal("The reviewed source map has no z_s=3.0 entry; z_s=3.0 is served by the v1 map")
    elif source_map == V1_SOURCE_MAP.resolve():
        producer = mapping.get("producer", {}).get("driver", {})
        if z != Z_PRODUCTION or producer != {"path": str(V1_DRIVER), "sha256": PINNED[V1_DRIVER]}:
            raise Refusal("The v1 source map serves only z_s=3.0 and must record the pinned v1 driver")
        if mapping.get("target_tag") != TAGS[Z_PRODUCTION] or digest(source_map) != PINNED[V1_SOURCE_MAP]:
            raise Refusal("The v1 source map differs from its pinned record")
    else:
        raise Refusal(f"Source map {source_map} is neither the reviewed map nor the v1 z_s=3.0 map")
    if Path(mapping["response_module"]) != R1 / "aligned_response.py" or mapping["response_attr"] != "response_R":
        raise Refusal("Source map response differs from aligned_response.response_R")
    hits = np.flatnonzero(np.isclose(mapping["z"], z, rtol=0, atol=1e-12))
    if len(hits) != 1:
        raise Refusal(f"Source map has no unique endpoint for z={z}")
    i = int(hits[0])
    entry = {k: float(mapping[k][i]) for k in ("z", "lambda_mpc", "chi_mpc", "D_mpc")}
    if derive_endpoint(z) != entry:
        raise Refusal(f"Source-map entry {entry} differs from the re-derived endpoint")
    return {**entry, "tag": TAGS[z], "source_map": record(source_map)}


def t_max_for(lam: float, dt: float) -> float:
    """Propagator horizon dt*ceil(lambda/dt): ceil(lambda) for FK (dt=1), the companion rule for FF (dt=2).

    The FK and FF configs use c_closed_form_only, so t_max enters only the cache key and
    the R probes, not the values.
    """
    return float(dt * math.ceil(lam / dt))


# -------------------------------------------------------------- accepted records
def main_grid() -> tuple[list, list]:
    """The accepted main reference direction and the 40 main directions, in main-grid order."""
    cfg = yaml.safe_load((MAIN_EXEC / "config_main.yaml").read_text(encoding="utf-8"))
    grid = cfg["sweep"]["positions_grid"]
    if grid["x"] != [[0.0, 0.0, 1.0]] or len(grid["y"]) != N_MAIN:
        raise Refusal("Accepted main grid is not one pole and 40 directions")
    return grid["x"][0], grid["y"]


def accepted_records() -> dict:
    """Accepted main FK, multiz FK and C23 FF records used for the field-by-field diffs."""
    main = {"id": "main_fk_gl36", "root": str(MAIN_EXEC),
            "config": yaml.safe_load((MAIN_EXEC / "config_main.yaml").read_text(encoding="utf-8")),
            "loader_text": (MAIN_EXEC / "callable_main.py").read_text(encoding="utf-8"),
            "contract": read_json(MAIN_EXEC / "main_contract.json"), "record": record(MAIN_EXEC / "config_main.yaml")}
    admissions = read_json(MULTIZ_V10)["admissions"]
    first, second = admissions[0]["target"], admissions[1]["target"]
    if (first["id"] != "main_fk_gl36" or first["config"] != main["config"]
            or first["loader_text"] != main["loader_text"]
            or second["id"] != "multiz_production_gl36"):
        raise Refusal("final_multiz_mixed_v10 admissions differ from the main execution files")
    multiz = {"id": second["id"], "root": second["root"], "config": second["config"],
              "loader_text": second["loader_text"], "contract": second["contract"], "record": record(MULTIZ_V10)}
    main_release, z3p2_release = read_json(FF_MAIN_RELEASE), read_json(FF_Z3P2_RELEASE)
    v11 = read_json(MAIN_V11)["admissions"][1]["target"]
    if (v11["id"] != "full_ff_main_gl24" or v11["config"] != main_release["job"]["config"] or v11["gl"] != 24
            or main_release["job"]["gl"] != 24 or z3p2_release["job"]["gl"] != 24):
        raise Refusal("The accepted main FF admission differs from the C23 main release record")
    ff = {"ff_main_release": {"id": main_release["job"]["id"], "config": main_release["job"]["config"],
                              "worker_count": main_release["worker_count"], "record": record(FF_MAIN_RELEASE)},
          "ff_main_v11": {"id": v11["id"], "config": v11["config"], "record": record(MAIN_V11)},
          "ff_z3p2_release": {"id": z3p2_release["job"]["id"], "config": z3p2_release["job"]["config"],
                              "worker_count": z3p2_release["worker_count"], "record": record(FF_Z3P2_RELEASE)}}
    return {"main_fk_gl36": main, "multiz_production_gl36": multiz, **ff}


def base_record(spec: FoldSpec) -> dict:
    return record(MAIN_EXEC / "config_main.yaml") if spec.base == "main_fk_gl36" else (
        record(FF_MAIN_RELEASE if spec.base == "ff_main_release" else FF_Z3P2_RELEASE))


def require_values(cfg: dict, expected: dict, label: str) -> None:
    for keys, value in expected.items():
        node: Any = cfg
        for key in keys:
            node = node[key]
        if node != value:
            raise Refusal(f"{label} {'.'.join(keys)}={node!r}, expected {value!r}")


def check_fk_base(cfg: dict) -> None:
    """The accepted main config: GL36, four workers, six pairs, one K vertex, z_s=5."""
    require_values(cfg, {("system", "linear", "R_time_module"): str(R1 / "aligned_response.py"),
                         ("system", "linear", "R_time_vectorized"): True,
                         ("propagators", "c_closed_form_module"): str(C_WRAPPER),
                         ("propagators", "c_closed_form_only"): True, ("propagators", "dt"): 1.0,
                         ("propagators", "t_max"): 2318.0, ("sweep", "t_final_grid"): [EXPECTED_LAMBDA[Z_CANARY]],
                         ("sweep", "n_gauss"): 36, ("sweep", "n_jobs"): 4, ("expand", "n_jobs"): 1,
                         ("propagators", "n_jobs"): 4, ("sweep", "vertex_types"): ["FK"],
                         ("expand", "orders"): [0, 2]}, "Accepted main FK config")
    if (len(cfg["sweep"]["component_pairs"]) != 6 or len(cfg["output"]) != 1
            or [v["name"] for v in cfg["system"]["nonlocal_vertices"]] != ["K"]):
        raise Refusal("Accepted main FK config must have six pairs, one output and one K vertex")


def check_ff_base(cfg: dict, name: str) -> None:
    """A reviewed C23 FF release config: F only, GL24, four workers, closed-form C."""
    lam, pairs, t_max = ((EXPECTED_LAMBDA[Z_CANARY], 6, 2318.0) if name == "ff_main_release"
                         else (EXPECTED_LAMBDA[Z_REGRESSION], 1, 2268.0))
    require_values(cfg, {("system", "linear", "R_time_module"): str(R1 / "aligned_response.py"),
                         ("system", "noise", "kappa2", "module"): str(SACHS / "scripts/kappa2_callable.py"),
                         ("system", "nonlocal_vertices"): [], ("sweep", "vertex_types"): ["F"],
                         ("propagators", "c_closed_form_only"): True, ("propagators", "dt"): 2.0,
                         ("propagators", "t_max"): t_max, ("sweep", "t_final_grid"): [lam],
                         ("sweep", "n_gauss"): 24, ("sweep", "n_jobs"): 4, ("expand", "n_jobs"): 1,
                         ("propagators", "n_jobs"): 4, ("expand", "orders"): [0, 2]}, f"C23 FF {name} config")
    if (cfg["system"]["vertices"] != [{"coupling_path": str(SACHS / "scripts/inputs/F_tensor.npy"), "name": "F"}]
            or len(cfg["sweep"]["component_pairs"]) != pairs or len(cfg["output"]) != 1
            or cfg["sweep"]["positions_grid"] != {"x": [main_grid()[0]], "y": main_grid()[1]}):
        raise Refusal(f"C23 FF {name} config has unexpected vertices, pairs, outputs or directions")


def check_ff_release_pair(accepted: dict) -> None:
    """The main and z=3.2 C23 releases differ only in endpoint, pairs, C module and owned paths."""
    left = dict(flatten(accepted["ff_main_release"]["config"]))
    right = dict(flatten(accepted["ff_z3p2_release"]["config"]))
    changed = {f for f in set(left) | set(right) if left.get(f, "<absent>") != right.get(f, "<absent>")}
    if not changed <= FF_RELEASE_FREE_FIELDS:
        raise Refusal(f"The C23 FF release configs differ beyond endpoint and paths: {sorted(changed)}")


# ------------------------------------------------------------------ pins
def pin_records(paths: set[Path]) -> list[dict]:
    unique = {str(p.resolve()): p for p in paths}
    return [{"path": key, "sha256": digest(unique[key])} for key in sorted(unique)]


def shared_pins() -> list[dict]:
    """Inputs common to every batch: reviewed modules, runners, Order-0 script, driver, framework, canoes."""
    paths = set(SHARED_INPUTS) | {DRIVER}
    paths |= {p for p in FW.rglob("*") if p.is_file()}
    paths |= set((CANOES / "src/canoes").rglob("*.py"))
    return pin_records(paths)


def stage_inputs(stage: str) -> set[Path]:
    """Stage-specific pinned inputs: source maps, base configs, tables, reference settings."""
    paths: set[Path] = set()
    for name in STAGES[stage].targets:
        if name in FOLDS:
            spec = FOLDS[name]
            paths |= {spec.source_map, Path(base_record(spec)["path"])}
            if spec.pipeline == "fk":
                paths |= {REFINED, MAIN_EXEC / "main_contract.json", MAIN_EXEC / "callable_main.py"}
        else:
            spec0 = ORDER0[name]
            paths |= {spec0.source_map or REVIEWED_SOURCE_MAP, spec0.reference, O0_PLAN}
    return paths


def batch_pins(stage: str, common: list[dict]) -> list[dict]:
    merged = {p["path"]: p for p in common + pin_records(stage_inputs(stage))}
    return [merged[key] for key in sorted(merged)]


def check_canoes_pins(pins: list[dict]) -> None:
    """canoes Python sources (tracked and untracked) equal the accepted main contract pins exactly."""
    prefix = str(CANOES / "src/canoes") + "/"
    reference = {p["path"]: p["sha256"] for p in read_json(MAIN_EXEC / "main_contract.json")["pins"]
                 if p["path"].startswith(prefix)}
    mine = {p["path"]: p["sha256"] for p in pins if p["path"].startswith(prefix)}
    if mine != reference:
        raise Refusal("canoes source files differ from the accepted T-001 main contract pins")


# ----------------------------------------------------- diffs against accepted
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


def loader_body(text: str) -> str:
    """Loader text before its single guard-install block."""
    hits = [text.find(marker) for marker in GUARD_MARKERS if marker in text]
    if len(hits) != 1:
        raise Refusal("Loader text must hold exactly one guard-install block")
    return text[:hits[0]]


def loader_diff(target: dict, ref: dict, same_guard: bool) -> dict:
    """Same guard: equal after contract binding. Other guard (multiz): equal before the guard block."""
    mine = target["loader_text"]
    lines = list(difflib.unified_diff(ref["loader_text"].splitlines(), mine.splitlines(), "accepted",
                                      "materialised", n=0, lineterm=""))
    if same_guard:
        def normal(text: str, path: str, sha: str) -> str:
            return text.replace(repr(path), "'<CONTRACT_PATH>'").replace(repr(sha), "'<CONTRACT_SHA256>'")
        equal = (normal(mine, str(Path(target["root"]) / "main_contract.json"), target["contract_sha256"])
                 == normal(ref["loader_text"], str(Path(ref["root"]) / "main_contract.json"),
                           sha_bytes(canonical(ref["contract"]))))
        return {"comparison": "equal after contract binding", "equal": equal, "raw_diff": lines}
    own = make_loader(Path(target["root"]) / target["files"]["contract"], target["contract_sha256"])
    return {"comparison": ("equal before the guard-install block (accepted multiz used continuation_guard); "
                           "the materialised guard block is the main_guard install"),
            "equal": loader_body(mine) == loader_body(ref["loader_text"]) and mine == own, "raw_diff": lines}


def fold_diff(target: dict, accepted: dict, accepted_id: str) -> dict:
    """Field-by-field diff with one accepted record; differences outside the allowance refuse."""
    ref = accepted[accepted_id]
    allowed = set(ALLOWED[(target["name"], accepted_id)])
    if target["n_jobs_source"] == "cli":
        allowed.add("workers")
    report: dict[str, Any] = {"accepted_id": accepted_id, "accepted_record": ref["record"],
                              "allowed_categories": sorted(allowed),
                              "config": leaf_diff(target["config"], ref["config"], CONFIG_CATEGORIES, allowed)}
    unexpected = [row["field"] for row in report["config"] if not row["allowed"]]
    if target["pipeline"] == "fk":
        strip = {k: v for k, v in target["contract"].items() if k != "pins"}
        report["contract"] = leaf_diff(strip, {k: v for k, v in ref["contract"].items() if k != "pins"},
                                       CONTRACT_CATEGORIES, allowed)
        report["pins"] = pin_diff(target["contract"]["pins"], ref["contract"]["pins"])
        report["loader"] = loader_diff(target, ref, accepted_id == "main_fk_gl36")
        unexpected += [row["field"] for row in report["contract"] if not row["allowed"]]
        if report["pins"]["changed_common"]:
            unexpected.append("pins changed: " + ", ".join(report["pins"]["changed_common"]))
        if not report["loader"]["equal"]:
            unexpected.append("loader text")
    else:
        report["contract"] = "no accepted contract: the C23 FF runner bound its inputs through its own closure"
    if unexpected:
        raise Refusal(f"{target['id']} differs from accepted {accepted_id} beyond the allowed fields: {unexpected}")
    return report


# ---------------------------------------------------------- fold targets
def materialise_fold_config(spec: FoldSpec, base: dict, root: Path, result: Path, lam: float,
                            n_jobs: int) -> dict:
    """Copy of the accepted config with owned paths, endpoint, pair (0,0), direction subset, GL and workers."""
    cfg = copy.deepcopy(base)
    cfg["expand"]["cache_path"] = str(root / "cache_expand")
    cfg["propagators"]["cache_path"] = str(root / "cache_propagators")
    cfg["output"][0]["path"] = str(result)
    cfg["propagators"]["t_max"] = t_max_for(lam, float(cfg["propagators"]["dt"]))
    cfg["sweep"]["t_final_grid"] = [lam]
    cfg["sweep"]["component_pairs"] = [[0, 0]]
    if spec.y_indices is not None:
        cfg["sweep"]["positions_grid"]["y"] = [base["sweep"]["positions_grid"]["y"][i] for i in spec.y_indices]
    cfg["sweep"]["n_gauss"] = spec.gl
    cfg["sweep"]["n_jobs"] = n_jobs
    if spec.pipeline == "fk":
        cfg["system"]["nonlocal_vertices"][0]["coupling_module"] = str(root / "callable_main.py")
    else:
        cfg["propagators"]["c_closed_form_module"] = str(C_WRAPPER)
    return cfg


def table_info() -> dict:
    """continuation_guard.validate_table on the refined main table, plus geometry and overlay checks."""
    try:
        info = continuation_guard().validate_table(REFINED, "multiz", "tree", 15360, "auto")
    except (ValueError, KeyError, OSError) as exc:
        raise Refusal(f"Table {REFINED} fails continuation_guard.validate_table: {exc}") from exc
    if info["geometry_sha256"] != MAIN_GEOMETRY or not info["refinement_overlay_present"]:
        raise Refusal(f"Unexpected geometry or missing refinement overlay in {REFINED}")
    return info


def make_fk_contract(info: dict, directions: list, t_max: float, pins: list[dict]) -> dict:
    """Same keys as the accepted main contract."""
    return {"pins": pins, "table": str(REFINED), "table_sha256": digest(REFINED), "shape": info["shape"],
            "shells": info["shells"], "geometry_sha256": info["geometry_sha256"], "directions": directions,
            "query_interval": [0.0, t_max]}


def make_ff_contract(directions: list, lam: float, pins: list[dict]) -> dict:
    return {"kind": "C23 FF fold inputs; no vertex guard (FF has no K vertex)", "pins": pins,
            "c_module": record(C_WRAPPER), "c_table": record(C23_TABLE), "directions": directions,
            "t_final": lam}


def make_loader(contract_path: Path, contract_sha: str) -> str:
    """Same construction as the accepted main loader: PCHIP source, TABLE_PATH, main_guard install."""
    source = PCHIP.read_text(encoding="utf-8")
    if source.count(MARKER) != 1:
        raise Refusal("Ambiguous TABLE_PATH marker in the PCHIP source")
    text = source.replace(MARKER, f"TABLE_PATH = Path({str(REFINED)!r})")
    text += "\nimport importlib.util as _guard_import\n"
    text += (f'_guard_spec = _guard_import.spec_from_file_location("candidate_main_guard", '
             f"{str(REVIEW / 'main_guard.py')!r})\n")
    text += "_guard_mod = _guard_import.module_from_spec(_guard_spec)\n_guard_spec.loader.exec_module(_guard_mod)\n"
    text += f"_guard_mod.install(globals(), {str(contract_path)!r}, {contract_sha!r})\n"
    return text


def fold_workers(spec: FoldSpec, base: dict, n_jobs_arg: int | str) -> tuple[int, str]:
    if n_jobs_arg == "accepted":
        return int(base["sweep"]["n_jobs"]), f"accepted:{spec.base}"
    return int(n_jobs_arg), "cli"


def build_fold_target(name: str, work: Path, args: dict, accepted: dict, pins: list[dict]) -> dict:
    """In-memory config, contract and loader of one fold, diffed against its accepted records."""
    spec = FOLDS[name]
    source = load_endpoint(spec.source_map, spec.z)
    ident = f"{name}_gl{spec.gl}"
    root = work / ident
    result = root / f"xi_{ident}.npz"
    base = accepted[spec.base]["config"]
    check_fk_base(base) if spec.pipeline == "fk" else check_ff_base(base, spec.base)
    n_jobs, n_source = fold_workers(spec, base, args["n_jobs"])
    cfg = materialise_fold_config(spec, base, root, result, source["lambda_mpc"], n_jobs)
    config_text = canonical(cfg).decode() + "\n"
    if yaml.safe_load(config_text) != cfg:
        raise Refusal("Canonical JSON config does not round-trip through YAML")
    directions = [cfg["sweep"]["positions_grid"]["x"][0], *cfg["sweep"]["positions_grid"]["y"]]
    if spec.pipeline == "fk":
        contract = make_fk_contract(table_info(), directions, cfg["propagators"]["t_max"], pins)
        files = {"config": "config_main.yaml", "contract": "main_contract.json", "loader": "callable_main.py"}
    else:
        contract = make_ff_contract(directions, source["lambda_mpc"], pins)
        files = {"config": "config_ff.yaml", "contract": "ff_contract.json", "loader": None}
    contract_sha = sha_bytes(canonical(contract))
    loader_text = make_loader(root / files["contract"], contract_sha) if spec.pipeline == "fk" else None
    target = {"id": ident, "name": name, "pipeline": spec.pipeline, "z": spec.z, "source": source,
              "root": str(root), "result": str(result), "files": files, "base_config": base_record(spec),
              "config": cfg, "config_sha256": sha_bytes(config_text.encode()), "contract": contract,
              "contract_sha256": contract_sha, "loader_text": loader_text,
              "loader_sha256": None if loader_text is None else sha_bytes(loader_text.encode()),
              "n_jobs": n_jobs, "n_jobs_source": n_source, "gl": spec.gl,
              "y_indices": None if spec.y_indices is None else list(spec.y_indices),
              "expected_rows": len(cfg["sweep"]["positions_grid"]["y"]),
              "guard": "main_guard" if spec.pipeline == "fk" else "none", "c_module": record(C_WRAPPER),
              "table_sha256": digest(REFINED) if spec.pipeline == "fk" else None}
    target["accepted_diff"] = {acc: fold_diff(target, accepted, acc) for acc in spec.accepted}
    return target


# ---------------------------------------------------------- Order-0 targets
def recorded_tokens() -> list[tuple[str, str | None]]:
    """The accepted z=3.2 production argv of source_folded_multiz_plan.json, without --out."""
    runs = [r for r in read_json(O0_PLAN)["production_runs"] if r["z"] == Z_REGRESSION]
    if len(runs) != 1:
        raise Refusal("source_folded_multiz_plan.json has no unique z=3.2 production run")
    argv = runs[0]["argv"]
    if argv[:2] != [str(PYCCL), str(O0_SCRIPT)] or argv[-2:] != ["--out", str(O0_Z3P2)]:
        raise Refusal("The recorded z=3.2 argv uses another interpreter, script or output")
    body, tokens, i = argv[2:-2], [], 0
    while i < len(body):
        flag = body[i]
        if not flag.startswith("--") or flag == "--out":
            raise Refusal(f"Unexpected token in the recorded argv: {flag!r}")
        if flag in BOOLEAN_FLAGS:
            tokens.append((flag, None))
            i += 1
        else:
            tokens.append((flag, body[i + 1]))
            i += 2
    return tokens


def order0_tokens(spec: Order0Spec) -> list[tuple[str, str | None]]:
    """Recorded tokens with the spec's value changes, z_s and an explicit --source-map."""
    changes = dict(spec.changes)
    if spec.z != Z_REGRESSION:
        changes["--z-source"] = f"{spec.z:g}"
    tokens = []
    for flag, value in recorded_tokens():
        tokens.append((flag, changes.pop(flag, value)))
        if flag == "--z-source" and spec.source_map is not None:
            tokens.append(("--source-map", str(spec.source_map)))
    if changes:
        raise Refusal(f"Order-0 changes name flags absent from the recorded argv: {sorted(changes)}")
    return tokens


def check_tokens(tokens: list[tuple[str, str | None]], expected: dict) -> None:
    """Every argv value parses to the expected recorded setting (type and value)."""
    for flag, value in tokens:
        key = flag[2:].replace("-", "_")
        if key not in expected:
            raise Refusal(f"Argv flag {flag} has no recorded setting")
        want = expected[key]
        if value is None:
            got: Any = True
        elif isinstance(want, bool):
            raise Refusal(f"Argv flag {flag} carries a value but its setting is boolean")
        else:
            got = type(want)(value)
        if got != want or type(got) is not type(want):
            raise Refusal(f"Argv {flag} {value!r} does not reproduce the recorded setting {key}={want!r}")


def build_order0_target(name: str, work: Path) -> dict:
    """In-memory argv, expected settings and job record of one Order-0 job."""
    spec = ORDER0[name]
    source_map = spec.source_map or REVIEWED_SOURCE_MAP
    source = load_endpoint(source_map, spec.z)
    root = work / name
    out = root / f"{name}.npz"
    tokens = order0_tokens(spec)
    reference = read_json(spec.reference)["settings"]
    expected = {**reference, "out": str(out), "z_source": spec.z, "source_map": str(source_map)}
    check_tokens(tokens, expected)
    diff = leaf_diff(expected, reference, SETTINGS_CATEGORIES, {"path", "endpoint"})
    unexpected = [row["field"] for row in diff if not row["allowed"]]
    if name == "o0_parity_zs3p2" and tokens != recorded_tokens():
        unexpected.append("argv differs from the recorded z=3.2 argv")
    if unexpected:
        raise Refusal(f"{name} differs from its accepted settings beyond out/z_source/source_map: {unexpected}")
    argv = [str(PYCCL), str(O0_SCRIPT)]
    for flag, value in tokens:
        argv += [flag] if value is None else [flag, value]
    argv += ["--out", str(out)]
    outputs = [str(out), str(out.with_suffix(".json")), str(out.with_suffix(".source.py"))]
    job = {"argv": argv, "env": ORDER0_ENV, "cwd": str(root), "expected_settings": expected, "outputs": outputs,
           "log": str(root / f"{name}.log"), "reference_settings": record(spec.reference), "z": spec.z,
           "lambda_mpc": source["lambda_mpc"]}
    return {"id": name, "name": name, "pipeline": "order0", "z": spec.z, "source": source, "root": str(root),
            "out": str(out), "outputs": outputs, "log": job["log"], "argv": argv, "expected_settings": expected,
            "arrays": list(spec.arrays), "job": job, "job_sha256": sha_bytes(canonical(job)),
            "accepted_diff": {"reference": record(spec.reference), "settings": diff}}


# ------------------------------------------------------------------ batch plan
def build_batch(args: dict) -> dict:
    """Deterministic in-memory batch plan; reads and hashes, never writes."""
    verify_pinned()
    stage = STAGES[args["stage"]]
    runtime, framework, canoes = runtime_state(), framework_state(), canoes_state()
    accepted = accepted_records()
    check_ff_release_pair(accepted)
    common = shared_pins()
    check_canoes_pins(common)
    pins = batch_pins(args["stage"], common)
    work = Path(args["work_dir"])
    if stage.pipeline == "order0":
        targets = [build_order0_target(name, work) for name in stage.targets]
    else:
        targets = [build_fold_target(name, work, args, accepted, pins) for name in stage.targets]
    return {"schema": 1, "kind": "zs3p0 T-001 Order-0/FF/main-FK replay batch", "scope": SCOPE,
            "status": "CANDIDATE_NOT_SCIENTIFICALLY_ACCEPTED", "args": args,
            "stage": {"name": args["stage"], "pipeline": stage.pipeline, "role": stage.role},
            "driver": record(DRIVER), "stage_test": bool(args["stage_test"]), "runtime": runtime,
            "framework": framework, "canoes": canoes,
            "order0_runtime": order0_runtime() if stage.pipeline == "order0" else None,
            "reference_pins": [record(p) for p in REFERENCE_FILES], "shared_pins": common, "pins": pins,
            "targets": targets}


def generated_records(target: dict) -> list[dict]:
    root = Path(target["root"])
    if target["pipeline"] == "order0":
        return [{"path": str(root / "job.json"), "sha256": target["job_sha256"]}]
    files = target["files"]
    out = [{"path": str(root / files["contract"]), "sha256": target["contract_sha256"]},
           {"path": str(root / files["config"]), "sha256": target["config_sha256"]}]
    if files["loader"] is not None:
        out.append({"path": str(root / files["loader"]), "sha256": target["loader_sha256"]})
    return sorted(out, key=lambda item: item["path"])


def loader_self_check(target: dict) -> str:
    """Import the FK loader in the fold environment: its guard must install (no vertex query)."""
    root = Path(target["root"])
    code = ("import importlib.util as u, sys\n"
            "s = u.spec_from_file_location('zs3v2_loader_self_check', sys.argv[1])\n"
            "m = u.module_from_spec(s)\ns.loader.exec_module(m)\n"
            "assert m.ALLOW_INTERPOLATION is False and m._CACHE is None\n"
            "print('GUARD_INSTALLED', m.coupling_fn_batch.__name__, str(m.TABLE_PATH))\n")
    done = subprocess.run([str(PYTHON), "-B", "-c", code, str(root / target["files"]["loader"])],
                          env=fold_environment(root), cwd=root, capture_output=True, text=True, check=False)
    expected = f"GUARD_INSTALLED batch {REFINED}"
    if done.returncode != 0 or done.stdout.strip() != expected:
        raise Refusal(f"Loader guard self-check failed for {target['id']}: {done.stdout}{done.stderr[-2000:]}")
    return done.stdout.strip()


def c_module_self_check(target: dict) -> str:
    """Import the C23 wrapper in the fold environment: its table hash check must pass (no fold)."""
    root = Path(target["root"])
    code = ("import importlib.util as u, sys\n"
            "s = u.spec_from_file_location('zs3v2_c_self_check', sys.argv[1])\n"
            "m = u.module_from_spec(s)\ns.loader.exec_module(m)\n"
            "assert m.C_fn_batch.vectorized is True and m.C_fn_batch.source_coord == 'lambda_project_mpc'\n"
            "print('C_LOADED', m.C_fn_batch.table_sha256)\n")
    done = subprocess.run([str(PYTHON), "-B", "-c", code, str(C_WRAPPER)], env=fold_environment(root), cwd=root,
                          capture_output=True, text=True, check=False)
    expected = f"C_LOADED {PINNED[C23_TABLE]}"
    if done.returncode != 0 or done.stdout.strip() != expected:
        raise Refusal(f"C wrapper self-check failed for {target['id']}: {done.stdout}{done.stderr[-2000:]}")
    return done.stdout.strip()


def order0_self_check(target: dict, plan: dict) -> str:
    """The PyCCL interpreter in the job directory matches the plan; the script bytes are pinned."""
    if order0_runtime(Path(target["root"])) != plan["order0_runtime"] or digest(O0_SCRIPT) != PINNED[O0_SCRIPT]:
        raise Refusal(f"Order-0 runtime self-check failed for {target['id']}")
    return f"ORDER0_RUNTIME_OK {plan['order0_runtime']['dependencies']}"


def refuse_production_location(target: dict, work: Path) -> None:
    """No output beside an accepted config or under a PRODUCTION-marked folder; checked before mkdir."""
    out = Path(target["result"] if target["pipeline"] != "order0" else target["out"])
    if not inside(out, work) or out.parent.parent != work:
        raise Refusal(f"Output {out} is not in its own target root inside {work}")
    markers = [folder / "PRODUCTION" for folder in (out.parent, work, HERE)]
    protected = {MAIN_EXEC.resolve(), R1A.resolve(), ACCEPTANCE.resolve(), V1.resolve()}
    if out.parent.resolve() in protected or any(m.exists() for m in markers):
        raise Refusal(f"Output would sit in a production or accepted directory: {out}")


# ---------------------------------------------------------------- prepare
def n_jobs_arg(text: str) -> int | str:
    """argparse type: a positive ASCII integer, or the literal 'accepted'."""
    if text == "accepted":
        return text
    if not re.fullmatch(r"[1-9][0-9]*", text):
        raise argparse.ArgumentTypeError("must be a positive integer or 'accepted'")
    return int(text)


def materialise_batch(work: Path, plan: dict) -> dict:
    """Write each target's records; check them, the fold import binding and the loader or C wrapper."""
    checks = {}
    for target in plan["targets"]:
        root = Path(target["root"])
        root.mkdir(exist_ok=False)
        if target["pipeline"] == "order0":
            write_new_bytes(root / "job.json", canonical(target["job"]))
            check_records(generated_records(target))
            checks[target["id"]] = order0_self_check(target, plan)
            continue
        files = target["files"]
        write_new_bytes(root / files["contract"], canonical(target["contract"]))
        if files["loader"] is not None:
            write_new_text(root / files["loader"], target["loader_text"])
        write_new_text(root / files["config"], canonical(target["config"]).decode() + "\n")
        check_records(generated_records(target))
        written = yaml.safe_load((root / files["config"]).read_text(encoding="utf-8"))
        try:
            replay().validate_config_paths(written, work)
        except (RuntimeError, OSError) as exc:
            raise Refusal(f"{target['id']} config fails replay_suite.validate_config_paths: {exc}") from exc
        subprocess_binding(fold_environment(root), root)
        checks[target["id"]] = loader_self_check(target) if target["pipeline"] == "fk" else c_module_self_check(target)
    return checks


def record_prepare_failure(work: Path, exc: BaseException) -> None:
    """Mark a batch directory whose preparation failed; it never becomes runnable."""
    note = {"schema": 1, "status": "PREPARE_FAILED", "error": repr(exc),
            "note": "Not a batch: there is no prepared.json. Prepare again with a new --work-dir."}
    try:
        write_new_json(work / "prepare_failed.json", note)
    except (OSError, Refusal) as marker_error:
        print(f"WARNING: could not write {work / 'prepare_failed.json'}: {marker_error}", file=sys.stderr)


def stage_args(ns: argparse.Namespace) -> dict:
    stage = STAGES[ns.stage]
    if stage.pipeline == "order0" and ns.n_jobs != "accepted":
        raise Refusal("Order-0 jobs are single processes; --n-jobs applies only to folds")
    work = new_dir_arg(ns.work_dir, stage.default_dir)
    return {"stage": ns.stage, "work_dir": str(work), "n_jobs": ns.n_jobs,
            "stage_test": bool(getattr(ns, "allow_stage_test", False))}


def plan_command(ns: argparse.Namespace) -> int:
    """Print a compact view of one stage's plan; writes nothing."""
    require_invocation()
    plan = build_batch(stage_args(ns))
    print(json.dumps({"stage": plan["stage"], "plan_sha256": sha_bytes(canonical(plan)),
                      "targets": [summarise_target(t) for t in plan["targets"]]}, indent=2))
    return 0


def summarise_target(target: dict) -> dict:
    view = {"id": target["id"], "z": target["z"], "lambda_mpc": target["source"]["lambda_mpc"]}
    if target["pipeline"] == "order0":
        view.update(argv=target["argv"], settings_diff=[(r["field"], r["category"])
                                                        for r in target["accepted_diff"]["settings"]])
    else:
        view.update(n_jobs=target["n_jobs"], n_jobs_source=target["n_jobs_source"], gl=target["gl"],
                    rows=target["expected_rows"], contract_sha256=target["contract_sha256"],
                    diffs={k: [(r["field"], r["category"]) for r in v["config"]]
                           for k, v in target["accepted_diff"].items()})
    return view


def prepare_command(ns: argparse.Namespace) -> int:
    """Build the plan in memory, then materialise it into a new --work-dir and write prepared.json."""
    require_invocation()
    if HERE != PLANNED_HOME and not ns.allow_stage_test:
        raise Refusal(f"The driver must sit at {PLANNED_HOME} (use --allow-stage-test only for test prepares)")
    args = stage_args(ns)
    work = Path(args["work_dir"])
    if not inside(work, HERE) or work.exists() or not work.parent.is_dir():
        raise Refusal(f"--work-dir must be a new directory inside {HERE} with an existing parent: {work}")
    plan = build_batch(args)
    for target in plan["targets"]:
        refuse_production_location(target, work)
    work.mkdir(exist_ok=False)
    try:
        checks = materialise_batch(work, plan)
        prepared = {"schema": 1, "status": "PREPARED_NOT_STARTED", "scope": SCOPE, "prepared_unix": time.time(),
                    "plan_sha256": sha_bytes(canonical(plan)), "self_check": checks,
                    "generated": {t["id"]: generated_records(t) for t in plan["targets"]}, "plan": plan}
        write_new_json(work / "prepared.json", prepared)
    except BaseException as exc:
        record_prepare_failure(work, exc)
        raise
    print(json.dumps({"prepared": str(work / "prepared.json"), "plan_sha256": prepared["plan_sha256"],
                      "stage": plan["stage"], "stage_test": plan["stage_test"],
                      "targets": [summarise_target(t) for t in plan["targets"]]}, indent=2))
    return 0


def load_prepared(work: Path) -> tuple[dict, dict]:
    """prepared.json of a runnable batch: intact plan hash, this driver, its own directory, not a stage test."""
    work = Path(os.path.abspath(work)).resolve()
    path = work / "prepared.json"
    if not path.is_file():
        raise Refusal(f"No prepared.json in {work}; a failed or missing prepare needs a new --work-dir")
    prepared = read_json(path)
    plan = prepared["plan"]
    if sha_bytes(canonical(plan)) != prepared["plan_sha256"]:
        raise Refusal("prepared.json plan hash mismatch")
    if plan["stage_test"] or HERE != PLANNED_HOME:
        raise Refusal("Stage-test batches, or a driver outside its planned directory, cannot run, gate or export")
    if plan["driver"]["sha256"] != digest(DRIVER):
        raise Refusal("The driver changed after preparation")
    if Path(plan["args"]["work_dir"]) != work:
        raise Refusal(f"prepared.json was prepared for {plan['args']['work_dir']}, not {work}")
    if plan["stage"]["name"] not in STAGES:
        raise Refusal(f"Unknown stage {plan['stage']['name']}")
    return prepared, plan


def require_current(plan: dict, prepared: dict) -> None:
    """Rebuild the plan from current inputs; any difference refuses."""
    if canonical(build_batch(plan["args"])) != canonical(plan):
        raise Refusal("Inputs, runtime or records changed since preparation")
    for records in prepared["generated"].values():
        check_records(records)


def require_stage(plan: dict, name: str) -> None:
    if plan["stage"]["name"] != name:
        raise Refusal(f"Expected a {name} batch, got {plan['stage']['name']}")


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


def memory_admission(key: tuple[str, int], min_gib: float | None) -> dict:
    """Free + inactive memory must cover the recorded peak for this pipeline and worker count plus 12 GiB."""
    if min_gib is None:
        if key not in MEASURED_PEAK_BYTES:
            raise Refusal(f"No recorded peak for {key}; pass --min-available-gib")
        required = MEASURED_PEAK_BYTES[key] + MEMORY_FLOOR_BYTES
    else:
        required = int(min_gib * 2**30)
    available = available_memory_bytes()
    if available < required:
        raise Refusal(f"Available memory {available / 2**30:.1f} GiB < required {required / 2**30:.1f} GiB")
    return {"method": "vm_stat free+inactive, observation not reservation", "available_bytes": available,
            "required_bytes": required, "key": list(key),
            "rule": "explicit" if min_gib is not None else "recorded peak + 12 GiB"}


def target_outputs(target: dict) -> list[Path]:
    root = Path(target["root"])
    if target["pipeline"] == "order0":
        return [*map(Path, target["outputs"]), Path(target["log"])]
    return [Path(target["result"]), root / "cache_expand", root / "cache_propagators"]


def refuse_existing(target: dict) -> None:
    root = Path(target["root"])
    names = [root / n for n in ("started.json", "completed.json", "failed.json")] + target_outputs(target)
    if any(p.exists() or p.is_symlink() for p in names):
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


def process_table() -> str:
    return subprocess.check_output(["ps", "-axo", "pid=,command="], text=True)


def require_no_active_folds() -> None:
    """The reviewed ps check plus Order-0 and suite drivers; ps being unavailable (sandbox) refuses."""
    try:
        text = process_table()
        shared().check_active_folds(text)
    except (RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        raise Refusal(f"Active-fold check refused: {exc}") from exc
    active = [line for line in text.splitlines() if any(name in line for name in EXTRA_ACTIVE)]
    if active:
        raise Refusal("Active-fold check refused: another job is active:\n" + "\n".join(active))


def gate_report_path(ns: argparse.Namespace, plan: dict) -> Path:
    if getattr(ns, "gate_report", None) is not None:
        return Path(ns.gate_report)
    return HERE / STAGES[GATE_STAGE[plan["stage"]["pipeline"]]].default_dir / "gate_report.json"


def run_command(ns: argparse.Namespace) -> int:
    """One fold or one Order-0 job. A production target first needs its pipeline's PASS gate report."""
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    target = next((t for t in plan["targets"] if t["id"] == ns.target), None)
    if target is None:
        raise Refusal(f"Unknown target {ns.target}")
    refuse_existing(target)
    gate = None if plan["stage"]["role"] == "gate" else gate_evidence(gate_report_path(ns, plan), plan)
    with shared_fold_lock():
        require_no_active_folds()
        require_current(plan, prepared)
        refuse_existing(target)
        root = Path(target["root"])
        if target["pipeline"] == "order0":
            memory = memory_admission(("order0", 1), ns.min_available_gib)
            order0_self_check(target, plan)
            execute_order0(work, prepared, plan, target, memory, gate)
            return 0
        memory = memory_admission((target["pipeline"], target["n_jobs"]), ns.min_available_gib)
        env = fold_environment(root)
        subprocess_binding(env, root)
        loader_self_check(target) if target["pipeline"] == "fk" else c_module_self_check(target)
        execute_fold(work, prepared, plan, target, env, memory, gate)
    return 0


def started_record(work: Path, prepared: dict, target: dict, memory: dict, gate: dict | None, env: dict) -> dict:
    return {"schema": 1, "status": "STARTED_NOT_COMPLETED", "target_id": target["id"], "pid": os.getpid(),
            "start_unix": time.time(), "plan_sha256": prepared["plan_sha256"],
            "prepared": record(work / "prepared.json"),
            "lock": str(LOCK), "memory": memory, "gate": gate,
            "environment": {k: env[k] for k in sorted(env) if k.startswith(("PYTHON", "CANOES", "MPL", "OMP",
                                                                             "OPENBLAS", "MKL", "VECLIB"))}}


def execute_fold(work: Path, prepared: dict, plan: dict, target: dict, env: dict, memory: dict,
                 gate: dict | None = None) -> None:
    """started.json, run_fk_variant.run in the fold environment, result validation, completed.json."""
    root, result = Path(target["root"]), Path(target["result"])
    config = root / target["files"]["config"]
    fold = variant().run
    write_new_json(root / "started.json", {**started_record(work, prepared, target, memory, gate, env),
                                           "n_jobs": target["n_jobs"], "gl": target["gl"], "orders": ORDERS})
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


def launch_order0(target: dict, env: dict, log_handle: Any) -> int:
    """Run source_folded_reference.py with the prepared argv; returns its exit code."""
    done = subprocess.run(target["argv"], env=env, cwd=target["root"], stdout=log_handle,
                          stderr=subprocess.STDOUT, timeout=ORDER0_TIMEOUT_S, check=False)
    return done.returncode


def execute_order0(work: Path, prepared: dict, plan: dict, target: dict, memory: dict,
                   gate: dict | None = None) -> None:
    """started.json, one source_folded_reference.py process, output validation, completed.json."""
    root = Path(target["root"])
    env = order0_environment(root)
    write_new_json(root / "started.json", {**started_record(work, prepared, target, memory, gate, env),
                                           "argv": target["argv"]})
    started = time.monotonic()
    try:
        with Path(target["log"]).open("xb") as log_handle:
            code = launch_order0(target, env, log_handle)
        if code != 0:
            raise Refusal(f"source_folded_reference.py exited with {code}; see {target['log']}")
        check_records(plan["pins"] + plan["reference_pins"] + prepared["generated"][target["id"]])
        validation = validate_order0_output(target)
    except BaseException as exc:
        write_new_json(root / "failed.json", {"schema": 1, "status": "failed", "target_id": target["id"],
                                              "error": repr(exc), "elapsed_seconds": time.monotonic() - started})
        raise
    write_new_json(root / "completed.json", {
        "schema": 1, "status": "completed", "target_id": target["id"], "plan_sha256": prepared["plan_sha256"],
        "prepared": record(work / "prepared.json"), "started": record(root / "started.json"),
        "elapsed_seconds": time.monotonic() - started, "returncode": 0,
        "outputs": [record(Path(p)) for p in target["outputs"]], "log": record(Path(target["log"])),
        "output_validation": validation, "scientific_acceptance": False})
    print(json.dumps({"completed": target["id"], "out": record(Path(target["out"]))}, indent=2))


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
    """Reviewed validate_result plus the (0,0) order-2 endpoint and direction check."""
    try:
        schema = shared().validate_result(result, target["config"])
        rows = read_rows(result)
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        raise Refusal(f"Result {result} fails the reviewed result validation: {exc}") from exc
    grid = target["config"]["sweep"]["positions_grid"]
    lam, count = target["config"]["sweep"]["t_final_grid"][0], target["expected_rows"]
    ok = (rows["value"].shape == (count,) and np.all(rows["a"] == 0) and np.all(rows["b"] == 0)
          and np.all(rows["order"] == 2) and np.all(rows["t_final"] == lam)
          and np.all(rows["x"] == np.asarray(grid["x"][0], float))
          and sorted(map(tuple, rows["y"])) == sorted(map(tuple, np.asarray(grid["y"], float))))
    if not ok:
        raise Refusal(f"Result {result} fails the {count}-row (0,0) order-2 endpoint check")
    return {**schema, "rows": count, "t_final": lam, "pair": [0, 0],
            "value_sha256": sha_bytes(np.ascontiguousarray(rows["value"]).tobytes())}


def completion_of(target: dict, plan_sha: str) -> dict:
    """completed.json of a target whose completion belongs to this plan and has no failure record."""
    root = Path(target["root"])
    if (root / "failed.json").exists():
        raise Refusal(f"{target['id']} has failed.json")
    if not (root / "completed.json").is_file():
        raise Refusal(f"{target['id']} has not completed (no completed.json in {root})")
    completed = read_json(root / "completed.json")
    if completed["target_id"] != target["id"] or completed["plan_sha256"] != plan_sha:
        raise Refusal(f"{target['id']} completion belongs to another plan")
    check_records([completed["started"]])
    return completed


def completed_rows(target: dict, plan_sha: str) -> tuple[dict, dict]:
    """Rows of a completed fold whose completion, started record and result are all intact."""
    completed = completion_of(target, plan_sha)
    check_records([completed["result"]])
    if completed["result"]["path"] != str(Path(target["result"]).resolve()):
        raise Refusal(f"{target['id']} completion points at another result")
    if validate_fold_result(Path(target["result"]), target) != completed["output_validation"]:
        raise Refusal(f"{target['id']} result no longer matches its completion record")
    root = Path(target["root"])
    return read_rows(Path(target["result"])), {"completed": record(root / "completed.json"),
                                               "started": completed["started"], "result": completed["result"],
                                               "elapsed_seconds": completed["elapsed_seconds"]}


def load_arrays(path: Path, keys: list[str] | tuple[str, ...]) -> dict:
    with np.load(path, allow_pickle=False) as data:
        missing = [k for k in keys if k not in data.files]
        if missing:
            raise Refusal(f"{path} lacks arrays {missing}")
        return {k: np.asarray(data[k]) for k in keys}


def validate_order0_output(target: dict) -> dict:
    """Outputs present; settings, endpoint, script and power hashes as planned; no warnings; finite arrays."""
    out = Path(target["out"])
    for path in map(Path, target["outputs"]):
        if not path.is_file():
            raise Refusal(f"Order-0 output missing: {path}")
    if digest(out.with_suffix(".source.py")) != PINNED[O0_SCRIPT]:
        raise Refusal("Order-0 script snapshot differs from the pinned script")
    with np.load(out, allow_pickle=False) as data:
        meta = json.loads(str(data["metadata_json"]))
    arrays = load_arrays(out, target["arrays"])
    expected = target["expected_settings"]
    problems = []
    if read_json(out.with_suffix(".json")) != meta:
        problems.append("JSON sidecar differs from the embedded metadata")
    if meta.get("settings") != expected:
        problems.append("recorded settings differ from the planned settings")
    for key, want in (("warnings", []), ("observable_operator", "full"), ("script_sha256", PINNED[O0_SCRIPT]),
                      ("power_table_sha256", PINNED[PCAMB]), ("source_lambda_mpc", target["source"]["lambda_mpc"]),
                      ("source_redshift", target["z"]), ("source_map", expected["source_map"])):
        if meta.get(key) != want:
            problems.append(f"metadata {key}={meta.get(key)!r}, expected {want!r}")
    problems += [f"non-finite {k}" for k, v in arrays.items() if not np.all(np.isfinite(v))]
    if problems:
        raise Refusal(f"Order-0 output {out} fails validation: {problems}")
    return {"settings_match": True, "warnings": [], "source_lambda_mpc": meta["source_lambda_mpc"],
            "source_redshift": meta["source_redshift"], "metadata_sha256": sha_bytes(canonical(meta)),
            "array_sha256": {k: sha_bytes(np.ascontiguousarray(v).tobytes()) for k, v in sorted(arrays.items())}}


def completed_order0(target: dict, plan_sha: str) -> tuple[dict, dict]:
    """Arrays of a completed Order-0 job whose completion, outputs and log are intact."""
    completed = completion_of(target, plan_sha)
    check_records(completed["outputs"] + [completed["log"]])
    if [item["path"] for item in completed["outputs"]] != [str(Path(p).resolve()) for p in target["outputs"]]:
        raise Refusal(f"{target['id']} completion points at other outputs")
    if validate_order0_output(target) != completed["output_validation"]:
        raise Refusal(f"{target['id']} outputs no longer match their completion record")
    root = Path(target["root"])
    return load_arrays(Path(target["out"]), target["arrays"]), {
        "completed": record(root / "completed.json"), "started": completed["started"],
        "outputs": completed["outputs"], "elapsed_seconds": completed["elapsed_seconds"]}


# ------------------------------------------------------------- comparisons
def row_keys(rows: dict) -> list[tuple]:
    return [(tuple(rows["x"][i]), tuple(rows["y"][i]), float(rows["t_final"][i]), int(rows["a"][i]),
             int(rows["b"][i]), int(rows["order"][i])) for i in range(rows["value"].shape[0])]


def compare_rows(mine: dict, ref: dict, label: str, reference: Path) -> dict:
    """Row-matched comparison: max relative difference, bitwise flag, pass at REL_TOL."""
    keys, ref_keys = row_keys(mine), row_keys(ref)
    if len(set(keys)) != len(keys) or sorted(keys) != sorted(ref_keys):
        return {"label": label, "reference": record(reference), "pass": False, "bitwise": False,
                "reason": "row sets differ"}
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


def compare_arrays(mine: np.ndarray, ref: np.ndarray) -> dict:
    """Bitwise flag, max |difference|, max |difference| over the reference peak, max pointwise relative."""
    a, b = np.asarray(mine, float), np.asarray(ref, float)
    if a.shape != b.shape:
        return {"shape_equal": False, "bitwise": False, "max_abs_over_peak": math.inf}
    diff, scale = np.abs(a - b), np.abs(b)
    peak, worst = float(np.max(scale)), float(np.max(diff))
    rel = np.where(scale > 0, diff / np.where(scale > 0, scale, 1.0), np.where(diff > 0, np.inf, 0.0))
    return {"shape_equal": True, "shape": list(a.shape), "bitwise": a.tobytes() == b.tobytes(),
            "max_abs_diff": worst, "max_abs_over_peak": worst / peak if peak > 0 else (0.0 if worst == 0 else math.inf),
            "max_rel_pointwise": float(np.max(rel))}


def separation_arcmin(x: Any, y: Any) -> float:
    """nhfp_fig_evidence.separation_arcmin."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    cosang = np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y))
    return float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))) * 60.0)


def main_pair_rows(path: Path) -> dict:
    """Pair (0,0) order-2 rows of an accepted main product, in main-grid order."""
    with np.load(path, allow_pickle=False) as data:
        sel = (data["a"] == 0) & (data["b"] == 0) & (data["order"] == 2)
        rows = {k: np.asarray(data[k])[sel] for k in ("x", "y", "t_final", "a", "b", "order", "value")}
    if rows["value"].shape != (N_MAIN,) or not np.array_equal(rows["y"], np.asarray(main_grid()[1], float)):
        raise Refusal(f"{path} does not hold the 40 (0,0) rows in main-grid order")
    return rows


def reference_rows(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {k: np.asarray(data[k]) for k in ("x", "y", "t_final", "a", "b", "order", "value")}


def subset_rows(rows: dict, indices: tuple[int, ...] | list[int]) -> dict:
    return {k: v[list(indices)] for k, v in rows.items()}


def multiz_rows(indices: tuple[int, ...]) -> dict:
    """Accepted multiz tree FK at z=3.2 (row 3) as (0,0) order-2 rows at main indices."""
    with np.load(ACCEPTED / "products/multiz.npz", allow_pickle=False) as data:
        multiz = {k: np.asarray(data[k]) for k in ("z", "lam", "gamma", "fk")}
    x, ys = main_grid()
    sep = np.array([separation_arcmin(x, y) for y in ys])
    if (multiz["z"][3] != Z_REGRESSION or multiz["lam"][3] != EXPECTED_LAMBDA[Z_REGRESSION]
            or np.max(np.abs(multiz["gamma"] - sep)) > NODE_ATOL_ARCMIN):
        raise Refusal("Accepted multiz row 3 is not z=3.2 on the main grid")
    idx = list(indices)
    n = len(idx)
    return {"x": np.tile(np.asarray(x, float), (n, 1)), "y": np.asarray(ys, float)[idx],
            "t_final": np.full(n, EXPECTED_LAMBDA[Z_REGRESSION]), "a": np.zeros(n, np.int64),
            "b": np.zeros(n, np.int64), "order": np.full(n, 2, np.int64), "value": multiz["fk"][3, idx]}


# -------------------------------------------------------------------- gates
def order0_gate_checks(plan: dict, prepared: dict) -> tuple[dict, dict]:
    """O1: z=3.2 recorded argv against source_folded_multiz_z3p2; O2: z=5 full channels against the z=5 source."""
    by_name = {t["name"]: t for t in plan["targets"]}
    with np.load(ACCEPTED / "products/multiz.npz", allow_pickle=False) as data:
        o0_row3, z_row3 = np.asarray(data["o0"][3]), float(data["z"][3])
    accepted_kk = load_arrays(O0_Z3P2, ["kk"])["kk"]
    precondition = {"accepted_z3p2_kk_equals_multiz_o0_row3_bitwise": accepted_kk.tobytes() == o0_row3.tobytes(),
                    "multiz_row3_z": z_row3}
    checks, completions = {}, {}
    for name, (reference, exact, tolerance, report) in ORDER0_GATES.items():
        arrays, completion = completed_order0(by_name[name], prepared["plan_sha256"])
        ref = load_arrays(reference, [*exact, *tolerance, *report])
        rows = {k: compare_arrays(arrays[k], ref[k]) for k in (*exact, *tolerance, *report)}
        passed = (all(rows[k]["bitwise"] for k in exact)
                  and all(rows[k]["max_abs_over_peak"] <= ORDER0_PEAK_TOL for k in tolerance))
        if name == "o0_parity_zs3p2":
            passed = passed and precondition["accepted_z3p2_kk_equals_multiz_o0_row3_bitwise"] and z_row3 == 3.2
        checks[name] = [{"label": f"{name} against {reference.name}", "reference": record(reference),
                         "exact_arrays": list(exact), "tolerance_arrays": list(tolerance), "arrays": rows,
                         "bitwise": all(r["bitwise"] for r in rows.values()), "pass": bool(passed)}]
        completions[name] = completion
    checks["precondition"] = [{"label": "accepted z=3.2 Order-0 equals multiz o0[3]", **precondition,
                               "bitwise": precondition["accepted_z3p2_kk_equals_multiz_o0_row3_bitwise"],
                               "pass": bool(precondition["accepted_z3p2_kk_equals_multiz_o0_row3_bitwise"])}]
    return checks, completions


def ff_gate_checks(plan: dict, prepared: dict) -> tuple[dict, dict]:
    """F1: FF canary rows against the accepted main FF pair (0,0) at main indices 0 and 6."""
    target = plan["targets"][0]
    rows, completion = completed_rows(target, prepared["plan_sha256"])
    reference = ACCEPTED / "products/main_ff.npz"
    check = compare_rows(rows, subset_rows(main_pair_rows(reference), CANARY_INDICES),
                         "accepted main FF (0,0) at main indices [0, 6]", reference)
    return {"F1_ff_canary_zs5": [check]}, {target["name"]: completion}


def fk_gate_checks(plan: dict, prepared: dict) -> tuple[dict, dict]:
    """K1: z=3.2 regression rows against the accepted multiz fk[3] at main indices 0,5,10,15,25,35."""
    target = plan["targets"][0]
    rows, completion = completed_rows(target, prepared["plan_sha256"])
    check = compare_rows(rows, multiz_rows(REGRESSION_INDICES), "accepted multiz fk[3] (z=3.2) at main indices "
                         "[0, 5, 10, 15, 25, 35]", ACCEPTED / "products/multiz.npz")
    return {"K1_fk_regression_zs3p2": [check]}, {target["name"]: completion}


def evaluate_gate(work: Path, prepared: dict, plan: dict) -> dict:
    """Gate report of a gate batch from its completed targets; reads only."""
    if plan["stage"]["role"] != "gate":
        raise Refusal("Gates take only a gate batch (order0_parity, ff_canary or fk_regression)")
    pipeline = plan["stage"]["pipeline"]
    checks, completions = {"order0": order0_gate_checks, "ff": ff_gate_checks, "fk": fk_gate_checks}[pipeline](
        plan, prepared)
    passed = all(c["pass"] for group in checks.values() for c in group)
    criterion = (f"bitwise grids; max |difference|/peak <= {ORDER0_PEAK_TOL:g} for kk (and xip, xim, kgt)"
                 if pipeline == "order0" else f"max |value/reference - 1| <= {REL_TOL:g}")
    workers = ({t["id"]: {"used": t["n_jobs"], "source": t["n_jobs_source"], "gl": t["gl"]} for t in plan["targets"]}
               if pipeline != "order0" else {t["id"]: {"processes": 1} for t in plan["targets"]})
    return {"schema": 1, "status": "PASS" if passed else "FAIL", "pipeline": pipeline, "criterion": criterion,
            "bitwise_all": all(c.get("bitwise", False) for group in checks.values() for c in group),
            "prepared": record(work / "prepared.json"), "plan_sha256": prepared["plan_sha256"],
            "completions": completions, "workers": workers, "checks": checks, "scientific_acceptance": False}


def report_summary(report: dict) -> dict:
    keys = ("label", "max_rel_diff", "bitwise", "pass")
    return {"status": report["status"], "bitwise_all": report.get("bitwise_all"),
            "checks": {k: [{x: c.get(x) for x in keys} for c in v] for k, v in report["checks"].items()}}


def gate_command(ns: argparse.Namespace) -> int:
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    out = work / "gate_report.json"
    if out.exists() or out.is_symlink():
        raise Refusal(f"{out} exists; a gate report is never replaced (a failed gate needs a new gate batch)")
    require_current(plan, prepared)
    report = evaluate_gate(work, prepared, plan)
    write_new_json(out, report)
    print(json.dumps(report_summary(report), indent=2))
    return 0 if report["status"] == "PASS" else 1


def bound_report(path: Path, name: str) -> tuple[dict, Path, dict, dict]:
    """A PASS report that sits in the batch it reports on, with that batch current and intact."""
    path = Path(os.path.abspath(path)).resolve()
    if not path.is_file():
        raise Refusal(f"No {name} at {path}")
    report = read_json(path)
    if not isinstance(report, dict) or report.get("status") != "PASS":
        raise Refusal(f"{name} {path} is not PASS")
    prepared_record = report.get("prepared")
    if not isinstance(prepared_record, dict) or not isinstance(prepared_record.get("path"), str):
        raise Refusal(f"{name} {path} has no prepared-batch record")
    work = Path(prepared_record["path"]).parent
    if path.parent != work:
        raise Refusal(f"{name} {path} does not sit in the batch it reports on ({work})")
    prepared, plan = load_prepared(work)
    if report["prepared"] != record(work / "prepared.json") or report["plan_sha256"] != prepared["plan_sha256"]:
        raise Refusal(f"{name} belongs to another prepared batch")
    require_current(plan, prepared)
    return report, work, prepared, plan


def gate_evidence(gate_report: Path, plan: dict) -> dict:
    """Re-verify a PASS gate report from its own batch and bind it to ``plan``'s inputs."""
    pipeline = plan["stage"]["pipeline"]
    report, gate_work, gate_prepared, gate_plan = bound_report(gate_report, "gate report")
    if gate_plan["stage"]["name"] != GATE_STAGE[pipeline]:
        raise Refusal(f"Gate batch is {gate_plan['stage']['name']}, not {GATE_STAGE[pipeline]}")
    keys = ("driver", "runtime", "framework", "canoes", "order0_runtime", "reference_pins", "shared_pins")
    if any(gate_plan[k] != plan[k] for k in keys):
        raise Refusal("Gate batch ran with another driver, runtime, framework, canoes state or shared inputs")
    if json.loads(json.dumps(evaluate_gate(gate_work, gate_prepared, gate_plan))) != report:
        raise Refusal(f"Gate report {gate_report} differs from a recomputation from its batch")
    if pipeline != "order0":
        reference = gate_plan["targets"][0]
        for target in plan["targets"]:
            if any(reference[k] != target[k] for k in GATE_MATCH_KEYS):
                raise Refusal(f"{target['id']} differs from gate target {reference['id']} in {GATE_MATCH_KEYS}")
    return {"gate_report": record(Path(os.path.abspath(gate_report)).resolve()), "status": report["status"],
            "bitwise_all": report["bitwise_all"], "workers": report["workers"],
            "gate_plan_sha256": gate_prepared["plan_sha256"]}


# ------------------------------------------------------------- consistency
def v1_tree_rows() -> dict:
    """v1 z_s=3.0 tree_main_15360 rows, tied to its completion record and assembled ladder member."""
    completed = read_json(V1_TREE / "completed.json")
    if (completed.get("status") != "completed" or completed["result"] != record(V1_TREE_RESULT)
            or completed["result"]["sha256"] != PINNED[V1_TREE_RESULT]):
        raise Refusal("v1 tree_main_15360 completion does not record the pinned result")
    rows = read_rows(V1_TREE_RESULT)
    member = reference_rows(V1_LADDER / LADDER_TOP)
    if not all(np.array_equal(np.asarray(rows[k]), np.asarray(member[k])) for k in rows):
        raise Refusal("v1 assembled tree l_max=15360 member differs from the v1 raw result")
    if not np.all(rows["t_final"] == EXPECTED_LAMBDA[Z_PRODUCTION]) or rows["value"].shape != (len(LADDER_INDICES),):
        raise Refusal("v1 tree_main_15360 is not the 11-row z_s=3.0 result")
    return rows


def evaluate_consistency(work: Path, prepared: dict, plan: dict) -> dict:
    """C1: z_s=3.0 main FK at main indices [0,5..14] against the v1 tree_main_15360 fold."""
    require_stage(plan, "fk_production")
    target = plan["targets"][0]
    rows, completion = completed_rows(target, prepared["plan_sha256"])
    ys = [tuple(y) for y in np.asarray(main_grid()[1], float)]
    mine = [tuple(y) for y in rows["y"]]
    if any(ys[i] not in mine for i in LADDER_INDICES):
        raise Refusal("The z_s=3.0 main FK lacks a ladder direction")
    order = [mine.index(ys[i]) for i in LADDER_INDICES]
    check = compare_rows(subset_rows(rows, order), v1_tree_rows(),
                         "v1 z_s=3.0 tree_main_15360 fold at main indices [0, 5..14]", V1_TREE_RESULT)
    return {"schema": 1, "status": "PASS" if check["pass"] else "FAIL",
            "criterion": f"max |value/v1 - 1| <= {REL_TOL:g}",
            "bitwise": check["bitwise"], "prepared": record(work / "prepared.json"),
            "plan_sha256": prepared["plan_sha256"], "completion": completion,
            "v1": {"completed": record(V1_TREE / "completed.json"), "result": record(V1_TREE_RESULT),
                   "ladder_member": record(V1_LADDER / LADDER_TOP)},
            "checks": {"C1_v1_tree_main_15360": [check]}, "scientific_acceptance": False}


def consistency_command(ns: argparse.Namespace) -> int:
    require_invocation()
    prepared, plan = load_prepared(ns.work_dir)
    work = Path(plan["args"]["work_dir"])
    out = work / "consistency_report.json"
    if out.exists() or out.is_symlink():
        raise Refusal(f"{out} exists; a consistency report is never replaced")
    require_current(plan, prepared)
    report = evaluate_consistency(work, prepared, plan)
    write_new_json(out, report)
    print(json.dumps(report_summary(report), indent=2))
    return 0 if report["status"] == "PASS" else 1


def consistency_evidence(path: Path, fk_plan: dict) -> dict:
    report, work, prepared, plan = bound_report(path, "consistency report")
    require_stage(plan, "fk_production")
    if canonical(plan) != canonical(fk_plan):
        raise Refusal("Consistency report belongs to another fk_production batch")
    if json.loads(json.dumps(evaluate_consistency(work, prepared, plan))) != report:
        raise Refusal(f"Consistency report {path} differs from a recomputation from its batch")
    return {"report": record(Path(os.path.abspath(path)).resolve()), "status": report["status"],
            "bitwise": report["bitwise"]}


# -------------------------------------------------------------------- export
def plot_directions(config: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """prepare_plot_inputs.directions: reference direction, directions sorted by separation, separations."""
    positions = yaml.safe_load(config.read_text(encoding="utf-8"))["sweep"]["positions_grid"]
    x = np.asarray(positions["x"], float)
    y = np.asarray(positions["y"], float)
    if x.shape != (1, 3) or y.ndim != 2 or y.shape[1] != 3:
        raise Refusal("Expected one reference direction and a set of three-vector directions")
    unit_x = x[0] / np.linalg.norm(x[0])
    unit_y = y / np.linalg.norm(y, axis=1)[:, None]
    gamma = np.arccos(np.clip(unit_y @ unit_x, -1, 1)) * 10800 / np.pi
    order = np.argsort(gamma)
    if np.any(np.diff(gamma[order]) <= 0):
        raise Refusal("Directions do not give distinct angular separations")
    return x[0], y[order], gamma[order]


def order0_transform(source: Path, source_map: Path, z: float) -> tuple[dict, dict]:
    """prepare_plot_inputs.order0, with the source plane z instead of the hard-coded z=5.

    Returns the payload (x, y, a, b, order, t_final, value, in that order) and the provenance.
    """
    x, y, gamma = plot_directions(O0_PLOT_CONFIG)
    with np.load(source, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        if metadata["observable_operator"] != "full":
            raise Refusal("The manuscript baseline requires the full Ricci source")
        if not np.allclose(gamma, data["theta_arcmin"], rtol=0, atol=1e-7):
            raise Refusal("Source and plot angular grids differ")
        lam = float(metadata["source_lambda_mpc"])
        expected = json.loads(Path(source_map).read_text(encoding="utf-8"))
        match = np.flatnonzero(np.asarray(expected["z"]) == z)
        if match.size != 1 or not np.isclose(lam, expected["lambda_mpc"][int(match[0])], rtol=0, atol=1e-10):
            raise Refusal(f"Order-0 is not on the z={z} source plane of {source_map}")
        kk, plus, minus, cross = [np.asarray(data[k], float) for k in ("kk", "xip", "xim", "kgt")]
    for array in (kk, plus, minus, cross):
        if array.shape != gamma.shape or not np.all(np.isfinite(array)):
            raise Refusal("Invalid source correlation function")
    groups = {(0, 0): kk, (0, 1): -cross, (0, 2): np.zeros_like(kk),
              (1, 1): (plus + minus) / 2, (1, 2): np.zeros_like(kk), (2, 2): (plus - minus) / 2}
    if not np.allclose(groups[(1, 1)] + groups[(2, 2)], plus, rtol=2e-15, atol=1e-20):
        raise Refusal("Spin-trace round trip failed")
    if not np.allclose(groups[(1, 1)] - groups[(2, 2)], minus, rtol=2e-10, atol=1e-19):
        raise Refusal("Spin-difference round trip failed")
    n = len(gamma)
    payload = {"x": np.tile(x, (n * len(groups), 1)), "y": np.tile(y, (len(groups), 1)),
               "a": np.repeat([pair[0] for pair in groups], n), "b": np.repeat([pair[1] for pair in groups], n),
               "order": np.zeros(n * len(groups), dtype=int), "t_final": np.full(n * len(groups), lam),
               "value": np.concatenate(list(groups.values()))}
    provenance = {"kind": "Order-0 plot-format adapter", "source": record(source), "config": record(O0_PLOT_CONFIG),
                  "source_map": record(source_map), "lambda_source_mpc": lam, "z_source": z,
                  "parity_policy": "The scalar Gaussian reference has zero kappa-cross and plus-cross correlations.",
                  "calculation": ("Continuous outer-response projection in the source product. "
                                  "This file is not a framework rerun.")}
    return payload, provenance


def order0_transform_self_test() -> dict:
    """The transformation applied to the accepted z=5 source reproduces figure_products/order0.npz bitwise."""
    payload, _ = order0_transform(O0_Z5, REVIEWED_SOURCE_MAP, 5)
    accepted = load_arrays(O0_FIG, list(payload))
    mismatched = [k for k, v in payload.items()
                  if v.dtype != accepted[k].dtype or v.shape != accepted[k].shape
                  or v.tobytes() != accepted[k].tobytes()]
    if mismatched:
        raise Refusal(f"Order-0 transformation does not reproduce {O0_FIG} for {mismatched}")
    return {"accepted": record(O0_FIG), "source": record(O0_Z5), "arrays_bitwise": sorted(payload),
            "pass": True}


def fold_export_arrays(rows: dict, meta: dict) -> dict:
    """Rows in main-grid order, laid out as the accepted main_fk/main_ff products.

    Members: a, b, metadata_json, order, t_final, value, x, y (x and y float64 (n, 3)).
    """
    ys = [tuple(y) for y in np.asarray(main_grid()[1], float)]
    try:
        index = [ys.index(tuple(y)) for y in rows["y"]]
    except ValueError as exc:
        raise Refusal("A fold direction is not on the accepted main grid") from exc
    if sorted(index) != list(range(N_MAIN)):
        raise Refusal("Fold rows do not cover the 40 main directions exactly once")
    order = np.argsort(index)
    out = {"a": np.asarray(rows["a"])[order].astype(np.int64), "b": np.asarray(rows["b"])[order].astype(np.int64),
           "metadata_json": np.array(canonical(meta).decode()),
           "order": np.asarray(rows["order"])[order].astype(np.int64),
           "t_final": np.asarray(rows["t_final"], float)[order], "value": np.asarray(rows["value"], float)[order],
           "x": np.asarray(rows["x"], float)[order], "y": np.asarray(rows["y"], float)[order]}
    for key in ("a", "b", "order", "t_final", "value"):
        if not np.array_equal(out[key], np.asarray(rows[key])[order]):
            raise Refusal(f"Casting {key} to the accepted dtype changed its values")
    return out


def nhfp_channels(arrays: dict, order: int) -> tuple[np.ndarray, dict]:
    """nhfp_fig_evidence.load_channels on in-memory arrays: separations and per-pair values sorted by separation."""
    a, b, o, v, xs, ys = (np.asarray(arrays[k]) for k in ("a", "b", "order", "value", "x", "y"))
    grid, out = None, {}
    for pair in sorted({(int(i), int(j)) for i, j in zip(a, b, strict=True)}):
        idx = np.flatnonzero((a == pair[0]) & (b == pair[1]) & (o == order))
        if idx.size == 0:
            continue
        sep = np.array([separation_arcmin(xs[i], ys[i]) for i in idx])
        s = np.argsort(sep)
        if grid is None:
            grid = sep[s]
        elif not np.array_equal(grid, sep[s]):
            raise Refusal(f"Separation grid differs for pair {pair}")
        out[pair] = v[idx][s]
    if grid is None:
        raise Refusal(f"No rows at order {order}")
    return grid, out


def ladder_checks(fk_grid: np.ndarray, fk00: np.ndarray) -> dict:
    """nhfp_fig_evidence.load_products ladder checks on the v1 z_s=3.0 ladder.

    The v1 ladder grids must also equal the accepted z_s=5 ladder grids.
    """
    grids, values = {}, {}
    for member in LADDER_MEMBERS:
        for label, root in (("v1_zs3p0", V1_LADDER), ("accepted_zs5", ACCEPTED / "cutoff_ladder")):
            grid, channels = nhfp_channels(reference_rows(root / member), 2)
            grids[f"{label}:{member}"], values[f"{label}:{member}"] = grid, channels[(0, 0)]
    g_lad = grids[f"v1_zs3p0:{LADDER_TOP}"]
    idx = np.array([int(np.argmin(np.abs(fk_grid - t))) for t in g_lad])
    ratio = np.abs(fk00[idx] / values[f"v1_zs3p0:{LADDER_TOP}"] - 1)
    return {"ladder_grids_equal_v1_and_accepted": all(np.array_equal(g, g_lad) for g in grids.values()),
            "ladder_nodes_on_main_grid": bool(np.max(np.abs(fk_grid[idx] - g_lad)) <= NODE_ATOL_ARCMIN),
            "ladder_main_indices": idx.tolist(), "main_indices_expected": list(LADDER_INDICES),
            "main_fk_vs_v1_tree_15360_max_rel": float(np.max(ratio)),
            "main_fk_vs_v1_tree_15360_bitwise": fk00[idx].tobytes() == values[f"v1_zs3p0:{LADDER_TOP}"].tobytes()}


def export_checks(fk: dict, ff: dict, o0: dict, lam: float) -> dict:
    """Exported grids against the accepted z_s=5 grids, the NHFP grid rules and the v1 ladder angles."""
    acc_fk = main_pair_rows(ACCEPTED / "products/main_fk.npz")
    acc_ff = main_pair_rows(ACCEPTED / "products/main_ff.npz")
    acc_o0 = load_arrays(O0_FIG, ["x", "y", "a", "b", "order"])
    g_fk, ch_fk = nhfp_channels(fk, 2)
    g_ff, _ = nhfp_channels(ff, 2)
    g_o0, _ = nhfp_channels(o0, 0)
    g_acc, _ = nhfp_channels(reference_rows(ACCEPTED / "products/main_fk.npz"), 2)
    pole = np.array([0.0, 0.0, 1.0])
    checks = {
        "fk_rows_equal_accepted_main_fk_00_grid": bool(np.array_equal(fk["y"], acc_fk["y"])
                                                       and np.array_equal(fk["x"], acc_fk["x"])),
        "ff_rows_equal_accepted_main_ff_00_grid": bool(np.array_equal(ff["y"], acc_ff["y"])
                                                       and np.array_equal(ff["x"], acc_ff["x"])),
        "order0_layout_equals_accepted_order0": all(np.array_equal(o0[k], acc_o0[k]) for k in acc_o0),
        "pole_reference_direction": bool(np.all(fk["x"] == pole) and np.all(ff["x"] == pole)
                                         and np.all(o0["x"] == pole)),
        "pair_00_order2_only": bool(np.all(fk["a"] == 0) and np.all(fk["b"] == 0) and np.all(fk["order"] == 2)
                                    and np.all(ff["a"] == 0) and np.all(ff["b"] == 0) and np.all(ff["order"] == 2)),
        "t_final_all_lambda_s": bool(np.all(fk["t_final"] == lam) and np.all(ff["t_final"] == lam)
                                     and np.all(o0["t_final"] == lam)),
        "nhfp_grids_fk_ff_order0_equal": bool(np.array_equal(g_fk, g_ff) and np.array_equal(g_fk, g_o0)),
        "nhfp_grid_equals_accepted_zs5": bool(np.array_equal(g_fk, g_acc)),
    }
    ladder = ladder_checks(g_fk, ch_fk[(0, 0)])
    checks.update({k: v for k, v in ladder.items() if isinstance(v, bool)})
    checks["ladder_indices_are_main_0_5_to_14"] = ladder["ladder_main_indices"] == list(LADDER_INDICES)
    checks["main_fk_vs_v1_tree_15360_within_1e-12"] = ladder["main_fk_vs_v1_tree_15360_max_rel"] <= REL_TOL
    failed = [k for k, v in checks.items() if v is not True and k != "main_fk_vs_v1_tree_15360_bitwise"]
    if failed:
        raise Refusal(f"Export grid checks failed: {failed}")
    return {"checks": checks, "ladder": ladder}


def order0_diagnostics(prod: dict, probe: dict, full: dict) -> dict:
    """Structural checks (refusing) and numerical diagnostics (flags only) of the z_s=3.0 Order-0."""
    with np.load(ACCEPTED / "products/multiz.npz", allow_pickle=False) as data:
        gamma, o0 = np.asarray(data["gamma"]), np.asarray(data["o0"])
    taper5 = load_arrays(O0_Z5, ["ell_taper"])["ell_taper"]
    structural = {"production_theta_equals_multiz_gamma": prod["theta_arcmin"].tobytes() == gamma.tobytes(),
                  "full_theta_equals_multiz_gamma": full["theta_arcmin"].tobytes() == gamma.tobytes(),
                  "production_taper_equals_zs5_taper": prod["ell_taper"].tobytes() == taper5.tobytes(),
                  "full_taper_equals_zs5_taper": full["ell_taper"].tobytes() == taper5.tobytes()}
    failed = [k for k, v in structural.items() if not v]
    if failed:
        raise Refusal(f"Order-0 structural checks failed: {failed}")
    # As in source_folded_multiz_validation.json: the generated (interpolated) production
    # spectrum on the full ell grid against the fine probe at its five multipoles.
    ells, rel = [int(v) for v in probe["ell"]], []
    for i, ell in enumerate(ells):
        hit = np.flatnonzero(prod["ell"] == ell)
        if hit.size != 1:
            raise Refusal(f"Production Order-0 has no multipole ell={ell}")
        rel.append(float(prod["cl_pp"][hit[0]] / probe["cl_pp"][i] - 1))
    worst = max(abs(v) for v in rel)
    lo, hi = np.minimum(o0[2], o0[3]), np.maximum(o0[2], o0[3])
    between = (lo < prod["kk"]) & (prod["kk"] < hi)
    resolution = compare_arrays(full["kk"], prod["kk"])
    return {"structural": structural,
            "fine_probe": {"ell": ells, "relative_production_minus_fine_probe": rel, "max_absolute_relative": worst,
                           "flag_above": FINE_PROBE_FLAG, "flagged": worst > FINE_PROBE_FLAG,
                           "historical_multiz_range": [2.3e-4, 3.1e-3]},
            "bracket_multiz_z2p5_z3p2": {"angles_strictly_between": int(between.sum()), "angles": int(between.size),
                                         "outside_indices": np.flatnonzero(~between).tolist(),
                                         "flagged": not bool(np.all(between))},
            "full_channel_vs_production_kk": {**resolution, "flag_above": RESOLUTION_FLAG,
                                              "flagged": resolution["max_abs_over_peak"] > RESOLUTION_FLAG,
                                              "note": "n_chi 4096/k 2048 full-channel run against n_chi 2048/k 1024"},
            "kk_positive_angles": int(np.sum(prod["kk"] > 0))}


def export_inputs(ns: argparse.Namespace) -> dict:
    """Production batches, gate evidence and completed results of all three pipelines."""
    batches = {}
    for pipeline in ("order0", "ff", "fk"):
        stage = PRODUCTION_STAGE[pipeline]
        work = new_dir_arg(getattr(ns, f"{pipeline}_dir"), STAGES[stage].default_dir)
        prepared, plan = load_prepared(work)
        require_stage(plan, stage)
        require_current(plan, prepared)
        gate_path = getattr(ns, f"{pipeline}_gate") or (HERE / STAGES[GATE_STAGE[pipeline]].default_dir
                                                        / "gate_report.json")
        gate = gate_evidence(Path(gate_path), plan)
        if pipeline == "order0":
            results = {t["name"]: completed_order0(t, prepared["plan_sha256"]) for t in plan["targets"]}
        else:
            results = {t["name"]: completed_rows(t, prepared["plan_sha256"]) for t in plan["targets"]}
        for name, (_, completion) in results.items():
            started = read_json(Path(completion["started"]["path"]))
            if (started.get("gate") or {}).get("gate_report") != gate["gate_report"]:
                raise Refusal(f"{name} started under another gate report than {gate['gate_report']['path']}")
        batches[pipeline] = {"work": work, "prepared": prepared, "plan": plan, "gate": gate, "results": results}
    fk = batches["fk"]
    consistency_path = ns.consistency_report or Path(fk["work"]) / "consistency_report.json"
    batches["consistency"] = consistency_evidence(Path(consistency_path), fk["plan"])
    return batches


def fold_meta(batch: dict, name: str, family: str) -> dict:
    plan, (_, completion) = batch["plan"], batch["results"][name]
    target = next(t for t in plan["targets"] if t["name"] == name)
    return {"schema_version": 1, "scope": SCOPE, "family": family, "source_redshift": target["z"],
            "lambda_mpc": target["source"]["lambda_mpc"], "gl": target["gl"], "n_jobs": target["n_jobs"],
            "component_pairs": [[0, 0]], "raw_result": completion["result"], "completed": completion["completed"],
            "plan_sha256": batch["prepared"]["plan_sha256"], "driver": record(DRIVER),
            "export": "raw numeric members copied; object x/y converted to float64 (n,3); main-grid row order"}


def write_npz_member(path: Path, arrays: dict) -> dict:
    """Write one npz (members in the given order) with 'xb' and check that it reloads exactly."""
    with path.open("xb") as handle:
        np.savez(handle, **arrays)
    with zipfile.ZipFile(path) as archive:
        members = {name: sha_bytes(archive.read(name)) for name in archive.namelist()}
    with np.load(path, allow_pickle=False) as check:
        for key, value in arrays.items():
            value = np.asarray(value)
            if check[key].dtype != value.dtype or not np.array_equal(check[key], value):
                raise Refusal(f"Exported {path}:{key} does not reload exactly")
    return {"sha256": digest(path), "member_sha256": members}


def export_sources(batches: dict) -> dict:
    sources = {"v1_source_map": record(V1_SOURCE_MAP), "v1_tree_main_15360": record(V1_TREE_RESULT),
               "v1_ladder": {m: record(V1_LADDER / m) for m in LADDER_MEMBERS},
               "accepted_main_fk": record(ACCEPTED / "products/main_fk.npz"),
               "accepted_main_ff": record(ACCEPTED / "products/main_ff.npz"),
               "accepted_order0": record(O0_FIG), "accepted_multiz": record(ACCEPTED / "products/multiz.npz"),
               "order0_plot_config": record(O0_PLOT_CONFIG), "prepare_plot_inputs": record(PREPARE_PLOT)}
    for pipeline in ("order0", "ff", "fk"):
        batch = batches[pipeline]
        sources[pipeline] = {"prepared": record(Path(batch["work"]) / "prepared.json"), "gate": batch["gate"],
                             "results": {name: completion for name, (_, completion) in batch["results"].items()}}
    sources["consistency"] = batches["consistency"]
    return sources


def layout_references() -> dict:
    refs = {str(p): (record(p) if p.is_file() else "missing") for p in NHFP_REFERENCES}
    return {"files": refs, "note": ("Reading conventions reproduced from these files (load_products, load_channels, "
                                    "window_check, run_ccl taper); recorded, not pinned.")}


def publish_export(out: Path, manifest: dict, writers: dict) -> dict:
    """Write files and manifest into a sibling staging directory, then rename it to ``out`` (RENAME_EXCL)."""
    staging = out.parent / f".{out.name}.{os.getpid()}.{secrets.token_hex(4)}.staging"
    staging.mkdir(exist_ok=False)
    try:
        outputs = {}
        for name, (writer, layout) in writers.items():
            outputs[name] = {"path": str(out / name), **writer(staging / name), "layout": layout}
        manifest = {**manifest, "outputs": outputs}
        write_new_json(staging / "export_manifest.json", manifest)
        try:
            exclusive_rename(staging, out)
        except FileExistsError as exc:
            raise Refusal(f"--out appeared during export: {out}") from exc
    except BaseException:
        print(f"Export failed; partial files remain in {staging}; --out {out} was not created", file=sys.stderr)
        raise
    check_records([{"path": item["path"], "sha256": item["sha256"]} for item in outputs.values()])
    return manifest


def export_command(ns: argparse.Namespace) -> int:
    """Validate gates, consistency and production results, then publish the NHFP package-style copies."""
    require_invocation()
    out = new_dir_arg(ns.out, "export_nhfp_zs3p0")
    if HERE != PLANNED_HOME or not inside(out, HERE) or out.exists() or out.is_symlink() or not out.parent.is_dir():
        raise Refusal(f"--out must be a new directory inside {PLANNED_HOME}: {out}")
    batches = export_inputs(ns)
    lam = EXPECTED_LAMBDA[Z_PRODUCTION]
    fk = fold_export_arrays(batches["fk"]["results"]["fk_main_zs3p0"][0],
                            fold_meta(batches["fk"], "fk_main_zs3p0", "tree-level main FK, l_max=15360, refined table"))
    ff = fold_export_arrays(batches["ff"]["results"]["ff_zs3p0"][0],
                            fold_meta(batches["ff"], "ff_zs3p0", "C23 FF, Gaussian driving fields"))
    o0_results = batches["order0"]["results"]
    full_target = next(t for t in batches["order0"]["plan"]["targets"] if t["name"] == "o0_full_zs3p0")
    self_test = order0_transform_self_test()
    payload, provenance = order0_transform(Path(full_target["out"]), V1_SOURCE_MAP, Z_PRODUCTION)
    provenance["export"] = {"driver": record(DRIVER), "scope": SCOPE,
                            "source_job": "o0_full_zs3p0 (full channels, z=5 [order0] source protocol)"}
    checks = export_checks(fk, ff, payload, lam)
    diagnostics = order0_diagnostics(o0_results["o0_prod_zs3p0"][0], o0_results["o0_probe5_zs3p0"][0],
                                     o0_results["o0_full_zs3p0"][0])
    full_bytes = Path(full_target["out"]).read_bytes()
    o0_arrays = {**payload, "metadata_json": json.dumps(provenance, sort_keys=True)}

    def copy_full(path: Path) -> dict:
        with path.open("xb") as handle:
            handle.write(full_bytes)
        if digest(path) != sha_bytes(full_bytes):
            raise Refusal("order0_comparison copy differs from its source")
        return {"sha256": digest(path), "copied_from": record(Path(full_target["out"]))}

    writers = {"main_fk_zs3p0.npz": (lambda p: write_npz_member(p, fk), "t001_accepted/products/main_fk.npz"),
               "main_ff_zs3p0.npz": (lambda p: write_npz_member(p, ff), "t001_accepted/products/main_ff.npz"),
               "order0_zs3p0.npz": (lambda p: write_npz_member(p, o0_arrays), "figure_products/order0.npz"),
               "order0_comparison_zs3p0.npz": (copy_full, "raw source_folded_reference.py product (ell, ell_taper)")}
    manifest = {"schema": 1, "kind": "NHFP package-style z_s=3.0 export", "scope": SCOPE,
                "status": "CANDIDATE_NO_SCIENTIFIC_ACCEPTANCE", "driver": record(DRIVER),
                "source": {"z": Z_PRODUCTION, "lambda_mpc": lam, "source_map": record(V1_SOURCE_MAP)},
                "sources": export_sources(batches), "checks": checks, "order0_transform_self_test": self_test,
                "order0_diagnostics": diagnostics, "layout_references": layout_references(),
                "nhfp_use": ("Read these files by explicit path and sha256; reproduce/active_products.json is "
                             "unchanged and does not select them."),
                "caveats": ["T-001 acceptance covers finite amplitude and sampled sign at z_s=5 only.",
                            "The 16-shell radial grid of the K table is concentrated near z=5; radial accuracy "
                            "at z_s=3 is not established.",
                            "The z_s=3.0 endpoint falls between the C23 nodes 2250 and 2280 Mpc; C23 source-grid "
                            "accuracy is not granted.",
                            "FF and FK carry the C23 ell <= 5000 tapered window; Order-0 is the windowed SFT "
                            "projection, not the PyCCL all-multipole baseline that NHFP plots.",
                            "Quote every FK value with its l_max (15360 here)."],
                "active_manifest_modified": False}
    manifest = publish_export(out, manifest, writers)
    print(json.dumps({"export_manifest": record(out / "export_manifest.json"),
                      "outputs": {k: v["sha256"] for k, v in manifest["outputs"].items()},
                      "order0_flags": {k: v.get("flagged") for k, v in diagnostics.items() if isinstance(v, dict)}},
                     indent=2))
    return 0


# -------------------------------------------------------------------- CLI
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("plan", "prepare"):
        p = sub.add_parser(command, help=f"{command} one stage")
        p.add_argument("--stage", choices=sorted(STAGES), required=True)
        p.add_argument("--work-dir", type=Path, default=None, help="default: <driver dir>/<stage default>")
        p.add_argument("--n-jobs", type=n_jobs_arg, default="accepted",
                       help="fold sweep workers: 'accepted' (default, the accepted record) or a positive integer")
        if command == "prepare":
            p.add_argument("--allow-stage-test", action="store_true",
                           help="TEST ONLY: allow a driver outside its planned directory; such batches cannot run")
    p = sub.add_parser("run", help="run one prepared fold or Order-0 job under the shared lock")
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--min-available-gib", type=float, default=None)
    p.add_argument("--gate-report", type=Path, default=None,
                   help="PASS gate report (default: <driver dir>/<gate stage default>/gate_report.json)")
    p = sub.add_parser("gate", help="evaluate a gate batch")
    p.add_argument("--work-dir", type=Path, required=True)
    p = sub.add_parser("consistency", help="compare the z_s=3.0 main FK with the v1 tree_main_15360 fold")
    p.add_argument("--work-dir", type=Path, required=True)
    p = sub.add_parser("export", help="write the NHFP package-style copies")
    p.add_argument("--out", type=Path, default=None, help="default: <driver dir>/export_nhfp_zs3p0")
    for pipeline in ("order0", "ff", "fk"):
        p.add_argument(f"--{pipeline}-dir", type=Path, default=None)
        p.add_argument(f"--{pipeline}-gate", type=Path, default=None)
    p.add_argument("--consistency-report", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    ns = parse_args(argv)
    try:
        return {"plan": plan_command, "prepare": prepare_command, "run": run_command, "gate": gate_command,
                "consistency": consistency_command, "export": export_command}[ns.command](ns)
    except Refusal as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    except (OSError, json.JSONDecodeError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        print(f"REFUSED ({type(exc).__name__}): {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
