"""Stubbed tests of driver.py.

No fold starts. run_fk_variant is replaced by a stub that raises unless a test
gives it a writer, and the writer only writes synthetic or accepted rows. The
shared r1_sft061 lock path, the ps check and vm_stat are replaced, and every
file goes into a pytest tmp home that stands in for the driver directory.
run_command is called only where it must refuse before the lock; folds go
through execute_fold with the stub; gate logic goes through evaluate_gates.
Tests run in file order and share the session batches.

Run with the sft-wick interpreter and the driver's pinned environment:
  PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$FW:$CANOES_ROOT/src CANOES_ROOT=$HOME/projects/angular_statistics/canoes \
    $SFT -B -m pytest -p no:cacheprovider --import-mode=importlib --basetemp=<scratch dir> test_driver.py
"""

from __future__ import annotations

import copy
import fcntl
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DRIVER_PATH = Path(__file__).resolve().parent / "driver.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("zs3_driver_under_test", DRIVER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


drv = _load_driver()
Refusal = drv.Refusal
SCALE = {"tree_main_15360": 0.80, "tree_auto_7680": 0.79, "bihalofit_auto_7680": 0.95, "bihalofit_auto_15360": 0.97}


class FoldStub:
    """Stands in for run_fk_variant: a fold happens only through an explicit test writer."""

    def __init__(self):
        self.calls = []
        self.writer = None

    def run(self, config, output, n_gauss, orders):
        self.calls.append(str(output))
        if self.writer is None:
            raise AssertionError("a fold was requested where none may start")
        self.writer(Path(config), Path(output), n_gauss, orders)


def tripwire(*args, **kwargs):
    raise AssertionError("this step must not be reached")


def by_direction(rows):
    return {tuple(float(c) for c in y): float(v) for y, v in zip(rows["y"], rows["value"], strict=True)}


def accepted_values(name):
    if name.startswith("bihalofit"):
        return by_direction(drv.reference_rows(
            drv.ACCEPTED / "cutoff_ladder/bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz"))
    return by_direction(drv.main_subset_rows())


def make_writer(values):
    """Raw fold output in run_FF_single's layout: object x/y vectors plus five numeric members."""
    def write(config, output, n_gauss, orders):
        cfg = json.loads(config.read_text())
        assert os.environ["PYTHONPATH"] == drv.EXPECTED_PYTHONPATH
        assert os.environ["PYTHONNOUSERSITE"] == "1" and os.environ["OMP_NUM_THREADS"] == "1"
        assert n_gauss == cfg["sweep"]["n_gauss"] and orders == drv.ORDERS
        grid, n = cfg["sweep"]["positions_grid"], len(cfg["sweep"]["positions_grid"]["y"])
        xs, ys = np.empty(n, dtype=object), np.empty(n, dtype=object)
        for i, y in enumerate(grid["y"]):
            xs[i], ys[i] = np.asarray(grid["x"][0], float), np.asarray(y, float)
        value = np.array([values[tuple(float(c) for c in y)] for y in grid["y"]])
        with output.open("xb") as handle:
            np.savez(handle, x=xs, y=ys, t_final=np.full(n, cfg["sweep"]["t_final_grid"][0]),
                     a=np.zeros(n, np.int64), b=np.zeros(n, np.int64), order=np.full(n, 2, np.int64), value=value)
    return write


def prepare(argv):
    assert drv.prepare_command(drv.parse_args(["prepare", *argv])) == 0


def target_of(plan, name):
    return next(t for t in plan["targets"] if t["name"] == name)


def run_stub_fold(home, work, name, values, gate=None):
    prepared, plan = drv.load_prepared(work)
    target = target_of(plan, name)
    drv.refuse_existing(target)
    home.stub.writer = make_writer(values)
    try:
        drv.execute_fold(Path(plan["args"]["work_dir"]), prepared, plan, target,
                         drv.fold_environment(Path(target["root"])), {"test": "memory admission stubbed"}, gate)
    finally:
        home.stub.writer = None


@pytest.fixture(scope="session")
def home(tmp_path_factory):
    patch = pytest.MonkeyPatch()
    root = tmp_path_factory.mktemp("home")
    stub = FoldStub()
    patch.setattr(drv, "HERE", root)
    patch.setattr(drv, "PLANNED_HOME", root)
    patch.setattr(drv, "LOCK", root / "test_fold.lock")
    patch.setattr(drv, "variant", lambda: stub)
    patch.setattr(drv, "available_memory_bytes", lambda: 512 * 2**30)
    patch.setattr(drv.shared(), "check_active_folds", lambda process_text=None: None)
    source_map = root / "source_map_zs3p0.json"
    assert drv.source_map_command(3.0, source_map) == 0
    gates, prod = root / "gates_zs5_gl36", root / "folds_zs3p0_gl36"
    prepare(["--source-map", str(drv.REVIEWED_SOURCE_MAP), "--z", "5.0", "--gl", "36", "--n-jobs", "8",
             "--targets", *drv.GATES, "--work-dir", str(gates)])
    prepare(["--source-map", str(source_map), "--z", "3.0", "--gl", "36", "--n-jobs", "8",
             "--targets", *drv.PRODUCTION, "--work-dir", str(prod)])
    yield SimpleNamespace(root=root, gates=gates, prod=prod, source_map=source_map, stub=stub)
    patch.undo()


# ------------------------------------------------------------------ helpers
def test_publish_new_is_complete_or_absent(tmp_path):
    path = tmp_path / "record.json"

    def reject(candidate):
        assert candidate.read_bytes() == b"payload" and candidate != path
        raise Refusal("validation failed")

    with pytest.raises(Refusal, match="validation failed"):
        drv.publish_new(path, b"payload", validate=reject)
    assert list(tmp_path.iterdir()) == []
    drv.publish_new(path, b"payload")
    assert path.read_bytes() == b"payload" and list(tmp_path.iterdir()) == [path]
    with pytest.raises(Refusal, match="replace"):
        drv.publish_new(path, b"other")
    assert path.read_bytes() == b"payload"


def test_exclusive_rename_never_replaces(tmp_path):
    source, existing_file, empty_dir = tmp_path / "source", tmp_path / "file", tmp_path / "empty"
    source.write_text("new")
    existing_file.write_text("old")
    empty_dir.mkdir()
    for target in (existing_file, empty_dir):
        with pytest.raises(FileExistsError):
            drv.exclusive_rename(source, target)
    assert existing_file.read_text() == "old" and list(empty_dir.iterdir()) == []
    drv.exclusive_rename(source, tmp_path / "fresh")
    assert (tmp_path / "fresh").read_text() == "new" and not source.exists()


def test_check_records_refuses_changed_or_missing(tmp_path):
    pinned = tmp_path / "pinned.txt"
    pinned.write_text("a")
    item = drv.record(pinned)
    drv.check_records([item])
    pinned.write_text("b")
    with pytest.raises(Refusal, match="changed"):
        drv.check_records([item])
    pinned.unlink()
    with pytest.raises(Refusal, match="unreadable"):
        drv.check_records([item])


def test_load_module_executes_checked_bytes_without_driver_future_flags(tmp_path, monkeypatch):
    path = tmp_path / "reviewed_probe.py"
    path.write_text("def probe(x: int) -> int:\n    return x\n\n\nHERE = __file__\n")
    monkeypatch.setitem(drv.PINNED, path, "0" * 64)
    with pytest.raises(Refusal, match="changed"):
        drv.load_module("zs3_probe_bad", path)
    monkeypatch.setitem(drv.PINNED, path, drv.digest(path))
    module = drv.load_module("zs3_probe_ok", path)
    try:
        assert module.probe.__annotations__["x"] is int
        here = module.HERE
        assert here == str(path) and sys.modules["zs3_probe_ok"] is module
    finally:
        drv._MODULES.pop("zs3_probe_ok", None)
        sys.modules.pop("zs3_probe_ok", None)


def test_cli_arguments():
    base = ["prepare", "--source-map", "m", "--z", "3", "--gl", "36", "--targets", "tree_main_15360",
            "--work-dir", "w", "--n-jobs"]
    for bad in ("²", "0", "-1", "8.0", "eight", "٨"):
        with pytest.raises(SystemExit):
            drv.parse_args([*base, bad])
    with pytest.raises(SystemExit):
        drv.parse_args([*base, "8", "--tag", "zs5"])
    assert drv.parse_args([*base, "accepted"]).n_jobs == "accepted"
    assert drv.parse_args([*base, "12"]).n_jobs == 12
    assert drv.require_tag("zs3p0") == "zs3p0"
    for bad in ("../zs3p0", "zs3/p0", "", None):
        with pytest.raises(Refusal):
            drv.require_tag(bad)


def test_require_invocation_refuses_wrong_environment(monkeypatch):
    drv.require_invocation()
    monkeypatch.setenv("PYTHONPATH", "/elsewhere")
    with pytest.raises(Refusal, match="PYTHONPATH"):
        drv.require_invocation()
    monkeypatch.undo()
    monkeypatch.setattr(sys, "executable", "/usr/bin/python3")
    with pytest.raises(Refusal, match="interpreter"):
        drv.require_invocation()


def test_memory_admission(monkeypatch):
    monkeypatch.setattr(drv, "available_memory_bytes", lambda: 40 * 2**30)
    admitted = drv.memory_admission(4, None)
    assert admitted["required_bytes"] == drv.MEASURED_PEAK_BYTES[4] + drv.MEMORY_FLOOR_BYTES
    with pytest.raises(Refusal, match="Available memory"):
        drv.memory_admission(8, None)
    with pytest.raises(Refusal, match="No measured peak"):
        drv.memory_admission(6, None)
    assert drv.memory_admission(6, 30.0)["rule"] == "explicit"
    with pytest.raises(Refusal, match="Cannot parse vm_stat"):
        drv.parse_vm_stat("garbage")
    sample = ("Mach Virtual Memory Statistics: (page size of 16384 bytes)\nPages free:   10.\n"
              "Pages active:  5.\nPages inactive:   6.\n")
    assert drv.parse_vm_stat(sample) == 16 * 16384


def test_active_fold_check_errors_refuse(monkeypatch):
    for error in (PermissionError("ps blocked"), RuntimeError("Another fold is active: run_FF_single.py")):
        def raise_error(process_text=None, error=error):
            raise error
        monkeypatch.setattr(drv.shared(), "check_active_folds", raise_error)
        with pytest.raises(Refusal, match="Active-fold check"):
            drv.require_no_active_folds()


def test_refuse_production_location(tmp_path, monkeypatch):
    monkeypatch.setattr(drv, "HERE", tmp_path)
    work = tmp_path / "batch"
    target = {"result": str(work / "t" / "xi_t.npz"),
              "template_config": {"path": str(drv.CAND / "configs/tree_7680.yaml")}}
    drv.refuse_production_location(target, work)
    with pytest.raises(Refusal, match="own target root"):
        drv.refuse_production_location({**target, "result": str(tmp_path / "elsewhere" / "xi.npz")}, work)
    work.mkdir()
    (work / "PRODUCTION").write_text("")
    with pytest.raises(Refusal, match="production"):
        drv.refuse_production_location(target, work)


# ------------------------------------------------------- source map, prepare
def test_source_map_published_and_restricted(home):
    mapping = json.loads(home.source_map.read_text())
    assert mapping["z"] == [3.0] and mapping["lambda_mpc"] == [drv.EXPECTED_LAMBDA[3.0]]
    assert mapping["target_tag"] == "zs3p0" and mapping["producer"]["driver"]["sha256"] == drv.digest(drv.DRIVER)
    assert not list(home.root.glob(".*.partial"))
    with pytest.raises(Refusal, match="new file"):
        drv.source_map_command(3.0, home.source_map)
    with pytest.raises(Refusal, match="only for z_s"):
        drv.source_map_command(2.0, home.root / "source_map_zs2p0.json")


def test_load_endpoint_serves_only_the_two_planned_endpoints(home):
    gate = drv.load_endpoint(drv.REVIEWED_SOURCE_MAP, 5.0)
    assert (gate["tag"], gate["t_max"], gate["lambda_mpc"]) == ("zs5", 2318.0, drv.REVIEWED_Z5_LAMBDA)
    prod = drv.load_endpoint(home.source_map, 3.0)
    assert (prod["tag"], prod["t_max"], prod["lambda_mpc"]) == ("zs3p0", 2256.0, 2255.3607197066854)
    for path, z in ((drv.REVIEWED_SOURCE_MAP, 3.0), (drv.REVIEWED_SOURCE_MAP, 2.5), (home.source_map, 5.0)):
        with pytest.raises(Refusal):
            drv.load_endpoint(path, z)


@pytest.mark.parametrize("index, field, value, message", [
    (0, "target_tag", "zs5", "target_tag"),
    (1, "target_tag", "../zs3p0", "target_tag"),
    (2, "lambda_mpc", [2255.3607197066854 + 1e-9], "re-derived"),
    (3, "producer", {"driver": {"path": "elsewhere/driver.py", "sha256": "0" * 64}}, "written by this driver"),
])
def test_load_endpoint_refuses_tampered_maps(home, index, field, value, message):
    mapping = json.loads(home.source_map.read_text())
    mapping[field] = value
    path = home.root / f"tampered_source_map_{index}.json"
    path.write_text(json.dumps(mapping))
    with pytest.raises(Refusal, match=message):
        drv.load_endpoint(path, 3.0)


def test_prepared_batches(home):
    batches = {"gate": drv.load_prepared(home.gates), "prod": drv.load_prepared(home.prod)}
    for key, names, z in (("gate", drv.GATES, 5.0), ("prod", drv.PRODUCTION, 3.0)):
        prepared, plan = batches[key]
        assert plan["source"]["z"] == z and plan["stage_test"] is False
        assert sorted(t["name"] for t in plan["targets"]) == sorted(names)
        for target in plan["targets"]:
            assert (target["n_jobs"], target["gl"], target["leg"]) == (8, 36, "auto")
            assert prepared["loader_self_check"][target["id"]] == f"GUARD_INSTALLED batch {target['table']}"
            drv.check_records(prepared["generated"][target["id"]])
    assert drv.is_gate_plan(batches["gate"][1]) and not drv.is_production_plan(batches["gate"][1])
    assert drv.is_production_plan(batches["prod"][1]) and not drv.is_gate_plan(batches["prod"][1])
    assert not list(home.root.rglob("*.partial"))


def test_diff_against_accepted_refuses_out_of_allowance_changes(home):
    _, plan = drv.load_prepared(home.prod)
    accepted, source = drv.accepted_records(), plan["source"]
    base, spec = target_of(plan, "bihalofit_auto_7680"), drv.TARGETS["bihalofit_auto_7680"]
    assert drv.diff_against_accepted(base, spec, accepted, source, 36)["allowed_categories"] == [
        "endpoint", "path", "workers"]
    edits = []
    changed = copy.deepcopy(base)
    changed["config"]["sweep"]["vertex_types"] = ["FF"]
    edits.append((changed, "sweep.vertex_types"))
    changed = copy.deepcopy(base)
    changed["loader_text"] += "\n# edited\n"
    edits.append((changed, "loader text"))
    changed = copy.deepcopy(base)
    changed["n_jobs_source"] = "accepted:bihalofit_auto_7680_gl36"
    edits.append((changed, "sweep.n_jobs"))
    common = {p["path"] for p in accepted["bihalofit_auto_7680_gl36"]["contract"]["pins"]}
    changed = copy.deepcopy(base)
    index = next(i for i, p in enumerate(changed["contract"]["pins"]) if p["path"] in common)
    changed["contract"]["pins"][index] = {**changed["contract"]["pins"][index], "sha256": "0" * 64}
    edits.append((changed, "pins changed"))
    for edited, message in edits:
        with pytest.raises(Refusal, match=message):
            drv.diff_against_accepted(edited, spec, accepted, source, 36)
    _, gate_plan = drv.load_prepared(home.gates)
    gate_target = copy.deepcopy(target_of(gate_plan, "bihalofit_auto_15360"))
    gate_target["config"]["sweep"]["t_final_grid"] = [source["lambda_mpc"]]
    with pytest.raises(Refusal, match="sweep.t_final_grid"):
        drv.diff_against_accepted(gate_target, drv.TARGETS["bihalofit_auto_15360"], accepted,
                                  gate_plan["source"], 36)


def test_load_prepared_refusals(home):
    with pytest.raises(Refusal, match="No prepared.json"):
        drv.load_prepared(home.root / "missing_batch")
    copied = home.root / "copied_batch"
    copied.mkdir()
    shutil.copy2(home.prod / "prepared.json", copied / "prepared.json")
    with pytest.raises(Refusal, match="prepared for"):
        drv.load_prepared(copied)
    assert drv.main(["assemble", "--work-dir", str(home.root / "missing_batch"),
                     "--out", str(home.root / "never_written")]) == 2
    assert not (home.root / "never_written").exists()


def test_require_current_refuses_changed_records(home):
    prepared, plan = drv.load_prepared(home.prod)
    drv.require_current(plan, prepared)
    tampered = copy.deepcopy(prepared)
    tampered["generated"][plan["targets"][0]["id"]][0]["sha256"] = "0" * 64
    with pytest.raises(Refusal, match="changed"):
        drv.require_current(plan, tampered)
    stale = copy.deepcopy(plan)
    stale["targets"][0]["n_jobs"] = 3
    with pytest.raises(Refusal, match="changed since preparation"):
        drv.require_current(stale, prepared)


def test_prepare_failure_leaves_a_marker_and_no_batch(home, monkeypatch):
    def failing_self_check(target):
        raise Refusal("injected self-check failure")

    monkeypatch.setattr(drv, "loader_self_check", failing_self_check)
    work = home.root / "failed_prepare"
    argv = ["--source-map", str(drv.REVIEWED_SOURCE_MAP), "--z", "5.0", "--gl", "36", "--n-jobs", "8",
            "--targets", "bihalofit_auto_15360", "--work-dir", str(work)]
    with pytest.raises(Refusal, match="injected"):
        prepare(argv)
    marker = json.loads((work / "prepare_failed.json").read_text())
    assert marker["status"] == "PREPARE_FAILED" and not (work / "prepared.json").exists()
    with pytest.raises(Refusal, match="No prepared.json"):
        drv.load_prepared(work)
    with pytest.raises(Refusal, match="new directory"):
        prepare(argv)


# ----------------------------------------------------------- run refusals
def test_run_refuses_without_gate_report_before_the_lock(home, monkeypatch):
    monkeypatch.setattr(drv, "shared_fold_lock", tripwire)
    ns = drv.parse_args(["run", "--work-dir", str(home.prod), "--target", "tree_auto_7680_zs3p0_gl36"])
    with pytest.raises(Refusal, match="No gate report"):
        drv.run_command(ns)
    _, plan = drv.load_prepared(home.prod)
    assert home.stub.calls == [] and not (Path(target_of(plan, "tree_auto_7680")["root"]) / "started.json").exists()


def test_run_refuses_a_held_lock(home):
    with drv.LOCK.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ns = drv.parse_args(["run", "--work-dir", str(home.gates), "--target", "bihalofit_auto_15360_zs5_gl36"])
        with pytest.raises(Refusal, match="shared fold lock"):
            drv.run_command(ns)
    _, plan = drv.load_prepared(home.gates)
    assert home.stub.calls == []
    assert not (Path(target_of(plan, "bihalofit_auto_15360")["root"]) / "started.json").exists()


# ------------------------------------------------------- gate (stubbed folds)
def test_completed_rows_refuses_an_unfinished_target(home):
    prepared, plan = drv.load_prepared(home.gates)
    with pytest.raises(Refusal, match="has not completed"):
        drv.completed_rows(plan["targets"][0], prepared["plan_sha256"])


def test_gate_folds_stubbed(home):
    caller_env = dict(os.environ)
    for name in drv.GATES:
        run_stub_fold(home, home.gates, name, accepted_values(name))
    assert dict(os.environ) == caller_env
    _, plan = drv.load_prepared(home.gates)
    for target in plan["targets"]:
        root = Path(target["root"])
        completed = json.loads((root / "completed.json").read_text())
        started = json.loads((root / "started.json").read_text())
        assert completed["output_validation"]["rows"] == 11 and completed["scientific_acceptance"] is False
        assert "max_single_process_rss_bytes" in completed and started["gate"] is None
        assert started["environment"]["PYTHONNOUSERSITE"] == "1"
        with pytest.raises(Refusal, match="existing or interrupted"):
            drv.refuse_existing(target)


def test_evaluate_gates_pass(home):
    prepared, plan = drv.load_prepared(home.gates)
    report = drv.evaluate_gates(home.gates, prepared, plan)
    assert report["status"] == "PASS" and report["bitwise_all"] is True
    assert {k: [c["pass"] for c in v] for k, v in report["checks"].items()} == {
        "G1_bihalofit_auto_15360": [True], "G2_tree_main_15360": [True, True]}
    drv.write_new_json(home.gates / "gate_report.json", report)


def test_evaluate_gates_fail_and_untrusted_reports(home, monkeypatch):
    prepared, plan = drv.load_prepared(home.gates)
    original = drv.reference_rows

    def perturbed(path):
        rows = original(path)
        return {**rows, "value": rows["value"] * (1 + 3e-12)}

    monkeypatch.setattr(drv, "reference_rows", perturbed)
    failing = drv.evaluate_gates(home.gates, prepared, plan)
    monkeypatch.undo()
    assert failing["status"] == "FAIL" and failing["checks"]["G1_bihalofit_auto_15360"][0]["pass"] is False
    drv.write_new_json(home.gates / "gate_report_fail.json", failing)
    drv.write_new_json(home.gates / "gate_report_forged.json", {**failing, "status": "PASS"})
    misplaced = home.root / "gate_report_misplaced.json"
    shutil.copy2(home.gates / "gate_report.json", misplaced)
    _, prod_plan = drv.load_prepared(home.prod)
    for path, message in ((home.gates / "gate_report_fail.json", "not PASS"),
                          (home.gates / "gate_report_forged.json", "recomputation"),
                          (misplaced, "does not sit"), (home.gates / "absent.json", "No gate report")):
        with pytest.raises(Refusal, match=message):
            drv.gate_evidence(path, prod_plan)


def test_gate_evidence_binds_workers_and_inputs(home):
    _, prod_plan = drv.load_prepared(home.prod)
    report = home.gates / "gate_report.json"
    evidence = drv.gate_evidence(report, prod_plan)
    assert evidence["status"] == "PASS" and evidence["bitwise_all"] is True
    mismatch = copy.deepcopy(prod_plan)
    target_of(mismatch, "bihalofit_auto_15360")["n_jobs"] = 12
    with pytest.raises(Refusal, match="differs from gate target"):
        drv.gate_evidence(report, mismatch)
    control = copy.deepcopy(prod_plan)
    for target in control["targets"]:
        target["gl"], target["n_jobs"] = 24, 4
    assert drv.gate_evidence(report, control)["status"] == "PASS"
    other = copy.deepcopy(prod_plan)
    other["driver"] = {"path": "elsewhere/driver.py", "sha256": "0" * 64}
    with pytest.raises(Refusal, match="another driver"):
        drv.gate_evidence(report, other)


def test_run_refuses_failed_gate_and_worker_mismatch_before_the_lock(home, monkeypatch):
    work = home.root / "folds_zs3p0_gl36_w4"
    prepare(["--source-map", str(home.source_map), "--z", "3.0", "--gl", "36", "--n-jobs", "4",
             "--targets", "bihalofit_auto_15360", "--work-dir", str(work)])
    monkeypatch.setattr(drv, "shared_fold_lock", tripwire)
    failed = drv.parse_args(["run", "--work-dir", str(home.prod), "--target", "tree_auto_7680_zs3p0_gl36",
                             "--gate-report", str(home.gates / "gate_report_fail.json")])
    with pytest.raises(Refusal, match="not PASS"):
        drv.run_command(failed)
    mismatch = drv.parse_args(["run", "--work-dir", str(work), "--target", "bihalofit_auto_15360_zs3p0_gl36"])
    with pytest.raises(Refusal, match="differs from gate target"):
        drv.run_command(mismatch)
    assert home.stub.calls == [str(home.gates / t / f"xi_{t}.npz") for t in
                               ("bihalofit_auto_15360_zs5_gl36", "tree_main_15360_zs5_gl36")]


# -------------------------------------------------- production and assembly
def test_production_folds_stubbed(home, monkeypatch):
    prepared, plan = drv.load_prepared(home.prod)
    gate = drv.gate_evidence(home.gates / "gate_report.json", plan)
    for name in drv.PRODUCTION:
        values = {k: v * SCALE[name] for k, v in accepted_values(name).items()}
        run_stub_fold(home, home.prod, name, values, gate)
        started = json.loads((Path(target_of(plan, name)["root"]) / "started.json").read_text())
        assert started["gate"] == gate
    monkeypatch.setattr(drv, "shared_fold_lock", tripwire)
    again = drv.parse_args(["run", "--work-dir", str(home.prod), "--target", "tree_auto_7680_zs3p0_gl36"])
    with pytest.raises(Refusal, match="existing or interrupted"):
        drv.run_command(again)


def test_assemble_refusals_leave_no_out(home, monkeypatch):
    out = home.root / "cutoff_ladder_zs3p0"
    with pytest.raises(Refusal, match="only the z_s=3"):
        drv.assemble_command(drv.parse_args(["assemble", "--work-dir", str(home.gates), "--out", str(out)]))
    ns = drv.parse_args(["assemble", "--work-dir", str(home.prod), "--out", str(out)])

    def broken_plausibility(*args):
        raise Refusal("injected plausibility failure")

    monkeypatch.setattr(drv, "plausibility", broken_plausibility)
    with pytest.raises(Refusal, match="injected"):
        drv.assemble_command(ns)
    assert not out.exists() and not list(home.root.glob(".cutoff_ladder_zs3p0.*"))
    monkeypatch.undo()
    original, calls = drv.write_member, []

    def flaky_member(path, rows, meta):
        calls.append(path)
        if len(calls) == 3:
            raise OSError("injected disk failure")
        return original(path, rows, meta)

    monkeypatch.setattr(drv, "write_member", flaky_member)
    with pytest.raises(OSError, match="injected"):
        drv.assemble_command(ns)
    assert not out.exists() and len(list(home.root.glob(".cutoff_ladder_zs3p0.*.staging"))) == 1


def test_assemble_writes_members_and_manifest(home):
    out = home.root / "cutoff_ladder_zs3p0"
    ns = drv.parse_args(["assemble", "--work-dir", str(home.prod), "--out", str(out)])
    assert drv.assemble_command(ns) == 0
    manifest = json.loads((out / "assembly_manifest.json").read_text())
    expected = sorted(drv.TARGETS[name].member for name in drv.PRODUCTION)
    assert sorted(manifest["outputs"]) == expected
    assert sorted(str(p.relative_to(out)) for p in out.rglob("*.npz")) == expected
    _, plan = drv.load_prepared(home.prod)
    for name in drv.PRODUCTION:
        target, member = target_of(plan, name), drv.TARGETS[name].member
        item = manifest["outputs"][member]
        assert item["path"] == str(out / member) and item["sha256"] == drv.digest(out / member)
        raw = drv.read_rows(Path(target["result"]))
        with np.load(out / member, allow_pickle=False) as data:
            meta = json.loads(str(data["metadata_json"]))
            assert (meta["leg_order"], meta["ell_max"], meta["source_redshift"]) == ("auto", target["cutoff"], 3.0)
            assert np.array_equal(data["value"], raw["value"]) and np.array_equal(data["y"], raw["y"])
            assert data["x"].dtype == np.float64 and data["x"].shape == (11, 3)
    assert manifest["gate"]["status"] == "PASS" and manifest["source"]["z"] == 3.0
    assert manifest["outputs"]["tree_auto/table_tree_cut15360_r4_xi.npz"]["target_id"] == "tree_main_15360_zs3p0_gl36"
    plaus = manifest["plausibility"]
    assert len(plaus["rows"]) == 11 and plaus["acceptance_threshold"] is None
    with pytest.raises(Refusal, match="new directory"):
        drv.assemble_command(ns)
