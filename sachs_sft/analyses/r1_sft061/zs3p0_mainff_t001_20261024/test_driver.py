"""Stubbed tests of driver.py (version 2).

No fold and no Order-0 job starts. run_fk_variant is replaced by a stub that raises
unless a test gives it a writer, and launch_order0 by a stub that writes copies of
accepted products. The shared r1_sft061 lock path, the ps table and vm_stat are
replaced, and every file goes into a pytest tmp home that stands in for the driver
directory. Tests run in file order and share the session batches.

Run with the sft-wick interpreter and the driver's pinned environment:
  PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$FW:$CANOES_ROOT/src CANOES_ROOT=$HOME/projects/angular_statistics/canoes \
    $SFT -B -m pytest -p no:cacheprovider --import-mode=importlib --basetemp=<scratch dir> test_driver.py
"""

from __future__ import annotations

import copy
import fcntl
import importlib.util
import json
import shutil
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DRIVER_PATH = Path(__file__).resolve().parent / "driver.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("zs3v2_driver_under_test", DRIVER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


drv = _load_driver()
Refusal = drv.Refusal
O0_PROBE_NPZ = drv.R1A / "source_folded_multiz_z3p2_probe.npz"


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


class Order0Stub:
    """Stands in for launch_order0: a job happens only through an explicit test writer."""

    def __init__(self):
        self.calls = []
        self.writer = None

    def __call__(self, target, env, log_handle):
        self.calls.append(target["id"])
        if self.writer is None:
            raise AssertionError("an Order-0 job was requested where none may start")
        return self.writer(target, env, log_handle)


def tripwire(*args, **kwargs):
    raise AssertionError("this step must not be reached")


def key(y):
    return tuple(float(c) for c in y)


def by_direction(ys, values):
    return {key(y): float(v) for y, v in zip(ys, values, strict=True)}


def main_values(product):
    rows = drv.main_pair_rows(drv.ACCEPTED / "products" / product)
    return by_direction(rows["y"], rows["value"])


def multiz_values():
    rows = drv.multiz_rows(tuple(range(drv.N_MAIN)))
    return by_direction(rows["y"], rows["value"])


def v1_values():
    rows = drv.v1_tree_rows()
    return by_direction(rows["y"], rows["value"])


def fold_writer(values):
    """Raw fold output in run_FF_single's layout: object x/y vectors plus five numeric members."""
    def write(config, output, n_gauss, orders):
        cfg = json.loads(config.read_text())
        assert drv.os.environ["PYTHONPATH"] == drv.EXPECTED_PYTHONPATH
        assert drv.os.environ["PYTHONNOUSERSITE"] == "1" and drv.os.environ["OMP_NUM_THREADS"] == "1"
        assert n_gauss == cfg["sweep"]["n_gauss"] and orders == drv.ORDERS
        grid, n = cfg["sweep"]["positions_grid"], len(cfg["sweep"]["positions_grid"]["y"])
        xs, ys = np.empty(n, dtype=object), np.empty(n, dtype=object)
        for i, y in enumerate(grid["y"]):
            xs[i], ys[i] = np.asarray(grid["x"][0], float), np.asarray(y, float)
        value = np.array([values[key(y)] for y in grid["y"]])
        with output.open("xb") as handle:
            np.savez(handle, x=xs, y=ys, t_final=np.full(n, cfg["sweep"]["t_final_grid"][0]),
                     a=np.zeros(n, np.int64), b=np.zeros(n, np.int64), order=np.full(n, 2, np.int64), value=value)
    return write


def order0_writer(source_npz, scale=None):
    """source_folded_reference.py outputs copied from an accepted product, with the planned settings."""
    def write(target, env, log_handle):
        assert env["OPENBLAS_NUM_THREADS"] == "1" and env["CANOES_SUPPRESS_METAL_WARNING"] == "1"
        assert not any(k.startswith(("PYTHONPATH", "SFT_WICK_")) for k in env)
        assert target["argv"][:2] == [str(drv.PYCCL), str(drv.O0_SCRIPT)]
        out = Path(target["out"])
        with np.load(source_npz, allow_pickle=False) as data:
            arrays = {k: np.asarray(data[k]) for k in data.files if k != "metadata_json"}
            meta = json.loads(str(data["metadata_json"]))
        if scale is not None:
            arrays["kk"] = arrays["kk"] * scale
        meta.update(settings=target["expected_settings"], source_lambda_mpc=target["source"]["lambda_mpc"],
                    source_redshift=target["z"], source_map=target["expected_settings"]["source_map"],
                    script_sha256=drv.PINNED[drv.O0_SCRIPT], warnings=[])
        np.savez(out, **arrays, metadata_json=json.dumps(meta, sort_keys=True))
        out.with_suffix(".json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
        out.with_suffix(".source.py").write_bytes(drv.O0_SCRIPT.read_bytes())
        log_handle.write(b"stub Order-0 job\n")
        return 0
    return write


def prepare(stage, extra=()):
    assert drv.prepare_command(drv.parse_args(["prepare", "--stage", stage, *extra])) == 0


def run(home, work, target_id, writer=None, order0=None, gate=None):
    argv = ["run", "--work-dir", str(work), "--target", target_id]
    if gate is not None:
        argv += ["--gate-report", str(gate)]
    home.stub.writer, home.order0.writer = writer, order0
    try:
        return drv.run_command(drv.parse_args(argv))
    finally:
        home.stub.writer = home.order0.writer = None


def gate(work):
    return drv.gate_command(drv.parse_args(["gate", "--work-dir", str(work)]))


def dirs(root):
    return {name: root / stage.default_dir for name, stage in drv.STAGES.items()}


@pytest.fixture(scope="session")
def home(tmp_path_factory):
    patch = pytest.MonkeyPatch()
    root = tmp_path_factory.mktemp("home")
    stub, order0 = FoldStub(), Order0Stub()
    patch.setattr(drv, "HERE", root)
    patch.setattr(drv, "PLANNED_HOME", root)
    patch.setattr(drv, "LOCK", root / "test_fold.lock")
    patch.setattr(drv, "variant", lambda: stub)
    patch.setattr(drv, "launch_order0", order0)
    patch.setattr(drv, "available_memory_bytes", lambda: 512 * 2**30)
    patch.setattr(drv, "process_table", lambda: "")
    for stage in drv.STAGES:
        prepare(stage)
    yield SimpleNamespace(root=root, dirs=dirs(root), stub=stub, order0=order0)
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
    assert path.read_bytes() == b"payload"
    with pytest.raises(Refusal, match="replace"):
        drv.publish_new(path, b"other")


def test_exclusive_rename_never_replaces(tmp_path):
    source, existing, empty = tmp_path / "source", tmp_path / "file", tmp_path / "empty"
    source.write_text("new")
    existing.write_text("old")
    empty.mkdir()
    for target in (existing, empty):
        with pytest.raises(FileExistsError):
            drv.exclusive_rename(source, target)
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


def test_cli_arguments():
    for bad in ("0", "-1", "8.0", "eight", "٨"):
        with pytest.raises(SystemExit):
            drv.parse_args(["prepare", "--stage", "fk_production", "--n-jobs", bad])
    with pytest.raises(SystemExit):
        drv.parse_args(["prepare", "--stage", "fk_ladder"])
    assert drv.parse_args(["prepare", "--stage", "fk_production"]).n_jobs == "accepted"
    assert drv.parse_args(["plan", "--stage", "ff_canary", "--n-jobs", "8"]).n_jobs == 8
    with pytest.raises(Refusal, match="single processes"):
        drv.stage_args(drv.parse_args(["plan", "--stage", "order0_parity", "--n-jobs", "4"]))


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
    assert drv.memory_admission(("fk", 4), None)["required_bytes"] == 27544076288 + drv.MEMORY_FLOOR_BYTES
    assert drv.memory_admission(("ff", 4), None)["required_bytes"] == 28 * 2**30
    assert drv.memory_admission(("order0", 1), None)["rule"] == "recorded peak + 12 GiB"
    with pytest.raises(Refusal, match="Available memory"):
        drv.memory_admission(("fk", 8), None)
    with pytest.raises(Refusal, match="No recorded peak"):
        drv.memory_admission(("ff", 8), None)
    assert drv.memory_admission(("ff", 8), 30.0)["rule"] == "explicit"
    sample = ("Mach Virtual Memory Statistics: (page size of 16384 bytes)\nPages free:   10.\n"
              "Pages active:  5.\nPages inactive:   6.\n")
    assert drv.parse_vm_stat(sample) == 16 * 16384


def test_active_fold_checks_refuse(monkeypatch):
    for line in ("1 python run_FF_single.py a b", "2 python source_folded_reference.py --z-source 3",
                 "3 python .../zs3p0_ladder_t001_20261023/driver.py run"):
        monkeypatch.setattr(drv, "process_table", lambda line=line: line)
        with pytest.raises(Refusal, match="Active-fold check"):
            drv.require_no_active_folds()

    def blocked():
        raise PermissionError("ps blocked")

    monkeypatch.setattr(drv, "process_table", blocked)
    with pytest.raises(Refusal, match="Active-fold check"):
        drv.require_no_active_folds()


def test_t_max_rule():
    lam = drv.EXPECTED_LAMBDA
    assert drv.t_max_for(lam[3.0], 1.0) == 2256.0 and drv.t_max_for(lam[3.0], 2.0) == 2256.0
    assert drv.t_max_for(lam[3.2], 1.0) == 2267.0 and drv.t_max_for(lam[3.2], 2.0) == 2268.0
    assert drv.t_max_for(lam[5.0], 1.0) == 2318.0 and drv.t_max_for(lam[5.0], 2.0) == 2318.0


# ------------------------------------------------------- endpoints and prepare
def test_load_endpoint_serves_only_three_endpoints(tmp_path):
    served = {(drv.REVIEWED_SOURCE_MAP, 3.2): "zs3p2", (drv.REVIEWED_SOURCE_MAP, 5.0): "zs5",
              (drv.V1_SOURCE_MAP, 3.0): "zs3p0"}
    for (path, z), tag in served.items():
        endpoint = drv.load_endpoint(path, z)
        assert endpoint["tag"] == tag and endpoint["lambda_mpc"] == drv.EXPECTED_LAMBDA[z]
    for path, z in ((drv.REVIEWED_SOURCE_MAP, 3.0), (drv.REVIEWED_SOURCE_MAP, 2.5), (drv.V1_SOURCE_MAP, 5.0),
                    (drv.V1_SOURCE_MAP, 3.2)):
        with pytest.raises(Refusal):
            drv.load_endpoint(path, z)
    copied = tmp_path / "source_map_zs3p0.json"
    shutil.copy2(drv.V1_SOURCE_MAP, copied)
    with pytest.raises(Refusal, match="neither"):
        drv.load_endpoint(copied, 3.0)


def test_prepared_batches(home):
    expected = {"order0_parity": 2, "order0_production": 3, "ff_canary": 1, "ff_production": 1,
                "fk_regression": 1, "fk_production": 1}
    for stage, count in expected.items():
        prepared, plan = drv.load_prepared(home.dirs[stage])
        assert plan["stage"]["name"] == stage and plan["stage_test"] is False
        assert len(plan["targets"]) == count
        for target in plan["targets"]:
            drv.check_records(prepared["generated"][target["id"]])
            check = prepared["self_check"][target["id"]]
            if target["pipeline"] == "order0":
                assert check.startswith("ORDER0_RUNTIME_OK") and target["argv"][0] == str(drv.PYCCL)
            else:
                assert target["n_jobs"] == 4 and target["n_jobs_source"].startswith("accepted:")
                assert check == (f"GUARD_INSTALLED batch {drv.REFINED}" if target["pipeline"] == "fk"
                                 else f"C_LOADED {drv.PINNED[drv.C23_TABLE]}")
    _, fk = drv.load_prepared(home.dirs["fk_production"])
    target = fk["targets"][0]
    assert (target["gl"], target["expected_rows"], target["config"]["propagators"]["t_max"]) == (36, 40, 2256.0)
    _, reg = drv.load_prepared(home.dirs["fk_regression"])
    assert reg["targets"][0]["y_indices"] == list(drv.REGRESSION_INDICES)
    _, ff = drv.load_prepared(home.dirs["ff_canary"])
    assert ff["targets"][0]["config"]["sweep"]["positions_grid"]["y"] == [
        drv.main_grid()[1][i] for i in drv.CANARY_INDICES]
    assert not list(home.root.rglob("*.partial"))


def test_order0_argv_reproduces_recorded_settings(home):
    _, plan = drv.load_prepared(home.dirs["order0_parity"])
    parity = plan["targets"][0]
    recorded = [r for r in drv.read_json(drv.O0_PLAN)["production_runs"] if r["z"] == 3.2][0]["argv"]
    assert parity["argv"][:-1] == recorded[:-1] and parity["argv"][-1] == parity["out"]
    assert [r["field"] for r in parity["accepted_diff"]["settings"]] == ["out"]
    tokens = drv.order0_tokens(drv.ORDER0["o0_prod_zs3p0"])
    bad = dict(parity["expected_settings"], n_chi=4096)
    with pytest.raises(Refusal, match="does not reproduce"):
        drv.check_tokens(tokens, bad)


def test_fold_diff_refuses_out_of_allowance_changes(home):
    accepted = drv.accepted_records()
    _, fk = drv.load_prepared(home.dirs["fk_regression"])
    base = fk["targets"][0]
    report = drv.fold_diff(base, accepted, "multiz_production_gl36")
    assert report["loader"]["equal"] and "guard_scope" in report["allowed_categories"]
    edits = []
    changed = copy.deepcopy(base)
    changed["config"]["sweep"]["vertex_types"] = ["FF"]
    edits.append((changed, "sweep.vertex_types"))
    changed = copy.deepcopy(base)
    changed["loader_text"] += "\n# edited\n"
    edits.append((changed, "loader text"))
    changed = copy.deepcopy(base)
    index = next(i for i, p in enumerate(changed["contract"]["pins"]) if p["path"].endswith("aligned_response.py"))
    changed["contract"]["pins"][index] = {**changed["contract"]["pins"][index], "sha256": "0" * 64}
    edits.append((changed, "pins changed"))
    for edited, message in edits:
        for accepted_id in ("main_fk_gl36", "multiz_production_gl36"):
            with pytest.raises(Refusal, match=message):
                drv.fold_diff(edited, accepted, accepted_id)
    _, ff = drv.load_prepared(home.dirs["ff_canary"])
    canary = copy.deepcopy(ff["targets"][0])
    canary["config"]["sweep"]["t_final_grid"] = [drv.EXPECTED_LAMBDA[3.0]]
    with pytest.raises(Refusal, match="sweep.t_final_grid"):
        drv.fold_diff(canary, accepted, "ff_main_release")
    canary = copy.deepcopy(ff["targets"][0])
    canary["config"]["sweep"]["n_jobs"] = 8
    with pytest.raises(Refusal, match="sweep.n_jobs"):
        drv.fold_diff(canary, accepted, "ff_main_release")


def test_load_prepared_and_require_current_refusals(home):
    with pytest.raises(Refusal, match="No prepared.json"):
        drv.load_prepared(home.root / "missing_batch")
    copied = home.root / "copied_batch"
    copied.mkdir()
    shutil.copy2(home.dirs["fk_production"] / "prepared.json", copied / "prepared.json")
    with pytest.raises(Refusal, match="prepared for"):
        drv.load_prepared(copied)
    prepared, plan = drv.load_prepared(home.dirs["fk_production"])
    tampered = copy.deepcopy(prepared)
    tampered["generated"][plan["targets"][0]["id"]][0]["sha256"] = "0" * 64
    with pytest.raises(Refusal, match="changed"):
        drv.require_current(plan, tampered)
    stale = copy.deepcopy(plan)
    stale["targets"][0]["n_jobs"] = 3
    with pytest.raises(Refusal, match="changed since preparation"):
        drv.require_current(stale, prepared)
    assert drv.main(["gate", "--work-dir", str(home.root / "missing_batch")]) == 2


def test_prepare_refusals_and_failure_marker(home, monkeypatch):
    with pytest.raises(Refusal, match="new directory"):
        prepare("fk_production")
    monkeypatch.setattr(drv, "loader_self_check", lambda target: (_ for _ in ()).throw(Refusal("injected")))
    work = home.root / "failed_prepare"
    with pytest.raises(Refusal, match="injected"):
        prepare("fk_regression", ["--work-dir", str(work)])
    assert drv.read_json(work / "prepare_failed.json")["status"] == "PREPARE_FAILED"
    assert not (work / "prepared.json").exists()
    monkeypatch.setattr(drv, "PLANNED_HOME", home.root / "elsewhere")
    with pytest.raises(Refusal, match="must sit at"):
        prepare("fk_regression", ["--work-dir", str(home.root / "never")])


def test_stage_test_batches_cannot_run(home, monkeypatch):
    work = home.root / "stage_test_batch"
    prepare("ff_canary", ["--work-dir", str(work), "--allow-stage-test"])
    with pytest.raises(Refusal, match="Stage-test"):
        drv.load_prepared(work)
    monkeypatch.setattr(drv, "PLANNED_HOME", home.root / "elsewhere")
    with pytest.raises(Refusal, match="Stage-test"):
        drv.load_prepared(home.dirs["ff_canary"])


# ----------------------------------------------------------- run refusals
def test_production_runs_refuse_without_gate_before_the_lock(home, monkeypatch):
    monkeypatch.setattr(drv, "shared_fold_lock", tripwire)
    for stage in ("order0_production", "ff_production", "fk_production"):
        _, plan = drv.load_prepared(home.dirs[stage])
        with pytest.raises(Refusal, match="No gate report"):
            run(home, home.dirs[stage], plan["targets"][0]["id"])
    assert home.stub.calls == [] and home.order0.calls == []


def test_run_refuses_a_held_lock(home):
    with drv.LOCK.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(Refusal, match="shared fold lock"):
            run(home, home.dirs["fk_regression"], "fk_regress_zs3p2_gl36")
    assert home.stub.calls == []


def test_completed_rows_refuses_an_unfinished_target(home):
    prepared, plan = drv.load_prepared(home.dirs["fk_regression"])
    with pytest.raises(Refusal, match="has not completed"):
        drv.completed_rows(plan["targets"][0], prepared["plan_sha256"])


# ------------------------------------------------------- gates (stubbed jobs)
def test_gate_jobs_stubbed(home):
    assert run(home, home.dirs["fk_regression"], "fk_regress_zs3p2_gl36", writer=fold_writer(multiz_values())) == 0
    assert run(home, home.dirs["ff_canary"], "ff_canary_zs5_gl24", writer=fold_writer(main_values("main_ff.npz"))) == 0
    assert run(home, home.dirs["order0_parity"], "o0_parity_zs3p2", order0=order0_writer(drv.O0_Z3P2)) == 0
    assert run(home, home.dirs["order0_parity"], "o0_parity_full_zs5", order0=order0_writer(drv.O0_Z5)) == 0
    _, plan = drv.load_prepared(home.dirs["fk_regression"])
    root = Path(plan["targets"][0]["root"])
    completed, started = drv.read_json(root / "completed.json"), drv.read_json(root / "started.json")
    assert completed["output_validation"]["rows"] == 6 and started["gate"] is None
    assert started["environment"]["PYTHONNOUSERSITE"] == "1"
    with pytest.raises(Refusal, match="existing or interrupted"):
        run(home, home.dirs["fk_regression"], "fk_regress_zs3p2_gl36")


def test_gate_fail_reports_are_recorded_not_replaced(home, monkeypatch):
    reports = {}
    for stage, reference, label in (("fk_regression", "multiz_rows", "K1_fk_regression_zs3p2"),
                                    ("ff_canary", "main_pair_rows", "F1_ff_canary_zs5")):
        prepared, plan = drv.load_prepared(home.dirs[stage])
        original = getattr(drv, reference)
        monkeypatch.setattr(drv, reference, lambda *a, original=original: {
            **original(*a), "value": original(*a)["value"] * (1 + 3e-12)})
        reports[stage] = drv.evaluate_gate(home.dirs[stage], prepared, plan)
        monkeypatch.undo()
        assert reports[stage]["status"] == "FAIL" and reports[stage]["checks"][label][0]["pass"] is False
    drv.write_new_json(home.dirs["fk_regression"] / "gate_report_fail.json", reports["fk_regression"])
    prepared, plan = drv.load_prepared(home.dirs["order0_parity"])
    original_load = drv.load_arrays

    def perturbed(path, keys):
        arrays = original_load(path, keys)
        if Path(path) == drv.O0_Z3P2 and "kk" in arrays:
            arrays["kk"] = arrays["kk"] * (1 + 1e-4)
        return arrays

    monkeypatch.setattr(drv, "load_arrays", perturbed)
    failing = drv.evaluate_gate(home.dirs["order0_parity"], prepared, plan)
    monkeypatch.undo()
    assert failing["status"] == "FAIL" and failing["checks"]["o0_parity_zs3p2"][0]["pass"] is False
    assert failing["checks"]["o0_parity_full_zs5"][0]["pass"] is True


def test_gates_pass_bitwise(home):
    for stage in ("fk_regression", "ff_canary", "order0_parity"):
        assert gate(home.dirs[stage]) == 0
        report = drv.read_json(home.dirs[stage] / "gate_report.json")
        assert report["status"] == "PASS" and report["bitwise_all"] is True
        with pytest.raises(Refusal, match="never replaced"):
            drv.gate_command(drv.parse_args(["gate", "--work-dir", str(home.dirs[stage])]))
    report = drv.read_json(home.dirs["order0_parity"] / "gate_report.json")
    kk = report["checks"]["o0_parity_zs3p2"][0]["arrays"]["kk"]
    assert kk["bitwise"] is True and kk["max_abs_over_peak"] == 0.0
    with pytest.raises(Refusal, match="Gates take only"):
        prepared, plan = drv.load_prepared(home.dirs["fk_production"])
        drv.evaluate_gate(home.dirs["fk_production"], prepared, plan)


def test_gate_evidence_binding_and_forgeries(home):
    _, prod = drv.load_prepared(home.dirs["fk_production"])
    report = home.dirs["fk_regression"] / "gate_report.json"
    assert drv.gate_evidence(report, prod)["status"] == "PASS"
    mismatch = copy.deepcopy(prod)
    mismatch["targets"][0]["n_jobs"] = 8
    with pytest.raises(Refusal, match="differs from gate target"):
        drv.gate_evidence(report, mismatch)
    other = copy.deepcopy(prod)
    other["driver"] = {"path": "elsewhere/driver.py", "sha256": "0" * 64}
    with pytest.raises(Refusal, match="another driver"):
        drv.gate_evidence(report, other)
    forged = home.dirs["fk_regression"] / "gate_report_forged.json"
    drv.write_new_json(forged, {**drv.read_json(report), "bitwise_all": False})
    misplaced = home.root / "gate_report_misplaced.json"
    shutil.copy2(report, misplaced)
    _, ff_prod = drv.load_prepared(home.dirs["ff_production"])
    for path, plan, message in ((home.dirs["fk_regression"] / "gate_report_fail.json", prod, "not PASS"),
                                (forged, prod, "recomputation"), (misplaced, prod, "does not sit"),
                                (report, ff_prod, "not ff_canary")):
        with pytest.raises(Refusal, match=message):
            drv.gate_evidence(path, plan)


# -------------------------------------------------- production and export
def fk_production_values():
    values = {k: 0.5 * v for k, v in main_values("main_fk.npz").items()}
    values.update(v1_values())
    return values


def test_production_jobs_stubbed(home, monkeypatch):
    with pytest.raises(Refusal, match="not PASS"):
        run(home, home.dirs["fk_production"], "fk_main_zs3p0_gl36",
            gate=home.dirs["fk_regression"] / "gate_report_fail.json")
    assert run(home, home.dirs["fk_production"], "fk_main_zs3p0_gl36", writer=fold_writer(fk_production_values())) == 0
    ff = {k: 0.4 * v for k, v in main_values("main_ff.npz").items()}
    assert run(home, home.dirs["ff_production"], "ff_zs3p0_gl24", writer=fold_writer(ff)) == 0
    for name, source in (("o0_prod_zs3p0", drv.O0_Z3P2), ("o0_probe5_zs3p0", O0_PROBE_NPZ),
                         ("o0_full_zs3p0", drv.O0_Z5)):
        assert run(home, home.dirs["order0_production"], name, order0=order0_writer(source)) == 0
    _, plan = drv.load_prepared(home.dirs["fk_production"])
    started = drv.read_json(Path(plan["targets"][0]["root"]) / "started.json")
    assert started["gate"]["status"] == "PASS"
    assert started["gate"]["gate_report"]["path"] == str(home.dirs["fk_regression"] / "gate_report.json")
    monkeypatch.setattr(drv, "shared_fold_lock", tripwire)
    with pytest.raises(Refusal, match="existing or interrupted"):
        run(home, home.dirs["ff_production"], "ff_zs3p0_gl24")


def test_order0_output_validation_refuses_wrong_settings(home, tmp_path):
    _, plan = drv.load_prepared(home.dirs["order0_production"])
    target = copy.deepcopy(next(t for t in plan["targets"] if t["name"] == "o0_prod_zs3p0"))
    target["expected_settings"] = {**target["expected_settings"], "n_chi": 4096}
    with pytest.raises(Refusal, match="recorded settings differ"):
        drv.validate_order0_output(target)


def test_consistency_report(home, monkeypatch):
    prepared, plan = drv.load_prepared(home.dirs["fk_production"])
    original = drv.v1_tree_rows
    monkeypatch.setattr(drv, "v1_tree_rows", lambda: {**original(), "value": original()["value"] * (1 + 1e-11)})
    assert drv.evaluate_consistency(home.dirs["fk_production"], prepared, plan)["status"] == "FAIL"
    monkeypatch.undo()
    ns = drv.parse_args(["consistency", "--work-dir", str(home.dirs["fk_production"])])
    assert drv.consistency_command(ns) == 0
    report = drv.read_json(home.dirs["fk_production"] / "consistency_report.json")
    assert report["status"] == "PASS" and report["bitwise"] is True
    with pytest.raises(Refusal, match="never replaced"):
        drv.consistency_command(ns)
    with pytest.raises(Refusal, match="fk_production"):
        drv.evaluate_consistency(home.dirs["ff_production"], *drv.load_prepared(home.dirs["ff_production"]))


def test_export_refusal_leaves_no_out(home, monkeypatch):
    out = home.root / "export_nhfp_zs3p0"
    monkeypatch.setattr(drv, "order0_transform_self_test", lambda: (_ for _ in ()).throw(Refusal("injected")))
    with pytest.raises(Refusal, match="injected"):
        drv.export_command(drv.parse_args(["export"]))
    assert not out.exists()
    monkeypatch.undo()
    original, calls = drv.write_npz_member, []

    def flaky(path, arrays):
        calls.append(path)
        if len(calls) == 2:
            raise OSError("injected disk failure")
        return original(path, arrays)

    monkeypatch.setattr(drv, "write_npz_member", flaky)
    with pytest.raises(OSError, match="injected"):
        drv.export_command(drv.parse_args(["export"]))
    assert not out.exists() and len(list(home.root.glob(".export_nhfp_zs3p0.*.staging"))) == 1


def nhfp_load_channels(path, order):
    """nhfp_fig_evidence.load_channels, verbatim apart from the error type."""
    with np.load(path, allow_pickle=False) as d:
        a, b, o, v, xs, ys = (np.asarray(d[k]) for k in ("a", "b", "order", "value", "x", "y"))
    grid, out = None, {}
    for pair in sorted({(int(i), int(j)) for i, j in zip(a, b)}):  # noqa: B905 (verbatim NHFP code)
        idx = np.flatnonzero((a == pair[0]) & (b == pair[1]) & (o == order))
        if idx.size == 0:
            continue
        sep = np.array([drv.separation_arcmin(xs[i], ys[i]) for i in idx])
        s = np.argsort(sep)
        if grid is None:
            grid = sep[s]
        elif not np.array_equal(grid, sep[s]):
            raise ValueError(f"{path}: separation grid differs for pair {pair}")
        out[pair] = v[idx][s]
    return grid, out


def test_export_writes_nhfp_layout(home):
    out = home.root / "export_nhfp_zs3p0"
    assert drv.export_command(drv.parse_args(["export"])) == 0
    manifest = drv.read_json(out / "export_manifest.json")
    names = {"main_fk_zs3p0.npz", "main_ff_zs3p0.npz", "order0_zs3p0.npz", "order0_comparison_zs3p0.npz"}
    assert set(manifest["outputs"]) == names and {p.name for p in out.iterdir()} == names | {"export_manifest.json"}
    for name, item in manifest["outputs"].items():
        assert item["sha256"] == drv.digest(out / name)
    for name, accepted in (("main_fk_zs3p0.npz", drv.ACCEPTED / "products/main_fk.npz"),
                           ("main_ff_zs3p0.npz", drv.ACCEPTED / "products/main_ff.npz"),
                           ("order0_zs3p0.npz", drv.O0_FIG)):
        with zipfile.ZipFile(out / name) as mine, zipfile.ZipFile(accepted) as ref:
            assert mine.namelist() == ref.namelist()
        with np.load(out / name, allow_pickle=False) as mine, np.load(accepted, allow_pickle=False) as ref:
            for k in ref.files:
                assert mine[k].dtype.kind == ref[k].dtype.kind and mine[k].ndim == ref[k].ndim
            assert np.all(mine["t_final"] == drv.EXPECTED_LAMBDA[3.0])
    g_fk, fk = nhfp_load_channels(out / "main_fk_zs3p0.npz", 2)
    g_ff, _ = nhfp_load_channels(out / "main_ff_zs3p0.npz", 2)
    g_o0, o0 = nhfp_load_channels(out / "order0_zs3p0.npz", 0)
    assert np.array_equal(g_fk, g_ff) and np.array_equal(g_fk, g_o0) and sorted(o0) == [
        (0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    g_lad, ladder = nhfp_load_channels(drv.V1_LADDER / drv.LADDER_TOP, 2)
    idx = np.array([int(np.argmin(np.abs(g_fk - t))) for t in g_lad])
    assert np.max(np.abs(fk[(0, 0)][idx] / ladder[(0, 0)] - 1)) == 0.0
    expected = fk_production_values()
    rows = drv.reference_rows(out / "main_fk_zs3p0.npz")
    assert [expected[key(y)] for y in rows["y"]] == rows["value"].tolist()
    assert manifest["checks"]["checks"]["main_fk_vs_v1_tree_15360_bitwise"] is True
    assert manifest["order0_transform_self_test"]["pass"] is True
    assert manifest["sources"]["consistency"]["status"] == "PASS"
    assert set(manifest["order0_diagnostics"]) >= {"fine_probe", "bracket_multiz_z2p5_z3p2",
                                                   "full_channel_vs_production_kk"}
    meta = json.loads(str(np.load(out / "order0_zs3p0.npz", allow_pickle=False)["metadata_json"]))
    assert meta["z_source"] == 3.0 and meta["lambda_source_mpc"] == drv.EXPECTED_LAMBDA[3.0]
    with pytest.raises(Refusal, match="new directory"):
        drv.export_command(drv.parse_args(["export"]))
