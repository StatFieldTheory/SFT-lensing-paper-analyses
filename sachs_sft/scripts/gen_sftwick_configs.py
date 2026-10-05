"""Generate historical five-run YAML recipes into an explicit fresh directory.

Historical geometry and wiring are preserved, not accepted current candidates.
No science callables are imported. --dry-run performs validation without writes.
"""
from __future__ import annotations
import copy, argparse, hashlib, json, os
from pathlib import Path
import yaml

# The base geometry travels with this script since 2026-09-02; the
# canoes_pipeline archive it used to be read from was retired.
BASE = Path(__file__).resolve().parent / "inputs" / "config_FK_LIMBER_z5_archived_2026-05.yaml"
SACHS = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft")
CALL = SACHS / "callables"
OUT = SACHS / "sftwick_outputs"

# Absolute callable module paths (self-describing, bundled tables)
C_CORR_OP = CALL / "C_propagator/corr_op/corr_op_C_callable.py"
C_LIMBER = CALL / "C_propagator/limber/limber_C_callable.py"
K_RCONTRACTED = CALL / "kappa3_vertex/R_contracted/windowed_squeezed_kappa3_callable.py"
K_EQT_LIMBER = CALL / "kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py"
# Shared scripts (CWD-relative for run_FF_single which chdir's to scripts/; we give absolute)
SCRIPTS = SACHS / "scripts"
D_CALLABLE = SCRIPTS / "D_callable.py"
KAPPA2_CALLABLE = SCRIPTS / "kappa2_callable.py"
F_TENSOR = SCRIPTS / "inputs" / "F_tensor.npy"

base = None  # Loaded only by explicit main(), never at import.


def make_config(*, c_module, k_module=None, k_already_R, k_equal_time,
                k_vectorized=False, vertex_types, orders, out_npz, cache_label):
    cfg = copy.deepcopy(base)
    # --- system.linear (response R) + noise (kappa2) → absolute vendored scripts ---
    cfg["system"]["linear"]["R_time_module"] = str(D_CALLABLE)
    cfg["system"]["noise"]["kappa2"]["module"] = str(KAPPA2_CALLABLE)
    cfg["system"]["vertices"][0]["coupling_path"] = str(F_TENSOR)
    # --- C propagator (absolute sachs_sft callable) ---
    # Historical scalar recipe is retained, not a claim about modern batch support.
    # Reviewed current FK/FF candidates bind C_fn_batch/vectorized true separately.
    cfg["propagators"]["c_closed_form_module"] = str(c_module)
    cfg["propagators"]["c_closed_form_attr"] = "C_fn"
    cfg["propagators"]["c_closed_form_only"] = True
    cfg["propagators"]["c_closed_form_vectorized"] = False
    # --- K vertex (only for FK runs) ---
    if k_module is None:
        cfg["system"].pop("nonlocal_vertices", None)
    else:
        # equal_time_limber exposes a vectorized coupling_fn_batch (.vectorized=True);
        # R_contracted is scalar-only (coupling_fn). Use the batch attr when available.
        coupling_attr = "coupling_fn_batch" if k_vectorized else "coupling_fn"
        k = {
            "name": "K", "order": 3,
            "coupling_module": str(k_module),
            "coupling_attr": coupling_attr,
            "coupling_vectorized": bool(k_vectorized),
            "equal_time": bool(k_equal_time),
        }
        if k_already_R:
            k["already_R_contracted"] = True
        cfg["system"]["nonlocal_vertices"] = [k]
    # --- expand orders + sweep vertex_types ---
    cfg["expand"]["orders"] = list(orders)
    cfg["expand"]["cache_path"] = str((out_npz.parent / f".cache_expand_{cache_label}"))
    cfg["sweep"]["vertex_types"] = list(vertex_types)
    cfg["sweep"]["integrate_over"] = "all"   # user requirement: always "all"
    cfg["propagators"]["cache_path"] = str((out_npz.parent / f".cache_propagators_{cache_label}"))
    # --- output (folder-local) ---
    cfg["output"] = [{"type": "npz", "path": str(out_npz)}]
    return cfg


runs = {
    # 2PCF — Order-0, C only (no K vertex)
    "2PCF/C_corr_op_O0": dict(
        c_module=C_CORR_OP, k_module=None, k_already_R=False, k_equal_time=False,
        vertex_types=["F"], orders=[0], cache_label="C_corr_op_O0",
        out_npz=OUT / "2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz"),
    "2PCF/C_limber_O0": dict(
        c_module=C_LIMBER, k_module=None, k_already_R=False, k_equal_time=False,
        vertex_types=["F"], orders=[0], cache_label="C_limber_O0",
        out_npz=OUT / "2PCF/C_limber_O0/xi_C_limber_O0.npz"),
    # 2PCF Order-[0,2] — these still compute the lensing 2-POINT function (the
    # observable is [phi_a(x), phi_b(y)] = 2-point); the F/K vertices enter only as
    # Order-2 post-Born corrections. "FF"/"FK" name the sft-wick VERTEX types in the
    # O(2) diagrams (F=local cubic, K=non-local kappa3), NOT the observable point
    # count. So they live under 2PCF/, not 3PCF/.
    "2PCF/C_corr_op_K_contracted_FF": dict(
        c_module=C_CORR_OP, k_module=K_RCONTRACTED, k_already_R=True, k_equal_time=False,
        vertex_types=["F"], orders=[0, 2], cache_label="C_corr_op_K_contracted_FF",
        out_npz=OUT / "2PCF/C_corr_op_K_contracted_FF/xi_C_corr_op_K_contracted_FF.npz"),
    "2PCF/C_corr_op_K_contracted_FK": dict(
        c_module=C_CORR_OP, k_module=K_RCONTRACTED, k_already_R=True, k_equal_time=False,
        vertex_types=["FK"], orders=[0, 2], cache_label="C_corr_op_K_contracted_FK",
        out_npz=OUT / "2PCF/C_corr_op_K_contracted_FK/xi_C_corr_op_K_contracted_FK.npz"),
    # corr_op C + equal_time_limber K (NOT R-contracted; sft-wick folds R).
    # equal_time_limber exposes a vectorized coupling_fn_batch -> use it.
    "2PCF/C_corr_op_K_limber_FK": dict(
        c_module=C_CORR_OP, k_module=K_EQT_LIMBER, k_already_R=False, k_equal_time=True,
        k_vectorized=True,
        vertex_types=["FK"], orders=[0, 2], cache_label="C_corr_op_K_limber_FK",
        out_npz=OUT / "2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz"),
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    global base, D_CALLABLE, KAPPA2_CALLABLE, F_TENSOR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--sachs-root", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    template = args.base.resolve()
    sachs = args.sachs_root.resolve()
    output = args.out_root.absolute()
    if not template.is_file() or not sachs.is_dir():
        raise ValueError("Missing base template or Sachs root")
    if os.path.lexists(output):
        raise ValueError("Output root must be fresh, including dry-run")
    if not output.parent.is_dir():
        raise ValueError("Output parent must exist")
    template_raw = template.read_bytes()
    template_digest = hashlib.sha256(template_raw).hexdigest()
    base = yaml.safe_load(template_raw)
    if not isinstance(base, dict) or any(k not in base for k in ("system", "expand", "propagators", "sweep", "output")):
        raise ValueError("Invalid five-block historical template")
    D_CALLABLE = sachs / "scripts/D_callable.py"
    KAPPA2_CALLABLE = sachs / "scripts/kappa2_callable.py"
    F_TENSOR = sachs / "scripts/inputs/F_tensor.npy"
    required = [template, D_CALLABLE, KAPPA2_CALLABLE, F_TENSOR, Path(__file__).resolve()]
    prepared = []
    for sub, old in runs.items():
        kw = dict(old)
        for key in ("c_module", "k_module"):
            if kw[key] is not None:
                kw[key] = sachs / kw[key].relative_to(SACHS)
                required.append(kw[key])
        kw["out_npz"] = output / sub / old["out_npz"].name
        cfg = make_config(**kw)
        prepared.append((output / sub / "config_L2.yaml", cfg))
    for path in required:
        if not path.is_file():
            raise ValueError("Missing explicit dependency: " + str(path))
    pins = {str(path): sha(path) for path in required if path != template}
    pins[str(template)] = template_digest
    # Bind parsed content to the initial capture, before any output reservation.
    for path, expected in pins.items():
        if sha(path) != expected:
            raise RuntimeError("Input drift before reservation: " + path)
    record = dict(recipe="historical five-run geometry and scalar C wiring; not current candidate acceptance",
                  base=str(template), out_root=str(output), pins=pins,
                  destinations=[str(dest) for dest, _ in prepared])
    if args.dry_run:
        print(json.dumps(record, indent=2))
        return 0
    output.mkdir()  # Exclusive fresh root; never reuse production run folders.
    for dest, cfg in prepared:
        dest.parent.mkdir(parents=True, exist_ok=False)
        with dest.open("x") as f:
            yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False, width=100)
    for path, expected in pins.items():
        if sha(path) != expected:
            raise RuntimeError("Input drift; retain partial output, no prepared record: " + path)
    record["config_hashes"] = {str(dest): sha(dest) for dest, _ in prepared}
    with (output / "prepared_configs.json").open("x") as f:
        json.dump(record, f, indent=2)
        f.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
