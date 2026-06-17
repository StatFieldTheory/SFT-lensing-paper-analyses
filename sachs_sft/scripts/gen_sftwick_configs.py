"""Generate the 5 sachs_sft sft-wick L2 configs from the verified +3.08e-5 baseline.

Loads the archived config_FK_LIMBER_z5.yaml (the +3.08e-5 geometry), then for each
run produces a complete standalone YAML with ONLY the wiring changed (which C, which
K, FF vs FK). The shared geometry (40-pt positions_grid, t_final=2313.029,
component_pairs, n_gauss=24) is therefore byte-identical across all configs — any
diff between them is the intended physics wiring, nothing else.

Callable modules point at ABSOLUTE sachs_sft paths (the new self-describing callables
with bundled tables), sidestepping sft-wick's dual path-resolution (CWD vs YAML-parent).
"""
from __future__ import annotations
import copy
from pathlib import Path
import yaml

ARCHIVE = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/_archive/canoes_pipeline")
BASE = ARCHIVE / "scripts" / "config_FK_LIMBER_z5.yaml"
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

with BASE.open() as f:
    base = yaml.safe_load(f)


def make_config(*, c_module, k_module=None, k_already_R, k_equal_time,
                k_vectorized=False, vertex_types, orders, out_npz, cache_label):
    cfg = copy.deepcopy(base)
    # --- system.linear (response R) + noise (kappa2) → absolute vendored scripts ---
    cfg["system"]["linear"]["R_time_module"] = str(D_CALLABLE)
    cfg["system"]["noise"]["kappa2"]["module"] = str(KAPPA2_CALLABLE)
    cfg["system"]["vertices"][0]["coupling_path"] = str(F_TENSOR)
    # --- C propagator (absolute sachs_sft callable) ---
    # c_closed_form_vectorized MUST be False: the sachs_sft C_fn wrappers expose a
    # SCALAR coupling_fn (C_fn(n1,t1,n2,t2) -> (3,3)); they do NOT re-export the
    # underlying .batch attribute, so a vectorized array-t call would crash on
    # float(t_array). The C propagator is cached once, so the per-sample loop cost
    # is paid once.
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

for sub, kw in runs.items():
    cfg = make_config(**kw)
    dest = OUT / sub / "config_L2.yaml"
    with dest.open("w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False, width=100)
    print(f"wrote {dest}")
