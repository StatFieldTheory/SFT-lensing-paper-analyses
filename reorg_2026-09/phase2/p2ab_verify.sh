#!/bin/zsh
# Verification after groups (a) and (b): the promoted config loads through sft-wick
# and every path it names exists; the permutation-aware tests pass from the new
# location; figure 17 regenerates from the cutoff_ladder folder and matches the
# deployed PDF; the rebuild modules import.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
S=/private/tmp/claude-501/-Users-zzhang-Documents-MyDrafts-STF-lensing/58fd681a-480f-4860-a212-3d176f65b575/scratchpad
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
SFTW=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
CAN=/Users/zzhang/projects/angular_statistics/canoes
KV=$PKG/sachs_sft/callables/kappa3_vertex
NEW=$PKG/sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix

echo "=== 1. promoted config loads through sft-wick; resolved paths exist ==="
cd $NEW && $SFTW - <<'PYEOF'
from pathlib import Path
from sft_wick.workflow.config import load_workflow_config
cfg = load_workflow_config(Path("config_L2.yaml"))
print("c_closed_form_module:", cfg.propagators.c_closed_form_module, Path(cfg.propagators.c_closed_form_module).exists())
nv = cfg.system.nonlocal_vertices[0] if hasattr(cfg.system, "nonlocal_vertices") else None
print("nonlocal vertex:", nv if nv is None else {k: (str(v)[:90]) for k, v in (nv.items() if isinstance(nv, dict) else vars(nv).items())})
print("output:", cfg.output[0].path, "| expand cache:", cfg.expand.cache_path, "| n_jobs sweep/props:", cfg.sweep.n_jobs, cfg.propagators.n_jobs)
print("t_final:", cfg.sweep.t_final_grid, "| n_gamma:", len(cfg.sweep.positions_grid["y"]) if isinstance(cfg.sweep.positions_grid, dict) else "?")
PYEOF

echo "=== 2. permutation-aware tests from the new location (canoes venv) ==="
cd $KV/equal_time_limber_cut15360_permaware && PYTHONPATH=$CAN/src $CAN/.venv/bin/python -m pytest tests -q 2>&1 | tail -3

echo "=== 3. figure 17 from the new location with default arguments ==="
mkdir -p $S/p2/cutoff && cd $KV/rebuild && MPLBACKEND=Agg $PY fig_cutoff_paper.py --out $S/p2/cutoff/fk_cutoff_convergence.pdf 2>&1 | tail -3
$PY $S/pdfcmp.py $R/figures/fk_cutoff_convergence.pdf $S/p2/cutoff/fk_cutoff_convergence.pdf fig17_new_location | head -1

echo "=== 4. rebuild modules import; vendored lambda grid readable; b_model tests ==="
cd $KV/rebuild && PYTHONPATH=$CAN/src $CAN/.venv/bin/python - <<'PYEOF'
import importlib, sys
sys.path.insert(0, ".")
import _common
print("_REPO:", _common._REPO, "| _ETL exists:", _common._ETL.is_dir(), "| ARCHIVED_L2 exists:", _common.ARCHIVED_L2.exists())
import numpy as np
lam, lam_h = _common.load_lambda_grid(_common.ARCHIVED_L2, 0.6711)
print("lambda grid:", lam.size, "shells,", round(float(lam[0]), 3), "..", round(float(lam[-1]), 3))
for m in ("observables", "assemble", "build_low", "build_band", "run_queue", "fk_kernel_crosscheck", "mean_kappa_z", "fig_redshift", "probe_weight_identity", "run_fk_variant"):
    importlib.import_module(m); print("import ok:", m)
PYEOF
cd $KV && PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q 2>&1 | tail -2

echo "=== 5. grep for stale layout strings in the moved code ==="
cd $PKG && grep -rn -E 'driver_field_emulators|_archive/|parents\[3\]' --include='*.py' --include='*.sh' $KV/rebuild $KV/equal_time_limber_cut15360_permaware $KV/b_model sachs_sft/scripts/run_FF_single.py | grep -v -E '^\S+:\s*(#|"""|\.\.\.)|Usage|history|moved|came from|2026-09-02|record' | cut -c1-160 || echo "(none)"
echo "=== git status ==="; git status --short | head -5
