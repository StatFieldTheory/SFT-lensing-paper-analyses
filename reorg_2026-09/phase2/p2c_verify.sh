#!/bin/zsh
# Verification after group (c): the moved drivers reproduce figures 2 and 5 from
# their new locations, the Monte-Carlo bootstrap resolves, the number scripts run,
# the multi-z helper imports in both the package and the talk wrapper.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
S=/private/tmp/claude-501/-Users-zzhang-Documents-MyDrafts-STF-lensing/58fd681a-480f-4860-a212-3d176f65b575/scratchpad
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
SFTW=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
CAN=/Users/zzhang/projects/angular_statistics/canoes
AN=$PKG/sachs_sft/analyses
FK=$PKG/sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz

echo "=== 1. regenerate.py from revision_2026-08 reproduces figures 2, 3, 4 ==="
mkdir -p $S/p2/rev && cd $AN/revision_2026-08 && MPLBACKEND=Agg $PY regenerate.py --fk-npz $FK --figure nlo cl eb --out-dir $S/p2/rev --suffix check 2>&1 | grep -E '^\[fig\]' | head -4
$PY $S/pdfcmp.py $R/figures/analysis3_NLO_FFFK.pdf $S/p2/rev/analysis3_nlo_O0_FF_FK_check.pdf fig2_via_regenerate | head -1
$PY $S/pdfcmp.py $R/figures/analysis3_cl_full_vs_O0.pdf $S/p2/rev/analysis3_cl_O0_FF_FK_check.pdf fig3_via_regenerate | head -1
$PY $S/pdfcmp.py $R/figures/cl_EB_polarization.pdf $S/p2/rev/cl_EB_polarization_check.pdf fig4_via_regenerate | head -1

echo "=== 2. figure 5 plot stage with the relocated FF planes ==="
cd $S && P0_SCRATCH=$S/p2/env_pyccl MPLBACKEND=Agg $PY $PKG/reorg_2026-09/phase0/harness/p0_regen.py --figure multiz --mode standalone 2>&1 | grep -E 'harness|multiz\]' | tail -2
$PY $S/pdfcmp.py $R/figures/multiz_kappa_xi_cl.pdf $S/p2/env_pyccl/multiz/multiz_kappa_xi_cl_2x5.pdf fig5_relocated_ff | head -1

echo "=== 3. mc_fk_complete: bootstrap resolves, headline and the 2-12 numbers ==="
cd $AN/mc_fk_complete && PYTHONPATH=$CAN/src $PY - <<'PYEOF'
import _bootstrap as b
print("TABLE exists:", b.TABLE.exists(), "| FOLD exists:", b.FOLD.exists(), "| MC_DIR exists:", b.MC_DIR.is_dir(), "| FIX_DIR exists:", b.FIX_DIR.is_dir())
ds, bg, core = b.wire()
print("wired: driver_stats _k3 =", type(ds._k3).__name__, "table:", ds._k3.TABLE_PATH.name)
PYEOF
$PY headline.py 2>&1 | tail -2
$PY numbers_2to12.py 2>&1 | tail -9

echo "=== 4. multi-z helper imports (package and talk wrapper, sft-wick env) ==="
cd $AN/analysis3 && $SFTW -c "import multiz_sweep as m; print('multiz_sweep OK; OUT_NPZ =', m.OUT_NPZ.relative_to(m.REPO)); print('CONFIGS exist:', all(p.exists() for p in m.CONFIGS.values()))"
cd $R/talk/scripts && $SFTW -c "import run_multiz_components as r; print('talk wrapper OK; OUT_NPZ =', r.OUT_NPZ); print('re-exports:', all(hasattr(r, n) for n in ('extract_kk','LAM_BASELINE_Z5','BOUNDARY_NPZ','CONFIGS','build_lambda_of_z','run_vertex_sweep')))"
cd $AN/analysis3 && $SFTW -c "import ast,sys; ast.parse(open('compute_multiz_kappa_2pcf.py').read()); print('compute_multiz_kappa_2pcf.py parses')"

echo "=== 5. audit scripts: syntax and the two path-patched imports ==="
cd $AN/fk_audit_2026-08 && for f in *.py normalisation/*.py reduced_shear/*.py; do $PY -c "import ast; ast.parse(open('$f').read())" || echo "SYNTAX ERROR $f"; done; echo "all audit scripts parse"
$PY -c "
from pathlib import Path
import re
for f in ('redshift_fold_trends.py','redshift_pershell.py','plot_fk_variants.py','plot_fk_variants_cl.py','mc_crosscheck_dense.py'):
    s = Path(f).read_text(); print(f, '-> stale?', 'driver_field_emulators' in s or 'SFT-lensing-paper-analyses' in s)
"
echo "=== 6. leftovers under the imported driver_field_emulators copy (tracked) ==="
cd $PKG && git ls-files driver_field_emulators | sed 's#/[^/]*$##' | sort | uniq -c | sort -rn | head -12
echo "=== git status ==="; git status --short | head -5; git -C $R/talk status --short | head -3
