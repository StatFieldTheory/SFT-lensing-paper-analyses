#!/bin/zsh
# Phase 1, analyses repo: (1) gen_sftwick_configs reads a local copy of its base
# geometry, (2) verbatim import of driver_field_emulators by copy with an md5
# manifest, (3) gitignore rules for the large, exploratory and private material,
# (4) start tracking figure_manifest.md and the 2026-06-17 cleanup docs.
# Idempotent: every step checks whether it already happened.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
cd "$R/SFT-lensing-paper-analyses"

echo "=== 1. gen_sftwick_configs: local base config ==="
SRC=_archive/canoes_pipeline/scripts/config_FK_LIMBER_z5.yaml
DST=sachs_sft/scripts/inputs/config_FK_LIMBER_z5_archived_2026-05.yaml
if [ ! -f "$DST" ]; then cp "$SRC" "$DST"; echo "copied base config"; fi
$PY - <<'PYEOF'
from pathlib import Path
p = Path("sachs_sft/scripts/gen_sftwick_configs.py"); s = p.read_text()
old_a = 'ARCHIVE = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/_archive/canoes_pipeline")\n'
old_b = 'BASE = ARCHIVE / "scripts" / "config_FK_LIMBER_z5.yaml"\n'
if old_a in s and old_b in s:
    s = s.replace(old_a, '# The base geometry travels with this script since 2026-09-02; the\n# canoes_pipeline archive it used to be read from was retired.\n')
    s = s.replace(old_b, 'BASE = Path(__file__).resolve().parent / "inputs" / "config_FK_LIMBER_z5_archived_2026-05.yaml"\n')
    p.write_text(s); print("patched gen_sftwick_configs.py")
else:
    print("gen_sftwick_configs.py already patched")
PYEOF
git add sachs_sft/scripts/gen_sftwick_configs.py "$DST"
git diff --cached --quiet || git commit -m "refactor: gen_sftwick_configs reads its base geometry from scripts/inputs instead of the retired archive" -m "$CO" | tail -1

echo "=== 2. import driver_field_emulators by copy ==="
if [ ! -d driver_field_emulators ]; then cp -a ../driver_field_emulators ./driver_field_emulators; echo "copied"; fi
if diff -rq ../driver_field_emulators ./driver_field_emulators > /dev/null; then echo "copy identical to the original"; else echo "WARNING: copy differs from the original"; fi
mkdir -p reorg_2026-09
MAN=reorg_2026-09/dfe_manifest_2026-09-02.txt
if [ ! -s "$MAN" ]; then
  find ../driver_field_emulators -type f -not -path '*/__pycache__/*' -not -name .DS_Store | sort | while read f; do
    printf '%s\t%s\t%s\n' "$(md5 -q "$f")" "$(stat -f %z "$f")" "${f#../}"
  done > "$MAN"
fi
echo "manifest lines: $(wc -l < "$MAN")"

echo "=== 3. gitignore ==="
$PY - <<'PYEOF'
from pathlib import Path
p = Path(".gitignore"); s = p.read_text()
for line in ("/CLEANUP_HANDOFF.md\n", "/CLEANUP_REPORT.md\n", "/figure_manifest.md\n", "/docs/\n"):
    s = s.replace(line, "")
s = s.replace("# NOTE: /docs/ is anchored to the repo root so it does NOT catch the tracked\n# proof docs under mathematica/docs/.\n", "")
marker = "# driver_field_emulators, imported 2026-09-02"
if marker not in s:
    s += """
# driver_field_emulators, imported 2026-09-02 during the reproducibility reorganisation.
# Large regenerable products, exploratory variants and private material stay out of git;
# every excluded file is listed with its md5 in reorg_2026-09/dfe_manifest_2026-09-02.txt.
driver_field_emulators/products/pieces/
driver_field_emulators/products/_ff_caches/
driver_field_emulators/products/_ff_dt1/
driver_field_emulators/products/cutoff_vertex.npz
driver_field_emulators/products/cutoff_vertex_tree.npz
driver_field_emulators/products/cutoff_study_partial.npz
driver_field_emulators/products/_variant_collapsed_dense_*/
driver_field_emulators/products/collapsed_dense_*
driver_field_emulators/products/collapsed_triples_dense.npz
driver_field_emulators/products/_low_cache_collapsed_dense.npz
driver_field_emulators/products/table_tree_cut1000_r*
driver_field_emulators/products/_variant_table_tree_cut1000_*/
driver_field_emulators/products/table_permclosed_cut1000*
driver_field_emulators/products/_variant_table_permclosed_cut1000_*/
driver_field_emulators/products/table_reproduce_deployed*
driver_field_emulators/products/_variant_table_reproduce_deployed/
driver_field_emulators/products/triples_full_r1.npz
driver_field_emulators/products/triples_full_r2.npz
driver_field_emulators/products/triples_permclosed_r1.npz
driver_field_emulators/products/triples_permclosed_r2.npz
driver_field_emulators/products/triples_collapsed_only.npz
driver_field_emulators/products/triples_production_only.npz
driver_field_emulators/products/triples_permutation_probe.npz
driver_field_emulators/products/xi_FF_dt1.npz
driver_field_emulators/products/_l2_lambda_grid_*.npz
driver_field_emulators/products/xi_kappa_channels_preRevision.pdf
driver_field_emulators/products/multiz_kappa_2pcf_5z_preRevision.npz
driver_field_emulators/figures/corrected/
driver_field_emulators/papers/
driver_field_emulators/revision_plan/
driver_field_emulators/code/mc_fk_complete/psd_fix/**/*.npz
driver_field_emulators/note/*.aux
driver_field_emulators/note/*.log
driver_field_emulators/note/*.fdb_latexmk
driver_field_emulators/note/*.fls
driver_field_emulators/note/*.out
driver_field_emulators/note/*.toc
**/.benchmarks/
"""
p.write_text(s); print("gitignore updated")
PYEOF
git add .gitignore figure_manifest.md CLEANUP_REPORT.md CLEANUP_HANDOFF.md docs
git diff --cached --quiet || git commit -m "docs: track figure_manifest, the 2026-06-17 cleanup report and handoff, and the design records" -m "$CO" | tail -1

echo "=== 4. import commit ==="
git add driver_field_emulators "$MAN"
echo "staged files: $(git diff --cached --name-only | wc -l)"
git diff --cached --name-only | while read f; do [ -f "$f" ] && stat -f %z "$f"; done | awk '{s+=$1} END {printf "staged bytes: %.1f MB\n", s/1e6}'
git diff --cached --quiet || git commit -m "chore: import driver_field_emulators as of 2026-09-02 with the paper-dependent products, code and documentation; large regenerable, exploratory and private material gitignored and listed in the manifest" -m "$CO" | tail -1
if ! git tag -l | grep -q '^dfe-imported-2026-09$'; then git tag dfe-imported-2026-09; echo "tagged dfe-imported-2026-09"; fi
echo "=== log ==="; git log --oneline -6
echo "=== status ==="; git status --short | head
