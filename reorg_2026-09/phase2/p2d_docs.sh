#!/bin/zsh
# Phase 2, group (d): documents. The August audit notes and LaTeX note, the
# handover/summary/session files, the agent briefs, the arXiv PDFs (gitignored)
# and the private editor correspondence (gitignored) move under docs/. The two
# obsolete June analysis3 notes and the withdrawn FK Monte-Carlo note go
# (history at tag pre-reorg-2026-09; untracked ones to the Trash folder).
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
T=$HOME/.Trash/STF_lensing_reorg_2026-09-02
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
DFE=driver_field_emulators
cd "$PKG"

echo "=== (d1) audit notes, LaTeX note, handover files ==="
mkdir -p docs/fk_audit_2026-08 docs/agent_briefs docs/papers docs/private
git mv $DFE/notes docs/fk_audit_2026-08/notes
git mv $DFE/note docs/fk_audit_2026-08/note
rm -f docs/fk_audit_2026-08/note/main.aux docs/fk_audit_2026-08/note/main.log docs/fk_audit_2026-08/note/main.fdb_latexmk docs/fk_audit_2026-08/note/main.fls docs/fk_audit_2026-08/note/main.out docs/fk_audit_2026-08/note/main.toc
rm -rf docs/fk_audit_2026-08/note/.benchmarks
git mv $DFE/README.md docs/fk_audit_2026-08/README_driver_field_emulators_2026-08.md
git mv $DFE/SUMMARY.md $DFE/HANDOVER.md $DFE/SESSION_STATE_2026-08-26.md $DFE/CHANGES_OUTSIDE_THIS_FOLDER.md docs/fk_audit_2026-08/

echo "=== (d2) agent briefs (were gitignored), papers, private ==="
mv $DFE/PROMPT_normalisation_study.md docs/agent_briefs/dfe_PROMPT_normalisation_study.md
mv $DFE/PROMPT_placement_complete_FK_MC.md docs/agent_briefs/dfe_PROMPT_placement_complete_FK_MC.md
mv sachs_sft/sftwick_outputs/PROMPT_run_2PCF_three_runs.md docs/agent_briefs/sftwick_outputs_PROMPT_run_2PCF_three_runs.md
mv sachs_sft/analyses/shear3pcf_fastnc/PROMPT.md docs/agent_briefs/shear3pcf_fastnc_PROMPT.md
mv sachs_sft/analyses/shear3pcf_fastnc/PROMPT_J2_RESIDUAL.md docs/agent_briefs/shear3pcf_fastnc_PROMPT_J2_RESIDUAL.md
mv $DFE/papers docs/papers/arxiv_2026-08
mv $DFE/revision_plan docs/private/revision_plan_2026-08
rm -f docs/private/revision_plan_2026-08/letter.aux docs/private/revision_plan_2026-08/letter.fdb_latexmk docs/private/revision_plan_2026-08/letter.fls docs/private/revision_plan_2026-08/letter.log docs/private/revision_plan_2026-08/letter.out docs/private/revision_plan_2026-08/plan.aux docs/private/revision_plan_2026-08/plan.fdb_latexmk docs/private/revision_plan_2026-08/plan.fls docs/private/revision_plan_2026-08/plan.log docs/private/revision_plan_2026-08/plan.out

echo "=== (d3) obsolete notes ==="
git rm -q sachs_sft/analyses/analysis3/FK_GAMMA_TRENDS_ACROSS_PANELS.md sachs_sft/analyses/analysis3/WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md
mkdir -p "$T/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt"
[ -f sachs_sft/analyses/mc_sachs_2pt/FK_NOTES.md ] && mv sachs_sft/analyses/mc_sachs_2pt/FK_NOTES.md "$T/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/"

echo "=== (d4) gitignore and index ==="
$PY - <<'PYEOF'
from pathlib import Path
p = Path(".gitignore"); s = p.read_text()
add = """
# docs/: the agent briefs are tracked here although PROMPT*.md is ignored elsewhere;
# the arXiv PDFs and the editor correspondence stay out of git.
!docs/agent_briefs/*.md
docs/papers/
docs/private/
"""
if "docs/private/" not in s:
    p.write_text(s + add); print("gitignore: docs rules")
PYEOF
cat > docs/README.md <<'EOF'
# docs/

| folder | what it is |
|---|---|
| `fk_audit_2026-08/` | the 2026-08-25/28 FK audit that produced the revision: `notes/` (one file per finding), `note/` (the 22-page LaTeX audit note, `main.pdf`), `SUMMARY.md`, `HANDOVER.md`, `SESSION_STATE_2026-08-26.md`, `CHANGES_OUTSIDE_THIS_FOLDER.md` (the audit trail of every deployment out of the old `driver_field_emulators/` tree), `README_driver_field_emulators_2026-08.md`. Historical: several absolute FK numbers in the notes are cut1000 values; the dated correction blocks (Phase 3 of the 2026-09 reorganisation) say which. |
| `agent_briefs/` | self-contained prompts handed to agent sessions (historical). |
| `onepz4_localization.md`, `pz4_fix_applied.md` | design records of the 2026-06-09 (1+z)^4 fix. |
| `papers/` | arXiv PDFs consulted in the audit (gitignored, local only). |
| `private/` | editor correspondence and the revision plan (gitignored, local only). |
EOF
git add .gitignore docs sachs_sft/sftwick_outputs sachs_sft/analyses/shear3pcf_fastnc $DFE
git commit -q -m "docs: move the 2026-08 FK audit notes, LaTeX note and handover files under docs/fk_audit_2026-08; collect the agent briefs; retire the two obsolete June analysis3 notes" -m "$CO"
echo "committed group d: $(git rev-parse --short HEAD)"
echo "=== what is left under the imported driver_field_emulators copy ==="
git ls-files $DFE | wc -l | xargs echo "tracked files:"; du -sh $DFE 2>/dev/null; ls $DFE
echo "=== status ==="; git status --short | head -5
