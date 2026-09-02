#!/bin/zsh
# Phase 1, paper repo: rewrite .gitignore (no global markdown ignore; agent files
# and every build output ignored explicitly) and remove the archived and orphan
# figure PDFs, which stay in history at tag pre-reorg-2026-09. Idempotent.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
cd "$R"
cat > .gitignore <<'EOF'
# LaTeX auxiliary files
.DS_Store
*.aux
*.bbl
*.bcf
*.blg
*.fdb_latexmk
*.fls
*.lof
*.log
*.lot
*.out
*.run.xml
*.synctex.gz
*.synctex(busy)
*.toc
*.nav
*.snm
*.vrb
*.idx
*.ind
*.ilg
*.glo
*.gls
*.glg
*.acn
*.acr
*.alg
*.ist
*.xdy
*.mw
*.mtc*
*.mlf*
*.mlt*
*Notes.bib
main.pdf
main_clean.pdf
main *.*
meascol.*

# latexmk
.latexmkrc

# Build output (if using a build directory)
build/
out/
output*/
talk/

inputs/
logs/
.mplconfig/

# Figure source/preview duplicates; the paper includes PDFs explicitly.
figures/*.png

# Cache and benchmark directories
*cache*/
*benchmark*/

**/.cache_propagators/**
**/.cache_expand/**
**/.pytest_cache/**
**/__pycache__/**

# Agent instruction files and local tooling; not manuscript content.
/CLAUDE.md
/AGENTS.md
.claude
.playwright-mcp/

# Nested repositories and folders that are not part of the manuscript.
SFT-lensing-paper-analyses/
archived/
PerStat/
interview_erc/
driver_field_emulators/
EOF
git add .gitignore
git diff --cached --quiet || git commit -m "chore: stop ignoring markdown globally; ignore the agent instruction files and every build output explicitly" -m "$CO" | tail -1
for p in figures/archived_2026-08-26 figures/archived_2026-08-28 figures/appendix_3cumulant_fastnc.pdf figures/cumulant_hierarchy.pdf; do
  if git ls-files --error-unmatch "$p" > /dev/null 2>&1; then git rm -r "$p" > /dev/null; echo "removed $p"; fi
done
git diff --cached --quiet || git commit -m "chore: remove the archived and orphan figure PDFs; all remain in history at tag pre-reorg-2026-09" -m "$CO" | tail -1
echo "=== log ==="; git log --oneline -4
echo "=== status ==="; git status --short
echo "figure pdf count: $(ls figures/*.pdf | wc -l)"
