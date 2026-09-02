#!/bin/zsh
# Phase 1, talk repo (decision D12, part 1 of 3): stop tracking the compiled
# __pycache__ files and add a .gitignore. The wrapper import (part 3) lands in
# Phase 2 once multiz_sweep.py exists in the analysis package. Idempotent.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
cd "$R/talk"
if git ls-files | grep -q '__pycache__'; then
  git ls-files | grep '__pycache__' | while read f; do git rm --cached "$f" > /dev/null; done
  echo "untracked the compiled files"
fi
if [ ! -f .gitignore ]; then
cat > .gitignore <<'EOF'
__pycache__/
*.pyc
.ruff_cache/
.DS_Store
EOF
  git add .gitignore
  echo "added .gitignore"
fi
git diff --cached --quiet || git commit -m "chore: stop tracking compiled python files and add a gitignore" -m "$CO" | tail -1
echo "=== log ==="; git log --oneline -3
echo "=== status ==="; git status --short | head
