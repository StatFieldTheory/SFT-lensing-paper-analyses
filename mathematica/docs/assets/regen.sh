#!/usr/bin/env bash
# Regenerate SVG diagrams from the .mmd sources in this folder.
# Usage:  bash scripts/docs/assets/regen.sh
set -euo pipefail
cd "$(dirname "$0")"
for src in _*.mmd; do
  out="${src#_}"
  out="${out%.mmd}.svg"
  npx --yes -p @mermaid-js/mermaid-cli mmdc -i "$src" -o "$out"
done
