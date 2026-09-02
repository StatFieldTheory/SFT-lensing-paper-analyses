#!/bin/sh
# Assemble and fold the corrected vertex at a series of multipole cutoffs.
#
# The bands are disjoint in max leg, so the partial sum up to each octave
# edge IS the vertex built with that cutoff. Assembly is therefore cheap and
# only the fold costs anything. Folds run one at a time on purpose: sft-wick
# forks with loky and each worker reloads the inputs, so several folds at
# once is how a machine gets lost.
#
# Usage:  sh sweep_cutoff.sh <model> <cutoff> [<cutoff> ...]
set -e
CAN=/Users/zzhang/projects/angular_statistics/canoes
P="$(cd "$(dirname "$0")" && pwd)/products"
SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
MODEL="$1"; shift
for CUT in "$@"; do
  OUT="$P/table_${MODEL}_cut${CUT}_r4.npz"
  if [ ! -f "$OUT" ]; then
    PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1 \
      "$CAN/.venv/bin/python" -u assemble.py \
        --low "$P/pieces/low_r4.npz" \
        --bands "$P/pieces/band_${MODEL}_r4_"*.npz \
        --cutoff "$CUT" --out "$OUT"
  fi
  if [ ! -f "$P/table_${MODEL}_cut${CUT}_r4_xi.npz" ]; then
    "$SFT" "$(dirname "$0")/run_fk_variant.py" "$OUT" --n-jobs ${N_JOBS:-4}
  fi
done
echo "[sweep] $MODEL done for cutoffs: $*"
