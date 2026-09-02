#!/bin/zsh
# Phase 2, group (c): analyses. mc_fk_complete (FK Monte-Carlo markers), the
# 2026-08 revision drivers, the FK audit code and its small products, the
# multi-z sweep helper (from the talk repo) with the vendored talk cache, and
# the talk-side wrapper. Run once after groups (a) and (b).
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
DFE=driver_field_emulators
AN=sachs_sft/analyses
HERE=$PKG/reorg_2026-09/phase2
cd "$PKG"

echo "=== (c1) moves ==="
git mv $DFE/code/mc_fk_complete $AN/mc_fk_complete
git mv $DFE/code/figures_corrected $AN/revision_2026-08
mkdir -p $AN/fk_audit_2026-08/products $AN/fk_audit_2026-08/figures
for f in probe_ell_cut_seam probe_lambda_interp probe_uv_convergence redshift_fold_trends redshift_pershell sweep_ell_high_max plot_fk_variants plot_fk_variants_cl compare_fk_variants mc_crosscheck_dense verify_factorization_numeric; do git mv $DFE/code/$f.py $AN/fk_audit_2026-08/; done
git mv $DFE/code/verify_redshift_scalings.wl $DFE/code/verify_squeezed_expansion.wl $AN/fk_audit_2026-08/
git mv $DFE/code/normalisation $AN/fk_audit_2026-08/normalisation
git mv $DFE/code/reduced_shear $AN/fk_audit_2026-08/reduced_shear
git mv $DFE/code/README.md $AN/fk_audit_2026-08/README_code_2026-08.md
git mv $DFE/products/redshift $AN/fk_audit_2026-08/products/redshift
for f in cl_ratio_kk cl_ratio_kk_conv15360 weight_identity solveQ_roundtrip mc_fixedres_s4 mc_fixedres_s8 mc_fixedres_s16 mc_plateau_seeds mc_sigma_scan_s4 mc_sigma_scan_s8 mc_sigma_scan_s16 mc_smallgamma_s4 mc_smallgamma_s8 mc_smallgamma_s16 reduced_shear_zs5; do git mv $DFE/products/$f.npz $AN/fk_audit_2026-08/products/; done
git mv $DFE/products/grid_comparison_table.txt $DFE/products/uv_convergence_bands.txt $AN/fk_audit_2026-08/products/
for f in cutoff_vertex cutoff_vertex_tree cutoff_study_partial; do mv $DFE/products/$f.npz $AN/fk_audit_2026-08/products/; done
for f in fk_variant_effects fk_variants_cl fk_variants_kappa; do git mv $DFE/figures/$f.pdf $DFE/figures/$f.png $AN/fk_audit_2026-08/figures/; done
mkdir -p $AN/analysis3/inputs
cp ../talk/scripts/run_multiz_components.py $AN/analysis3/multiz_sweep.py
cp ../talk/assets/figures/_data/multiz_components.npz $AN/analysis3/inputs/multiz_components_talk_2026-06-10.npz

echo "=== (c2) patches ==="
$PY $HERE/p2_patch.py c
$PY - <<'PYEOF'
from pathlib import Path
p = Path(".gitignore"); s = p.read_text()
add = "\n# vertex-level cutoff study (audit-note data, 6 MB, regenerable from the band pieces)\n**/fk_audit_2026-08/products/cutoff_*.npz\n"
if "fk_audit_2026-08/products/cutoff_" not in s:
    p.write_text(s + add); print("gitignore: cutoff_* audit products")
PYEOF
cat > $AN/mc_fk_complete/numbers_2to12.py <<'EOF'
"""FK as a fraction of Order-0 over 2'-12', the range current cosmic-shear
analyses actually use: xi_kappa at z_s = 5 and at the five multi-z planes.

Paths resolve relative to this file (the script moved into the analysis package
on 2026-09-02; it used to be run from the repository root with relative strings).
"""
import math
from pathlib import Path

import numpy as np

SACHS = Path(__file__).resolve().parents[2]
O0 = SACHS / "sftwick_outputs" / "2PCF" / "C_corr_op_O0" / "xi_C_corr_op_O0.npz"
FK = (SACHS / "sftwick_outputs" / "2PCF" / "C_corr_op_K_limber_FK_cut15360_permfix"
      / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz")
MULTIZ = SACHS / "analyses" / "analysis3" / "outputs" / "multiz_kappa_2pcf_5z.npz"
GS = np.array([2.0, 5.0, 8.0, 12.0])


def fold(path, order, a, b):
    # allow_pickle: the sweep npz files are this pipeline's own outputs, whose
    # direction vectors are stored as object arrays (every loader in the package
    # reads them this way); nothing here is downloaded or user supplied.
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == a) & (d["b"] == b) & (d["order"] == order)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    o = np.argsort(g)
    return np.asarray(g)[o], np.asarray(d["value"], float)[m][o]


print("=== z_s = 5, the main analysis: FK / Order-0 [%] ===")
go, o0 = fold(O0, 0, 0, 0)
gk, fk = fold(FK, 2, 0, 0)
print(f"{'gamma':>7} {'FK':>12} {'FK/O0 %':>9}")
r5 = []
for t in GS:
    a = float(np.interp(t, gk, fk)); b = float(np.interp(t, go, o0))
    r5.append(100 * a / b)
    print(f"{t:>7.1f} {a:>12.4e} {r5[-1]:>9.2f}")
print(f"  -> kappa-kappa over 2'-12': {min(r5):.2f}% to {max(r5):.2f}%")

print("\n=== multi-z, kappa-kappa: FK / Order-0 [%] at each source plane ===")
d = np.load(MULTIZ, allow_pickle=True)
z, g = np.asarray(d["z"], float), np.asarray(d["gamma"], float)
print(f"{'z_s':>5} {'2 arcmin':>10} {'5':>8} {'8':>8} {'12':>8}   {'range':>14}")
for i, zi in enumerate(z):
    r = [100 * float(np.interp(x, g, d["fk"][i])) / float(np.interp(x, g, d["o0"][i])) for x in GS]
    print(f"{zi:>5.1f} " + " ".join(f"{v:>8.2f}" for v in r) + f"   {min(r):.2f}-{max(r):.2f}%")
EOF
cat > $AN/analysis3/inputs/PROVENANCE.md <<'EOF'
# analysis3/inputs

* `multiz_components_talk_2026-06-10.npz`: byte-identical copy of
  `talk/assets/figures/_data/multiz_components.npz` (talk repository, 2026-06-10 06:15), the 20-plane
  sweep (z_s = 0.5 to 5, keys z, gamma, o0, ff, fk, lam) that `run_multiz_components.py` produced for
  the conference-talk animation with the production configs. The Order-0 planes of
  `outputs/multiz_kappa_2pcf_5z.npz`, and hence of the deployed figure 5, are PCHIP-in-z interpolants
  of its `o0` array at z_s = 1, 1.7, 2.5, 3.2, 4 (maximum relative difference 0.0, checked 2026-09-02).
  Its `ff` planes are the talk-era placeholder (z = 5 FF scaled per gamma) and its `fk` planes are the
  corner-frozen cut1000 vertex; neither is used by the paper any more.
EOF
cat > $AN/fk_audit_2026-08/README.md <<'EOF'
# fk_audit_2026-08

Probe and analysis scripts of the 2026-08-25/26 FK audit (`docs/fk_audit_2026-08/notes/*.md`,
`docs/fk_audit_2026-08/note/main.pdf`), moved here from `driver_field_emulators/code/` on 2026-09-02
with their small data products (`products/`) and figures (`figures/`). None of them feeds a
manuscript figure or number; they document why the August revision was made. Most need the canoes
virtualenv with `PYTHONPATH=$CANOES/src:../../callables/kappa3_vertex` (for `b_model`); each
docstring says which. `README_code_2026-08.md` is the folder README of the August code tree.
EOF
cat > $AN/revision_2026-08/README.md <<'EOF'
# revision_2026-08

The drivers that produced the 2026-08 revision figures by importing the paper's own generators and
rebinding their inputs (moved from `driver_field_emulators/code/figures_corrected/` on 2026-09-02):

| driver | what it made |
|---|---|
| `regenerate.py` | figures 2, 3, 4 with the FK input rebound to the cut15360 perm-aware sweep |
| `regenerate_multiz.py --reuse` | the FK planes of `analysis3/outputs/multiz_kappa_2pcf_5z.npz` (figure 5) |
| `run_ff_one_lambda.py`, `run_ff_multiz_serial.py` | the five real FF planes of figure 5 (`sftwick_outputs/2PCF/multiz/`) |
| `make_val_figure.py` | figure 6 with the corrected FK line and the `mc_fk_complete` markers |
| `regenerate_zeta_slices.py` | figure 7 at the converged cutoff (`equal_time_limber/outputs/zeta_bands_cut15360*.npz`) |
| `regenerate_mc.py`, `extend_ff_markers*.py` | superseded Monte-Carlo routes, kept for the record |

Since the generators now default to the corrected inputs (Phase 2 of the 2026-09 reorganisation),
`regenerate.py` is no longer needed to reproduce figures 2 to 4; it is kept as the record of how the
deployed files were made. Outputs go to `outputs/` here.
EOF
cat > ../talk/scripts/run_multiz_components.py <<'EOF'
"""Thin wrapper: the multi-redshift sweep code lives in the paper's analysis
package since 2026-09-02 (``SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/multiz_sweep.py``).

This module re-exports everything the talk scripts import (``extract_kk``,
``LAM_BASELINE_Z5``, ``BOUNDARY_NPZ``, ``CONFIGS``, ``build_lambda_of_z``,
``run_vertex_sweep``) and, when run as a script, writes the talk's own cache
``talk/assets/figures/_data/multiz_components.npz`` exactly as before.
"""
from __future__ import annotations

import sys
from pathlib import Path

_A3 = (Path(__file__).resolve().parents[2] / "SFT-lensing-paper-analyses" / "sachs_sft"
       / "analyses" / "analysis3")
if str(_A3) not in sys.path:
    sys.path.insert(0, str(_A3))

import multiz_sweep as _m  # noqa: E402
from multiz_sweep import *  # noqa: E402,F401,F403

OUT_NPZ = Path(__file__).resolve().parents[1] / "assets" / "figures" / "_data" / "multiz_components.npz"
_m.OUT_NPZ = OUT_NPZ

if __name__ == "__main__":
    _m.main()
EOF

echo "=== (c3) commits ==="
git add .gitignore $AN/mc_fk_complete $AN/revision_2026-08 $AN/fk_audit_2026-08 $AN/analysis3/multiz_sweep.py $AN/analysis3/inputs $AN/analysis3/compute_multiz_kappa_2pcf.py $AN/analysis3/plot_multiz_kappa_2x5.py $AN/mc_sachs_2pt/fig_xi_channels.py $DFE
git commit -q -m "refactor: move the FK Monte-Carlo, the 2026-08 revision drivers and the FK audit into sachs_sft/analyses; multi-z sweep helper and the talk cache vendored into analysis3; path arithmetic patched" -m "$CO"
echo "committed group c: $(git rev-parse --short HEAD)"
cd ../talk && git add scripts/run_multiz_components.py && git commit -q -m "refactor: run_multiz_components.py is a wrapper around the paper package's multiz_sweep.py" -m "$CO" && echo "committed talk wrapper: $(git rev-parse --short HEAD)"
cd "$PKG"; echo "=== status ==="; git status --short | head -8
