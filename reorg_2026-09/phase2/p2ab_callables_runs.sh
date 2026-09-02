#!/bin/zsh
# Phase 2, groups (a) callables and (b) run folders. Idempotent where cheap;
# run once. Every git mv keeps history; every patch is fail-loud (p2_patch.py).
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
DFE=driver_field_emulators
KV=sachs_sft/callables/kappa3_vertex
PA=$KV/equal_time_limber_cut15360_permaware
RB=$KV/rebuild
ETL=$KV/equal_time_limber
RUNS=sachs_sft/sftwick_outputs/2PCF
HERE=$PKG/reorg_2026-09/phase2
cd "$PKG"

echo "=== (a0) vendor the L2 lambda grid the rebuild and the vertex build read ==="
mkdir -p $ETL/inputs
if [ ! -f $ETL/inputs/l2_lambda_grid_z5covgrid16.npz ]; then
$PY - <<'PYEOF'
import numpy as np
from pathlib import Path
src = Path.home() / ".Trash/STF_lensing_reorg_2026-09-02/SFT-lensing-paper-analyses/_archive/canoes_pipeline/inputs/l2_callables_2026_05_26/windowed_squeezed_kappa3_lmax1000_z5_covgrid16_omega0316_h06711_physMpc.npz"
d = np.load(src, allow_pickle=False)
lam = np.asarray(d["lambda_grid"], dtype=np.float64)
out = Path("sachs_sft/callables/kappa3_vertex/equal_time_limber/inputs/l2_lambda_grid_z5covgrid16.npz")
np.savez(out, lambda_grid=lam, source=np.array(
    "lambda_grid of windowed_squeezed_kappa3_lmax1000_z5_covgrid16_omega0316_h06711_physMpc.npz "
    "(canoes_pipeline archive, built 2026-05-28, md5 in reorg_2026-09/trash_manifest_2026-09-02.txt); "
    "extracted 2026-09-02, the only array the builds read"))
print(f"lambda grid: {lam.size} shells, {lam[0]:.3f}..{lam[-1]:.3f} Mpc -> {out}")
PYEOF
fi

echo "=== (a1) moves: permutation-aware callable, rebuild, b_model, zeta bands ==="
git mv $DFE/code/callable_fixed $PA
git mv $DFE/products/table_permclosed_cut15360.npz $PA/table_permclosed.npz
git mv $DFE/code/rebuild $RB
git mv $DFE/code/make_dense_triples.py $DFE/code/make_collapsed_triples.py $DFE/code/run_fk_variant.py $RB/
mkdir -p $RB/products
for f in triples_permclosed_r4 triples_permclosed_extra triples_full_r4 fk_kernel_crosscheck mean_kappa_plateaus; do git mv $DFE/products/$f.npz $RB/products/; done
for m in tree bihalofit; do for c in 960 1920 3840 7680 15360; do git mv $DFE/products/table_${m}_cut${c}_r4.npz $RB/products/; done; done
git mv $DFE/products/logs $RB/logs
mv $DFE/products/pieces $RB/products/pieces
git mv $DFE/code/b_model $KV/b_model
mkdir -p $ETL/outputs
git mv $DFE/products/zeta_bands_cut15360.npz $DFE/products/zeta_bands_cut15360_gmax85.npz $ETL/outputs/

echo "=== (a2) patches ==="
$PY $HERE/p2_patch.py a
$PY - <<'PYEOF'
from pathlib import Path
p = Path(".gitignore"); s = p.read_text()
add = "\n# vertex-rebuild intermediates (octave-band pieces, hours of canoes builds; regenerable, listed with md5 in reorg_2026-09/dfe_manifest_2026-09-02.txt)\n**/products/pieces/\n**/_ff_caches/\n"
if "**/products/pieces/" not in s:
    p.write_text(s + add); print("gitignore: pieces and ff caches")
PYEOF
cat > $PA/PROVENANCE.md <<'EOF'
# PROVENANCE

* `perm_aware_kappa3_callable.py`, `compare_callables.py`, `README.md`, `tests/test_perm_aware.py`
  came from `driver_field_emulators/code/callable_fixed/` (2026-08-26, moved here 2026-09-02 with
  history). The callable and the README are byte-identical to the copies SFT-WL-B vendored on
  2026-08-31; `tests/test_perm_aware.py` and `compare_callables.py` had their path arithmetic
  updated for the new location (old contents at tag `dfe-imported-2026-09`).
* `table_permclosed.npz` is `driver_field_emulators/products/table_permclosed_cut15360.npz`,
  byte-identical (sha256 prefix `b400dc5b11e1033b`), renamed to the name the callable resolves.
  Built 2026-08-26 09:46 by `rebuild/assemble.py` from `rebuild/products/pieces/low_permclosed_r4.npz`
  plus the `band_tree_permclosed_*` and `band_tree_extra_*` octave windows (tree-level bispectrum,
  ell_max = 15360, permutation-closed r4 grid plus extra collapsed rows); logs in `rebuild/logs/`.
* This folder is the production vertex for every FK result in the manuscript since the 2026-08
  revision; the superseded cut1000 folder is `../equal_time_limber/` (see its DEPRECATED.md).
EOF
cat > $ETL/DEPRECATED.md <<'EOF'
# Superseded as the production vertex (2026-09-02)

`equal_time_limber_kappa3_callable.py` and `equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz`
are the cut1000, corner-frozen vertex of the June 2026 draft. They are kept byte-identical because
SFT-WL-B vendored them (2026-08-31) and because `build_equal_time_limber_table.py`,
`ell_band_decomp.py` and `plot_zeta_figures.py` here are still the live build and figure-7 code.

The production vertex of the manuscript is `../equal_time_limber_cut15360_permaware/`
(permutation-aware lookup, ell_max = 15360). The callable here sorts the query triple and rounds
the mesh key to 8 decimals, which is wrong for the spin channels (23% on zeta_TPP, 70% on
zeta_Dmod at the collapsed configurations) and cannot represent separations below 0.344 arcmin;
see `../equal_time_limber_cut15360_permaware/README.md`.

The "+3.086e-5" baseline quoted in this folder's README and docstrings is a pre-2026-06-09 stamp;
see `sftwick_outputs/2PCF/C_corr_op_K_limber_FK/README.md` for the history of that number.
EOF
git add $PA/PROVENANCE.md $ETL/DEPRECATED.md $ETL/inputs/l2_lambda_grid_z5covgrid16.npz .gitignore $RB $PA $KV/b_model $ETL
git commit -q -m "refactor: move the permutation-aware callable, the factored vertex rebuild and b_model into sachs_sft/callables/kappa3_vertex; vendor the L2 lambda grid; patch path arithmetic" -m "$CO"
echo "committed group a: $(git rev-parse --short HEAD)"

echo "=== (b1) the cut15360 perm-aware FK run folder ==="
NEW=$RUNS/C_corr_op_K_limber_FK_cut15360_permfix
mkdir -p $NEW
git mv $DFE/products/table_permclosed_cut15360_permfix_xi.npz $NEW/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz
git mv $DFE/products/_variant_table_permclosed_cut15360_permfix/config__variant_table_permclosed_cut15360_permfix.yaml $NEW/config_record_2026-08-26.yaml
$PY - <<'PYEOF'
import yaml
from pathlib import Path
src = Path("sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK/config_L2.yaml")
cfg = yaml.safe_load(src.read_text())
S = "../../../"
cfg["system"]["linear"]["R_time_module"] = S + "scripts/D_callable.py"
cfg["system"]["noise"]["kappa2"]["module"] = S + "scripts/kappa2_callable.py"
cfg["system"]["vertices"][0]["coupling_path"] = S + "scripts/inputs/F_tensor.npy"
cfg["system"]["nonlocal_vertices"][0]["coupling_module"] = S + "callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py"
cfg["expand"]["cache_path"] = ".cache_expand_C_corr_op_K_limber_FK_cut15360_permfix"
cfg["propagators"]["cache_path"] = ".cache_propagators_C_corr_op_K_limber_FK_cut15360_permfix"
cfg["propagators"]["c_closed_form_module"] = S + "callables/C_propagator/corr_op/corr_op_C_callable.py"
cfg["propagators"]["n_jobs"] = 6
cfg["sweep"]["n_jobs"] = 6
cfg["output"][0]["path"] = "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"
out = Path("sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/config_L2.yaml")
header = ("# FK sweep at ell_max = 15360 with the permutation-aware vertex (the manuscript's FK).\n"
          "# Derived 2026-09-02 from ../C_corr_op_K_limber_FK/config_L2.yaml: same 40-point gamma grid,\n"
          "# t_final, component pairs and n_gauss; only the vertex callable, the worker count (6) and the\n"
          "# paths differ. All paths are relative to this folder (module and data paths are resolved by\n"
          "# sft-wick against the YAML directory; output and cache paths by the runner, which chdirs here).\n"
          "# The verbatim config of the 2026-08-26 run is config_record_2026-08-26.yaml.\n")
out.write_text(header + yaml.safe_dump(cfg, sort_keys=False))
print("wrote", out)
PYEOF
cat > $NEW/README.md <<'EOF'
# 2PCF / C_corr_op_K_limber_FK_cut15360_permfix

**The FK two-point sweep of the manuscript** (figures 2, 3, 4, 6 and every FK number, ell_max = 15360).

* `xi_C_corr_op_K_limber_FK_cut15360_permfix.npz`: byte-identical to
  `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz`, produced 2026-08-26
  09:48 by `run_fk_variant.py table_permclosed_cut15360.npz --callable perm_aware_kappa3_callable.py
  --tag permfix --n-jobs 6` (fold time about 2.5 min, sft-wick env). Order-2 kappa-kappa at 0.5' is
  `+1.9480e-5`, 2.31% of Order-0.
* `config_L2.yaml`: the same run expressed with paths relative to this folder (vertex =
  `../../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/`); `config_record_2026-08-26.yaml`
  is the verbatim config of the August run (absolute paths of that machine).
* `PRODUCTION`: marker read by `scripts/run_FF_single.py`, which refuses to overwrite the result here
  unless `SFT_WICK_FORCE_OVERWRITE=1`.

Re-run (writes into this folder; the runner clears the two `.cache_*` directories first):

    cd sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix
    SFT_WICK_FORCE_OVERWRITE=1 SFT_WICK_CONFIG_YAML=$PWD/config_L2.yaml SFT_WICK_EXPAND_ORDERS=0,2 \
      SFT_WICK_SWEEP_N_GAUSS=24 /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
      ../../../scripts/run_FF_single.py ../../../scripts/inputs/kappa2_grid_lcut0_low_empty.npz \
      ./xi_C_corr_op_K_limber_FK_cut15360_permfix.npz

The superseded cut1000 sweep of the June draft is `../C_corr_op_K_limber_FK/`.
EOF
touch $NEW/PRODUCTION $RUNS/C_corr_op_O0/PRODUCTION $RUNS/C_corr_op_K_limber_FF/PRODUCTION $RUNS/C_corr_op_K_limber_FK/PRODUCTION
git rm -q $RUNS/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK_PROD_REF_2026-06-10.npz
$PY - <<'PYEOF'
from pathlib import Path
banner = """
## SUPERSEDED (2026-09-02)

This is the cut1000 sweep of the June 2026 draft, made with the corner-frozen vertex table and the
sorting callable. It is kept byte-identical because SFT-WL-B vendored it (sha256 prefix
`c3eae8cb8515668`) and because its config is the template `run_fk_variant.py` copies. The
manuscript's FK is `../C_corr_op_K_limber_FK_cut15360_permfix/`. The `+3.086e-5` figure above is a
pre-2026-06-09 stamp (before the h^4 and (1+z)^-4 fixes); this file's own kappa-kappa value at 0.5'
is `+4.2877e-6`, and the manuscript's converged value is `+1.9480e-5`. The byte-identical twin
`xi_C_corr_op_K_limber_FK_PROD_REF_2026-06-10.npz` was removed (history at tag pre-reorg-2026-09).
"""
p = Path("sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK/README.md")
if "SUPERSEDED (2026-09-02)" not in p.read_text():
    p.write_text(p.read_text().rstrip("\n") + "\n" + banner)
for name in ("C_corr_op_O0", "C_corr_op_K_limber_FF"):
    p = Path(f"sachs_sft/sftwick_outputs/2PCF/{name}/README.md")
    note = ("\n## Status (2026-09-02)\n\nLive production sweep of the manuscript (unchanged since 2026-06-01). "
            "The `+3.086e-5 baseline geometry` wording refers to the 40-point gamma grid, t_final = 2313.029 "
            "and n_gauss = 24 shared by every run; the number itself is a pre-2026-06-09 stamp of the FK "
            "kappa-kappa value. `PRODUCTION` marks this folder for the runner's overwrite guard.\n")
    if "Status (2026-09-02)" not in p.read_text():
        p.write_text(p.read_text().rstrip("\n") + "\n" + note)
print("run READMEs updated")
PYEOF

echo "=== (b2) cutoff ladder ==="
LAD=$RUNS/cutoff_ladder
mkdir -p $LAD/configs
for m in tree bihalofit; do for c in 960 1920 3840 7680 15360; do
  git mv $DFE/products/table_${m}_cut${c}_r4_xi.npz $LAD/
  git mv $DFE/products/_variant_table_${m}_cut${c}_r4/config__variant_table_${m}_cut${c}_r4.yaml $LAD/configs/config_${m}_cut${c}_r4_record.yaml
done; done
cat > $LAD/README.md <<'EOF'
# 2PCF / cutoff_ladder

The ten folded FK sweeps behind figure 17 (`fk_cutoff_convergence.pdf`) and the convergence numbers of the
cutoff subsection: `table_<model>_cut<N>_r4_xi.npz` for model = tree, bihalofit and N = 960, 1920, 3840,
7680, 15360. Each is the production FK fold (sorting callable, 40-point gamma grid, t_final = 2313.029,
n_gauss = 24, n_jobs = 4) of the vertex table `callables/kappa3_vertex/rebuild/products/table_<model>_cut<N>_r4.npz`,
assembled from the octave-band pieces by `rebuild/assemble.py` (2026-08-26; logs in `rebuild/logs/`).

`configs/` holds the verbatim run configs of 2026-08-26 (absolute paths of that machine, callable =
a materialised copy of the production callable with `TABLE_PATH` rewritten to the table). To regenerate:

    cd sachs_sft/callables/kappa3_vertex/rebuild
    sh sweep_cutoff.sh tree 960 1920 3840 7680 15360        # assemble (seconds) + fold (minutes each)
    sh sweep_cutoff.sh bihalofit 960 1920 3840 7680 15360

which writes `products/table_<model>_cut<N>_r4_xi.npz`; copy the folds here. Tree at 0.5' as a
percentage of Order-0: 0.48, 0.97, 1.53, 2.02, 2.31; BiHalofit: 0.78, 2.56, 8.01, 21.5, 41.1.
EOF

echo "=== (b3) multi-z sweeps ==="
MZ=$RUNS/multiz
mkdir -p $MZ/configs
git mv $DFE/products/order0_multiz_xi.npz $DFE/products/source_distances.npz $MZ/
for f in multiz_ff_real multiz_ff_real_all5 multiz_ff_real_z1p7 multiz_ff_real_z2p5 multiz_ff_real_z3p2 multiz_ff_real_z4p0 multiz_ff_real_z5p0ctl; do git mv $DFE/products/$f.npz $MZ/; done
git mv $DFE/products/_variant_multiz_conv/config__variant_multiz_conv.yaml $MZ/configs/config_fk_cut15360_permfix_5planes_record_2026-08-26.yaml
git mv $DFE/products/_variant_order0_multiz/config__variant_order0_multiz.yaml $MZ/configs/config_order0_7planes_record_2026-08-26.yaml
mv $DFE/products/_ff_caches $MZ/_ff_caches
cat > $MZ/README.md <<'EOF'
# 2PCF / multiz

Source-redshift sweeps behind figure 5 (`multiz_kappa_xi_cl.pdf`) and the redshift-trend notes.

* `order0_multiz_xi.npz` + `source_distances.npz`: exact Order-0 at seven source planes
  (z_s = 0.5, 1, 1.7, 2.5, 3.2, 4, 5; lambda = 1330.7 to 2318.0 Mpc), 2026-08-26, config
  `configs/config_order0_7planes_record_2026-08-26.yaml`. Note the figure's own planes (via the
  analysis3 `multiz_sweep.py` mapping) sit at lambda = 1822.721, 2094.894, 2216.793, 2266.582,
  2297.289 Mpc, so the z = 1.7 plane differs by 0.14 Mpc from this file's 2095.038.
* `multiz_ff_real_all5.npz` (and the per-plane files, plus the z = 5 control which reproduced the
  June production FF bit for bit): genuine FF sweeps at the five figure planes, one single-threaded
  process each, about 3.1 h per plane, 2026-08-27 (`analyses/revision_2026-08/run_ff_one_lambda.py`;
  logs `callables/kappa3_vertex/rebuild/logs/ff_z*.log`). `_ff_caches/` are their joblib caches (gitignored).
* The FK planes of the figure were folded with the perm-aware cut15360 vertex at the five figure planes
  (config `configs/config_fk_cut15360_permfix_5planes_record_2026-08-26.yaml`, 1034 s) directly into
  `analyses/analysis3/outputs/multiz_kappa_2pcf_5z.npz` by `analyses/revision_2026-08/regenerate_multiz.py --reuse`.
* The Order-0 planes of that npz are PCHIP-in-z interpolants of the 20-plane talk-era sweep vendored at
  `analyses/analysis3/inputs/multiz_components_talk_2026-06-10.npz`.
EOF

echo "=== (b4) runner: resolve against the YAML folder, refuse to overwrite PRODUCTION ==="
$PY $HERE/p2_patch.py b
git add $NEW $RUNS/C_corr_op_O0/PRODUCTION $RUNS/C_corr_op_K_limber_FF/PRODUCTION $RUNS/C_corr_op_K_limber_FK $LAD $MZ sachs_sft/scripts/run_FF_single.py $RUNS/C_corr_op_O0/README.md $RUNS/C_corr_op_K_limber_FF/README.md
git commit -q -m "refactor: promote the cut15360 perm-aware FK sweep, the cutoff ladder and the multi-z sweeps into sftwick_outputs run folders with relative configs; the runner resolves paths against the YAML folder and refuses to overwrite PRODUCTION results" -m "$CO"
echo "committed group b: $(git rev-parse --short HEAD)"
echo "=== status ==="; git status --short | head -10
