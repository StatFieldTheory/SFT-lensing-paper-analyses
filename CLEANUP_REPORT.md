# `SFT-lensing-paper-analyses/` cleanup report — 2026-06-17

Reorganized `SFT-lensing-paper-analyses/` into a documented, reproducible tree separating live production
code from archive, without breaking the paper. Executed per `SFT-lensing-paper-analyses/CLEANUP_HANDOFF.md`
with user sign-off (2026-06-17): archive `forecast/`, full archive sweep.

## Folder rename (2026-06-17, same session)
The top-level folder was renamed **`scripts/` → `SFT-lensing-paper-analyses/`**. Updated:
the 3 live `config_L2.yaml` absolute paths, `figures/make_response_corr_operator_figures.py`
(`ROOT/"scripts"` → `ROOT/"SFT-lensing-paper-analyses"`), the `equal_time_limber` build
scripts, `gen_sftwick_configs.py`, the shear3pcf backup scripts, the `talk/scripts/*`
generators that import from this folder, the repo-root `.gitignore` entry, `CLAUDE.md` folder
tokens, these docs, and the project memory. `__file__`-relative generators and the INNER
`sachs_sft/scripts/` dir were unaffected. Frozen `_archive/` content and historical `.log`
files keep their original `scripts/...` paths (provenance). Verified: build green, 3 live
config paths resolve, `make_response_corr_operator_figures.py` runs (then figures/ restored
byte-identical), zero stale top-level `scripts/` references in live/tracked files.

## Headline results
- **Live `sachs_sft/` shrank ~900M → 86M**; live `callables/` 820M → 8.7M.
- Archived (reversible `mv`) **815M** into `SFT-lensing-paper-analyses/_archive/cleanup_2026-06-17/`.
- Deleted only regenerable caches: 32 `__pycache__`, 11 `.ruff_cache`, 2 `.pytest_cache`,
  117 `*.pyc`, 5 `.DS_Store`.
- **`latexmk -g -pdf -halt-on-error main.tex` → exit 0, 30 pages, 0 missing-figure/LaTeX errors.**
- All **16 live paper figures** present and md5-identical before/after (figures/ untouched).

## Key correction to the handoff
The paper renders **16 figures, not 17**. `driving_field_cumulants.pdf` is in a
**commented-out** block (`sections/cosmology.tex:775`) and the file does not exist; the
handoff's `grep includegraphics` over-counted it. The comment-stripped live set is 16.

## What was deleted (regenerable caches only)
`__pycache__/`, `.ruff_cache/`, `.pytest_cache/`, `*.pyc`, `.DS_Store` throughout `SFT-lensing-paper-analyses/`.
Nothing else was deleted.

## What was archived → `SFT-lensing-paper-analyses/_archive/cleanup_2026-06-17/`
(Full catalogue in that folder's `README.md`. All reversible via `mv` back.)

| Bucket | Items | Size |
|---|---|---|
| `forecast/` | whole `SFT-lensing-paper-analyses/forecast/` (no `discussion.tex`, no `.tex` reference) | ~1.7M |
| `callables_superseded/C_propagator/limber/` | superseded Limber C-table (not in any live config) | **773M** |
| `callables_superseded/kappa3_vertex/` | `R_contracted`, `R_contracted_lcap48`, `R_contracted_lcap48_fulllimber`, `R_contracted_prod_fullell`, `equal_time_exact3j`, `equal_time_limber_fulllimber` | ~37M |
| `sftwick_outputs_superseded/2PCF/` | 8 non-live sweep dirs (contracted*, exact3j, prod_fullell, FK_fulllimber, C_limber_O0) | ~1M |
| `sachs_sft_session_scratch/` | `BUG2_MATH_FLOW.*`, `HANDOVER_*.md`, 3× `PROMPT_*.md`, 3× `_draft_*` | small |
| `shared_scripts_scratch/` | `analyze_mirror_*`, `build_mirror_*_resid.py` | small |
| `analysis_scratch/analysis1/` | `plot_analysis1_kk_vs_pyccl.py`, extended-ref npz, kk outputs | small |
| `analysis_scratch/analysis3/` | `probe_fk_*.py`+pngs, `export_kappa_z1_z1p7*`, `_diag_*`, `_compute_multiz.log` | small |
| `analysis_scratch/mc_sachs_2pt/` | `compare_to_sft.py`, `_fk_converge*.py`, `fk_mc_vs_analytic.py`, `ff_channel_check.py`, `fig_fk_matrix.py`, `PROMPT_NEXT_SESSION.md`, `figures_scratch/`, `outputs_scratch/` | small |

## What was kept in place (live production)
- All 13 in-repo figure generators + `_plot_style.py` in each folder.
- Live callables: `C_propagator/corr_op/`, `kappa3_vertex/equal_time_limber/`.
- Live sft-wick outputs: `2PCF/{C_corr_op_O0, C_corr_op_K_limber_FF, C_corr_op_K_limber_FK}`.
- Shared: `scripts/{D_callable, kappa2_callable, c_fast_gl, spin_rotation, gen_sftwick_configs, run_FF_single}.py`, `scripts/inputs/`.
- Live intermediates kept so figures can be re-plotted cheaply: `analysis3/outputs/multiz_kappa_2pcf_5z.npz`, `mc_sachs_2pt/outputs/appendix_mc_curve.npz`, each analysis's production `outputs/`/`figures/` PDFs.
- `mathematica/` (proofs), `docs/` (design records).
- `sachs_sft/analyses/shear3pcf_fastnc/` — **deliberate backup** (now labeled `BACKUP_NOT_IN_PAPER.md`).
- Pre-existing archives untouched: `_archive/canoes_pipeline/` (161M), `sachs_sft/_cleanup_archive_2026-06-09/` (71M).

## New documentation
- `SFT-lensing-paper-analyses/README.md` — repo map, envs, how to regenerate figures.
- `SFT-lensing-paper-analyses/figure_manifest.md` — 16-figure dependency table + regen commands.
- `SFT-lensing-paper-analyses/.gitignore` — caches + sft-wick `.cache_*` + large regenerable tables (note: repo-root `.gitignore` already ignores all of `SFT-lensing-paper-analyses/`).
- `SFT-lensing-paper-analyses/_archive/cleanup_2026-06-17/README.md` — archive catalogue.

## Verification performed (all PASS)
1. `latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex` → exit 0, 30 pages, no missing-figure/LaTeX errors.
2. All 16 live figure PDFs present in `figures/`; md5 identical before/after the spot-regen (figures/ untouched).
3. Spot-regenerated 3 figures with the PyCCL env: `make_cumulant_hierarchy.py`→/tmp (schematic + `_plot_style`), `plot_cl_EB_polarization.py` and `plot_analysis3_cl_decomposition.py`→ own `outputs/` (exercise the live `xi_*.npz` + PyCCL path). All exit 0, sensible numbers.
4. Boundary safety-check before moving: no live generator/config references any archived path (the `driver_stats.py` "limber" hit is a docstring explaining it does NOT use that table; `gen_sftwick_configs.py` names variants by design).
5. Every absolute path in the 3 live `config_L2.yaml` resolves to an existing file post-move; archived callables confirmed gone from the live tree (no leftovers).
6. Every archived item catalogued in `_archive/cleanup_2026-06-17/README.md`.

## AMBIGUOUS / FYI items (no destructive action taken — for your awareness)
1. **`FeynDiag_2pt_order2.pdf`, `FeynDiag_3pt_order1.pdf`** — no generator anywhere in the repo (likely hand-made vector art, or a deleted generator; `talk/scripts/make_feyndiag_dark.py` makes differently-named talk SVGs). PDFs are present; the paper builds. **Not reproducible from a script** — keep the PDFs safe.
2. **`cosmology_coord_maps.pdf`** — generator lives only in `_archive/canoes_pipeline/scripts/plot_coord_maps.py` and needs the external `~/projects/canoes/src`. PDF present; reproducible only with canoes checked out.
3. **`driving_field_cumulants.pdf`** — referenced only in a commented-out block (`cosmology.tex:775`); file absent. A canoes data figure dropped from the paper. No action.
4. **`figures/cumulant_hierarchy.pdf`** + `SFT-lensing-paper-analyses/figures/make_cumulant_hierarchy.py` — generator and PDF exist but **no `.tex` references the PDF** (orphan/alternative schematic). Kept the generator; left the orphan PDF in `figures/` (I do not edit `figures/`).
5. **`gen_sftwick_configs.py`** kept in the live tree still names the archived variant callables — regenerating the *variant* configs would need those folders restored from the archive first (the 3 live configs already exist as files, so the paper is unaffected).
