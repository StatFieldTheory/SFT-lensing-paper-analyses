# Handoff: reorganize `scripts/` into a clean, reproducible repo

You are a fresh agent. The paper at `/Users/zzhang/Documents/MyDrafts/STF_lensing`
is **finalized**. The `scripts/` folder is messy: it contains the current
production code that generates the paper's figures and data, mixed with stale,
superseded, and incorrect material. **Your job: turn `scripts/` into an organized,
documented, reproducible repo, separating live production code from archive, WITHOUT
breaking the paper.** After the reorganization, you MUST run a verification/review
phase (the user explicitly asked for a post-cleanup review).

This is a large, careful task. Plan first, classify before you move anything, move
incrementally, verify continuously, and ask the user about genuinely ambiguous
items rather than guessing.

---

## 0. Read these first (sources of truth — do not skip)

1. `/Users/zzhang/Documents/MyDrafts/STF_lensing/CLAUDE.md` — repo layout,
   conventions, env paths, authoritative formula sources, and the rule
   **"Archive, do not delete, old artifacts."**
2. Auto-memory index and entries:
   `~/.claude/projects/-Users-zzhang-Documents-MyDrafts-STF-lensing/memory/MEMORY.md`
   then read the entries relevant to the pipeline and figures, especially:
   `project_sachs_sft_pipeline_reorg`, `project_cl_fullvso0_fig_2026-06-10`,
   `project_multiz_kappa_2x5_fig_2026-06-15`, `project_cl_EB_parity_null_fig_2026-06-16`,
   `project_fig14_fastnc_validates_open_triangles_2026-06-14`, the kappa3 fix entries,
   and the figure-generation entries.
   **Caveat:** memory is point-in-time; verify every file:line / behavior claim
   against the current code before acting on it.

## 1. Hard constraints (NON-NEGOTIABLE)

- **The paper must still build.** `latexmk -pdf -interaction=nonstopmode main.tex`
  must stay green, and every figure the paper `\includegraphics` must remain present
  in `figures/` and remain reproducible from a documented generator.
- **The paper depends on exactly these 17 figures** (verify with
  `grep -rhoE 'includegraphics\[[^]]*\]\{figures/[^}]+\}' sections/*.tex main.tex`):
  `FeynDiag_2pt_order2.pdf`, `FeynDiag_3pt_order1.pdf`, `analysis1_O0_vs_pyccl.pdf`,
  `analysis3_NLO_FFFK.pdf`, `analysis3_cl_full_vs_O0.pdf`, `appendix_mc_workflow.pdf`,
  `cl_EB_polarization.pdf`, `corr_operator_slices_draft.pdf`, `cosmology_coord_maps.pdf`,
  `driving_field_cumulants.pdf`, `driving_field_portrait.pdf`, `multiz_kappa_xi_cl.pdf`,
  `null_geodesic_congruence.pdf`, `response_operator_draft.pdf`,
  `screen_basis_spin2_pattern.pdf`, `selection_rule_schematic.pdf`,
  `zeta_driving_field_slices_draft.pdf`.
  (`figures/` lives at the repo root, NOT under `scripts/`. Do not move or edit it.)
- **ARCHIVE, never delete** production-or-maybe-production material. Move stale items
  into a dated archive dir (e.g. `scripts/_archive/cleanup_2026-06-17/...`) with a
  README cataloguing what and why. The ONLY things safe to hard-delete are
  regenerable caches: `__pycache__/`, `.pytest_cache/`, `.ruff_cache/`, `*.pyc`
  (and ideally add them to a `.gitignore`).
- **Do NOT edit** `sections/*.tex`, `main.tex`, `biblio.bib`, or any `figures/*.pdf`.
  Consume them as ground truth.
- When unsure whether something is stale: **archive (reversible), don't delete**,
  and list it as AMBIGUOUS for the user.
- Keep `sections/_archive/input_validation_fastnc.tex` and the fastnc backup intact:
  the input-validation (fastnc) figure was deliberately removed from the paper but
  kept as backup; its scripts (`analyses/shear3pcf_fastnc/`) and figure
  (`figures/appendix_3cumulant_fastnc.pdf`) are NOT in the live paper tree but must
  be PRESERVED as backup, clearly labeled.

## 2. Environments (subagent Bash cannot activate conda; use absolute interpreters)

- PyCCL (analysis/plotting): `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
- sft-wick (L2 framework / 2PCF sweeps): `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
- LaTeX: `latexmk -pdf main.tex`
- Mathematica: `wolframscript -file <path>.wl`

## 3. Known current state (verify, then build on it)

- **Current pipeline:** `scripts/sachs_sft/` (since 2026-05-30). Structure:
  `callables/{C_propagator,kappa3_vertex}/<impl>/` (one self-contained folder per
  implementation), `analyses/{analysis1,analysis3,mc_sachs_2pt,shear3pcf_fastnc}/`,
  `scripts/` (shared code + `inputs/`), `sftwick_outputs/{2PCF,3PCF}/` (the produced
  correlation data + sft-wick `.cache_*` dirs).
- **Already-archived old pipelines (leave archived, just catalogue):**
  `scripts/_archive/canoes_pipeline/` (old tangled pipeline) and
  `scripts/sachs_sft/_cleanup_archive_2026-06-09/`.
- **Other top-level dirs:** `scripts/forecast/` (diagonal-Fisher forecast consumed by
  `sections/discussion.tex` — PRODUCTION), `scripts/mathematica/` (`.wl` symbolic
  proofs — PRODUCTION/reference), `scripts/figures/` (some figure generators, e.g.
  the driving-field glyph figures), `scripts/docs/`.
- **CRITICAL gotcha — copy/rename steps.** Several generators write a PDF into an
  `analyses/.../outputs/` dir under a different name, which is then copied into
  `figures/`. Your manifest MUST record both the generator output path and the copy.
  Confirmed examples (verify each):
  - `analyses/analysis3/plot_analysis3_cl_decomposition.py` → `outputs/analysis3_cl_O0_FF_FK.pdf` → copied to `figures/analysis3_cl_full_vs_O0.pdf` (PyCCL env; curved-sky Wigner-d transform).
  - `analyses/analysis3/plot_cl_EB_polarization.py` → `outputs/cl_EB_polarization.pdf` → copied to `figures/cl_EB_polarization.pdf` (PyCCL).
  - `analyses/analysis3/{compute_multiz_kappa_2pcf.py (sft-wick), plot_multiz_kappa_2x5.py (PyCCL)}` → `figures/multiz_kappa_xi_cl.pdf` (note: the 5-z 2PCF is PCHIP-interpolated in z from `talk/assets/figures/_data/multiz_components.npz`; document this external input).
  - `analyses/mc_sachs_2pt/fig_xi_channels.py` → `outputs/xi_kappa_channels.pdf` (and cache `outputs/appendix_mc_curve.npz`) → copied to `figures/appendix_mc_workflow.pdf` (PyCCL).
  - `analyses/analysis3/plot_analysis3_nlo_decomposition.py` → `figures/analysis3_NLO_FFFK.pdf`.
  - `analyses/analysis1/...` → `figures/analysis1_O0_vs_pyccl.pdf`.
  - zeta/corr slices → from `callables/kappa3_vertex/.../plot_zeta_figures.py` and the C_propagator slice generators (`zeta_driving_field_slices_draft.pdf`, `corr_operator_slices_draft.pdf`, `driving_field_cumulants.pdf`).
  - `screen_basis_spin2_pattern.pdf`, `driving_field_portrait.pdf` → `scripts/figures/` glyph generators (PyCCL; driving_field_portrait needs camb, seed 42).
  - The schematic/diagram figures (`FeynDiag_2pt_order2`, `FeynDiag_3pt_order1`,
    `null_geodesic_congruence`, `selection_rule_schematic`, `cosmology_coord_maps`,
    `response_operator_draft`) — trace these: likely TikZ standalones or Mathematica
    (`scripts/mathematica/`). Find the true source for each.

## 4. Method (phased; do not jump to moving files)

**Phase A — Inventory & dependency map (read-only).** For each of the 17 paper
figures, find its generator (`grep -rl "<figure-basename>" scripts/`), the generator's
inputs (npz/tables/configs it reads), the env, and any copy/rename step. Produce
`scripts/figure_manifest.md`: a table `paper figure | generator (path) | inputs |
env | one-line regeneration command | notes (copy step, external input)`. This map
defines the **live dependency tree** (every script + data file reachable from a paper
figure, plus `forecast/` for discussion.tex and `mathematica/` proofs).
Consider dispatching a few parallel read-only Explore/`general-purpose` subagents to
trace independent figure families (analysis1, analysis3, mc_sachs_2pt, schematics,
glyphs) in parallel, then synthesize.

**Phase B — Classify every path** under `scripts/` into exactly one of:
PRODUCTION (in the live tree), SHARED-UTIL, DATA (produced outputs/caches needed by
production), ALREADY-ARCHIVED, REGENERABLE-CACHE (safe delete), STALE/SUPERSEDED
(archive), or AMBIGUOUS (ask user). Use git history / memory / "PRE…"/"_archive"/
"_tmp"/dated names / superseded-impl folders as stale signals — but VERIFY a file is
not in the live tree before archiving it.

**Phase C — Propose target structure to the user BEFORE large moves.** Suggested:
```
scripts/
  README.md            # overview, envs, "how to regenerate every paper figure", repo map
  figure_manifest.md   # the Phase-A table
  sachs_sft/           # production pipeline (callables/, analyses/, scripts/, sftwick_outputs/)
  forecast/            # Fisher forecast (discussion.tex)
  mathematica/         # symbolic proofs (.wl) + docs
  figures/             # (clarify: what is scripts/figures vs repo-root figures/?)
  _archive/            # ALL stale/superseded, in dated subdirs, each with a README
```
Within `sachs_sft/analyses/*`, separate each analysis's production generators from
its `_tmp`/`probe_*`/`_archive_*`/superseded scratch (archive the latter). Keep the
fastnc backup clearly labeled.

**Phase D — Execute** conservatively: delete only regenerable caches; `mv` stale to
the dated `_archive/`; add a `.gitignore` for caches/large regenerable outputs.
Never `mv` anything in the live dependency tree.

## 5. Review / verification (MANDATORY — the user asked for this)

1. `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` → must be green,
   and 0 missing-figure errors. (`figures/` must be untouched.)
2. Confirm all 17 paper figures still exist in `figures/` and that
   `figure_manifest.md` lists a working generator + command for each.
3. **Spot-regenerate** the cheap PyCCL figures to a scratch dir and confirm they
   still run and produce a PDF (do NOT overwrite the paper's `figures/` copies;
   compare or write elsewhere). For expensive/sft-wick figures, verify the generator
   + all declared inputs still exist and the documented command is correct (don't
   necessarily re-run multi-hour sweeps).
4. Verify nothing in the live dependency tree was moved into `_archive/`
   (cross-check the manifest's paths against their current locations).
5. Verify every archived item is catalogued in an `_archive/.../README` with the
   reason it is stale.
6. Produce `scripts/CLEANUP_REPORT.md`: what moved where, what caches were deleted,
   the before/after tree, and the AMBIGUOUS list for the user to adjudicate.
7. Update the project memory with the new `scripts/` structure (new entry; update
   `project_sachs_sft_pipeline_reorg` if appropriate) so future sessions know the
   layout.

## 6. Deliverables
- `scripts/README.md`, `scripts/figure_manifest.md`, `scripts/CLEANUP_REPORT.md`.
- A dated `scripts/_archive/cleanup_2026-06-17/` (or similar) with a cataloguing README.
- A `.gitignore` covering caches and large regenerable outputs.
- Green `latexmk`; all 17 figures intact and reproducible; memory updated.
- A summary message listing AMBIGUOUS items you did NOT move, for the user to decide.

## 7. Working style
Plan-first with a TODO list. Be conservative and reversible. Do the inventory and
get user sign-off on the target structure before bulk moves. Prefer archiving over
deleting. End with the mandatory review. English in code/comments/commits; do not
git-commit unless the user asks.
