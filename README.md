# `SFT-lensing-paper-analyses/` — code & data behind the STF_lensing paper

Path-integral weak-lensing pipeline. This folder holds every generator, callable, and
data product reachable from the finalized paper, plus a dated archive of superseded
material. The manuscript and final figure PDFs live at the **repo root** (`main.tex`,
`sections/`, `figures/`); this folder is git-ignored (repo-root `.gitignore`).

Start here:
- **`figure_manifest.md`** — the live dependency map: each of the 16 paper figures →
  generator, inputs, env, copy step, one-line regen command. Read this first.
- **`CLEANUP_REPORT.md`** — what the 2026-06-17 reorganization moved/deleted (before/after).

## Repo map

```
SFT-lensing-paper-analyses/
  figure_manifest.md     # live figure -> generator/inputs/env/command map
  CLEANUP_REPORT.md       # 2026-06-17 reorg report + AMBIGUOUS list
  figures/                # Python generators for schematic/glyph/operator figures
                          #   (NOT the final PDFs; those are at repo-root figures/)
  sachs_sft/              # production sft-wick pipeline (since 2026-05-30)
    callables/            # the 2 live callables that drive sft-wick:
      C_propagator/corr_op/          #   correlation propagator C = <Phi Phi>
      kappa3_vertex/equal_time_limber/ #  non-local kappa3 vertex K (+ zeta-slice figure)
    analyses/             # one folder per analysis
      analysis1/          #   -> analysis1_O0_vs_pyccl.pdf
      analysis3/          #   -> analysis3_NLO_FFFK / analysis3_cl_full_vs_O0 / cl_EB_polarization / multiz_kappa_xi_cl
      mc_sachs_2pt/       #   -> appendix_mc_workflow.pdf (Monte-Carlo Sachs 2pt cross-check)
      shear3pcf_fastnc/   #   BACKUP (removed from paper; see BACKUP_NOT_IN_PAPER.md)
    scripts/              # shared production code (D_callable, kappa2_callable, c_fast_gl,
                          #   spin_rotation) + gen_sftwick_configs + inputs/
    sftwick_outputs/2PCF/ # the 3 live sweep outputs (O0, K_limber_FF, K_limber_FK) + caches
  mathematica/            # symbolic proofs (.wl) backing the paper equations (CLAUDE.md SoT)
  docs/                   # design records
  _archive/
    canoes_pipeline/                # old tangled pipeline (archived 2026-05-30, 161M)
    cleanup_2026-06-17/             # superseded material archived in this reorg (see its README)
  sachs_sft/_cleanup_archive_2026-06-09/   # earlier h^4-fix archive
```

## Environments

Subagent shells cannot `conda activate`; use the absolute interpreter:
- **PyCCL** (analysis/plotting, scipy, camb): `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`
- **sft-wick** (L2 framework / 2PCF sweeps): `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
- LaTeX: `latexmk -pdf main.tex` (from repo root)
- Mathematica: `wolframscript -file SFT-lensing-paper-analyses/mathematica/<name>.wl`

## Regenerating the paper figures

See `figure_manifest.md` for the full per-figure table and commands. Quick notes:
- The `SFT-lensing-paper-analyses/figures/` schematic/glyph generators and the `equal_time_limber/plot_zeta_figures.py`
  zeta-slice generator **write straight into the repo-root `figures/`** — rerunning them
  overwrites the paper PDFs. Redirect (`--out-stem`) or work on a copy if you only want to test.
- `analysis*` and `mc_sachs_2pt` generators write to their own `outputs/`/`figures/` dir under
  a slightly different name, which was then **copied** into `figures/` (the paper renames drop
  suffixes like `_fkem`, `_nlo_O0_FF_FK`, `_2x5`). `fig_xi_channels.py --replot` re-plots from
  the cached `outputs/appendix_mc_curve.npz` without redoing the Monte Carlo.
- `multiz_kappa_xi_cl.pdf` is 2-stage (sft-wick compute → PyCCL plot); the expensive
  intermediate `analysis3/outputs/multiz_kappa_2pcf_5z.npz` is kept so the plot can be redone cheaply.
- 2 figures have **no in-repo generator** (`FeynDiag_2pt_order2.pdf`, `FeynDiag_3pt_order1.pdf`,
  likely hand-made vector art) and `cosmology_coord_maps.pdf`'s generator lives only in
  `_archive/canoes_pipeline/` and needs the external `~/projects/canoes/src`. The PDFs are present.

## Conventions
English in code/comments/commits. Archive, do not delete, old artifacts. Do not edit
`sections/*.tex`, `main.tex`, or `figures/*.pdf` — consume them as ground truth.
