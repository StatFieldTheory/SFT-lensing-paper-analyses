# REPRODUCE.md

How to reproduce every figure and every quoted number of *Statistical Field Theory for
Weak Gravitational Lensing* from this repository. The original reproduction record
dates from 2026-09-02. The revision section below describes the subsequent numerical
update and distinguishes its validation from that historical record.

## 2026-09-29 referee revision: product selection and validation

The revision uses the actual source plane `z_s = 5`, with
`lambda_s = 2317.9695958203943 Mpc`, sft-wick **0.6.1** at commit
`6945815ad9e6efa2b9e768891d6bf92b1214670d`, and corrected kappa3 inputs with
the pair-phase fix, `n_phi = 512`, and the permutation-aware vertex callable.
The main FK result has `ell_max = 15360`. Lower source planes are evaluated at
their corresponding affine endpoints. The eight affected reproduction IDs are
**1, 2, 3, 4, 5, 6, 7, and 17** in section 3.

`reproduce/active_products.json` selects the complete revision bundle. The figure
and number readers verify its recorded SHA256 hashes, including the members of
product directories. A selected missing or altered product raises an error rather
than falling back to historical data. Without that manifest, the generators use
their historical inputs. Activate the manifest only after all required products
and acceptance checks are complete. The candidate manifest under
`sachs_sft/analyses/r1_sft061/` is not an activation or deployment record.

Revision products, source maps, run configurations, input hashes, numerical checks,
and execution records are under `sachs_sft/analyses/r1_sft061/`. In particular,
`framework_provenance.json` records the framework and canoes versions and working
state, and `RUN_STATE.md` records completion and outstanding gates. The aligned
Order-0/PyCCL comparison is under `sachs_sft/analyses/analysis1/r1_aligned/`.
The exact tracked canoes patch matching `framework_provenance.json` is saved as
`sachs_sft/analyses/r1_sft061/canoes_worktree.patch`; apply it with `git apply`
to the recorded canoes commit in a separate checkout when reconstructing that
runtime. The external repositories are not changed by this submission.
Historical products are retained. The pre-update figure and manuscript copies
are recorded in
`sachs_sft/analyses/r1_sft061/deployment_archive_2026-09-28/manifest.json`.

The validation scope is specific:

- Order-0 uses the source-folded calculation accepted in the aligned PyCCL check.
  FF retains the original 21-source-node correlation-propagator table, with its
  evolution endpoint updated to the actual source plane. Its quadrature probes
  test the outer integration, not convergence of that underlying table.
- The FF Monte Carlo comparison uses its matched local covariance model. It does
  not independently validate the full non-local cosmological FF input. The FK
  comparison records finite-sample, source-grid and injection effects separately.
- Physical FK BB vanishes at this order for the scalar inputs and observable map
  used here. The recovered FK BB curve is a numerical residual. The finite-angle
  transform produces leakage, but this check does not isolate every source of
  numerical error. The pure-E recovery check is input-specific, not a universal error
  bound. Finite cutoff increments likewise do not establish an infinite-cutoff
  limit.

The replacement cutoff ladder and all 12 active products were accepted on
2026-09-29. The eight affected figures were regenerated, visually inspected and
deployed through `reproduce/deploy.py`. All 17 figure comparisons passed: eight
updated files are `IDENTICAL(raw)` and nine unchanged files are
`IDENTICAL(content)`, with zero raster differences. The number checker completed
with 33 entries. See `final_figure_check.txt` and `final_numbers.md` under
`sachs_sft/analyses/r1_sft061/`. The revised manuscript compiled to 33 pages,
with no undefined citations or references or multiply defined labels. Rendered
page 1 and pages 19--25 and 30 were inspected. Template/layout and bibliography
warnings remain, so this is not a warning-free build. The manuscript response
and detailed build records are in the parent project's ignored `referee-r1/`.

Regeneration and numerical checks use these commands:

```bash
python reproduce/regen_figure.py --all      # rebuild all 17 figures into their outputs/ folders
python reproduce/check_figures.py           # compare each with the deployed figures/*.pdf
python reproduce/check_numbers.py           # recompute every quantitative claim
```

> **Historical warning, recorded 2026-09-06.** The kappa3 vertex rebuilt at
> `n_phi = 512` with the `zeta_Bmod` pair-phase fix passed an acceptance test that
> the table deployed at that time failed (864/864 cells). That historical table
> was left unchanged for reproduction. Its FK/Order-0 at 0.5' changed from 2.31
> to 1.74 per cent in the comparison then performed. These are historical setup
> numbers, not the actual-z_s=5, sft-wick 0.6.1 revision values above. Read
> [HANDOFF_2026-09-06_nphi512_vertex_rebuild.md](HANDOFF_2026-09-06_nphi512_vertex_rebuild.md)
> before quoting any FK number. `check_numbers.py` verifies claim-vs-product
> agreement, not product convergence, so it passes on both tables.

In the 2026-09-02 record, `regen_figure.py --all` took about 90 seconds and
`check_figures.py` printed `IDENTICAL(raw)` for all 17 PDFs.
`reproduce/deployed_md5.txt` records the deployed-byte contract. Repeat the checks
after any revision deployment rather than treating the historical result as current.

## 1. Environments

The following environment inventory and pins are the 2026-09-02 record. The
revision's framework and input provenance is recorded separately above.

Two conda environments, named here by their absolute interpreter. Override with the
environment variables `STF_PYCCL`, `STF_SFTWICK` and `CANOES_ROOT`.

| name | interpreter | what it has |
|---|---|---|
| PyCCL | `/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python` | Python 3.11.15, numpy 2.4.2, scipy 1.17.1, matplotlib 3.10.8, pyccl 3.3.1, camb 1.6.5 |
| sft-wick | `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python` | Python 3.14.3, matplotlib 3.10.8, `sft_wick` from `~/projects/SFT/sft-wick/src` |

Exact package lists: `reorg_2026-09/env/pip_freeze_{PyCCL,sft-wick}.txt`.

**Which environment matters.** The two agree on every number but differ at the 1 pt
bounding-box level, because matplotlib's tight bounding box depends on the font metrics
each environment resolves. A figure reproduces byte-for-byte only in the environment
listed for it in section 3; `regen_figure.py` picks the right one automatically.

**External repositories.** Two are needed and both are pinned in `reorg_2026-09/env/`:

| repository | role | pin |
|---|---|---|
| `sft-wick` (`~/projects/SFT/sft-wick`) | the L2 expansion engine: every sweep, and the Feynman diagrams | commit in `sft-wick_HEAD.txt`; the diagrams were first drawn at `dd03748` and reproduce unchanged at the pinned commit |
| `canoes` (`~/projects/angular_statistics/canoes`) | the correlation-propagator table loader and the vertex build | commit in `canoes_HEAD.txt`, branch `feat/nonlinear-bdelta-vertex`, working-tree patch in `canoes_worktree.patch` |

The callables locate canoes through `CANOES_ROOT`, then an importable `canoes`, then two
known locations, so no `PYTHONPATH` is needed for the figures. The vertex *build* and the
audit scripts do need `PYTHONPATH=$CANOES_ROOT/src` and the canoes virtualenv.

## 2. Repository map

This is the historical production layout. Revision replacements are kept under
`analyses/r1_sft061/` and selected by the manifest described above.

```
SFT-lensing-paper-analyses/
  reproduce/            regen_figure.py, check_figures.py, check_numbers.py, deploy.py,
                        pdfcmp.py, deployed_md5.txt, numbers.json, logs/
  figures/              generators for the schematic and operator figures; outputs/ holds their renders
  sachs_sft/
    callables/
      C_propagator/corr_op/                    the correlation propagator C = <Phi Phi> and its table
      kappa3_vertex/
        equal_time_limber_cut15360_permaware/  PRODUCTION vertex: permutation-aware callable + table_permclosed.npz
        equal_time_limber/                     the superseded cut1000 vertex (see its DEPRECATED.md), the build
                                               script, the ell-band decomposition and the figure-7 generator
        rebuild/                               factored vertex rebuild: LOW branch, octave bands, assembly,
                                               the fold runner, probes, products/ and logs/
        b_model/                               matter-bispectrum models (tree, BiHalofit) for the build
    sftwick_outputs/2PCF/
      C_corr_op_O0/                            Order-0 sweep            (PRODUCTION)
      C_corr_op_K_limber_FF/                   FF sweep                 (PRODUCTION)
      C_corr_op_K_limber_FK_cut15360_permfix/  FK sweep of the paper    (PRODUCTION)
      C_corr_op_K_limber_FK/                   the superseded June FK sweep (PRODUCTION marker, SUPERSEDED banner)
      cutoff_ladder/                           the ten folds behind figure 17
      multiz/                                  the source-redshift sweeps behind figure 5
    analyses/           analysis1, analysis3, mc_sachs_2pt, mc_fk_complete,
                        revision_2026-08 (the drivers that made the 2026-08 figures),
                        fk_audit_2026-08 (the audit probes), shear3pcf_fastnc (backup)
    scripts/            shared code and the sweep runner
  docs/                 fk_audit_2026-08 (notes + the 22-page audit note), agent_briefs,
                        the (1+z)^4 design records; papers/ and private/ are local only
  mathematica/          symbolic proofs backing the manuscript's equations
  reorg_2026-09/        the 2026-09 reorganisation: plan, evidence, scripts, environment pins
```

## 3. The 17 figures

The command table records the historical generators. With the revision manifest
active, `regen_figure.py` routes figure 1 through the aligned comparison generator,
and the affected readers select their hash-verified revision products.

`Env` is the environment that reproduces the deployed bytes. `+canoes` means the run
imports canoes (the resolver finds it; no `PYTHONPATH` needed). Times are from
2026-09-02 on a 103 GB machine.

| # | paper figure | command (`python` = the listed environment) | env | time |
|---|---|---|---|---|
| 1 | `analysis1_O0_vs_pyccl.pdf` | `sachs_sft/analyses/analysis1/plot_analysis1_O0_vs_pyccl_fkem.py` | sft-wick | 2 s |
| 2 | `analysis3_NLO_FFFK.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_nlo_decomposition.py` | PyCCL | 2 s |
| 3 | `analysis3_cl_full_vs_O0.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_cl_decomposition.py` | PyCCL | 4 s |
| 4 | `cl_EB_polarization.pdf` | `sachs_sft/analyses/analysis3/plot_cl_EB_polarization.py` | PyCCL | 2 s |
| 5 | `multiz_kappa_xi_cl.pdf` | `sachs_sft/analyses/analysis3/plot_multiz_kappa_2x5.py` | PyCCL | 3 s |
| 6 | `appendix_mc_workflow.pdf` | `sachs_sft/analyses/mc_sachs_2pt/fig_xi_channels.py --from-cache` | PyCCL +canoes | 4 s |
| 7 | `zeta_driving_field_slices_draft.pdf` | `sachs_sft/callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py` | PyCCL | 1 s |
| 8, 9 | `corr_operator_slices_draft.pdf`, `response_operator_draft.pdf` | `figures/make_response_corr_operator_figures.py` | PyCCL +canoes | 65 s |
| 10 | `screen_basis_spin2_pattern.pdf` | `figures/make_screen_basis_spin2_figure.py` | PyCCL | 1 s |
| 11 | `driving_field_portrait.pdf` | `figures/plot_driving_field_portrait.py --seed 42` | PyCCL | 3 s |
| 12 | `selection_rule_schematic.pdf` | `figures/make_selection_rule_schematic.py --seed 42` | PyCCL | 1 s |
| 13 | `null_geodesic_congruence.pdf` | `figures/make_null_congruence_schematic.py` | PyCCL | 1 s |
| 14 | `cosmology_coord_maps.pdf` | `figures/make_cosmology_coord_maps.py` | PyCCL | 1 s |
| 15, 16 | `FeynDiag_2pt_order2.pdf`, `FeynDiag_3pt_order1.pdf` | `figures/make_feyndiag.py` | sft-wick | 1 s |
| 17 | `fk_cutoff_convergence.pdf` | `sachs_sft/callables/kappa3_vertex/rebuild/fig_cutoff_paper.py` | PyCCL | 1 s |

Each generator writes into its own `outputs/` folder and never into the manuscript.
`reproduce/deploy.py` is the only step that copies into `../figures/`; it prints the old
and new md5 and records the result in `reproduce/deployed_md5.txt`.

Expected md5 sums are in `reproduce/deployed_md5.txt`. They change whenever a figure is
re-rendered, because matplotlib stamps the creation date into the PDF;
`check_figures.py` reports `IDENTICAL(content)` in that case and `IDENTICAL(raw)` when
the bytes match exactly.

## 4. Data products the figures read

The table below records the historical products. Their run-folder READMEs give the
original provenance. The active revision selects replacements through
`reproduce/active_products.json`, rather than overwriting these arrays.

| product | figures and numbers it feeds | how it was made | cost to remake |
|---|---|---|---|
| `sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz` | 1, 2, 3, 5, 6, 17 and every ratio | sft-wick L2 sweep, Order-0, 40-point gamma grid, `t_final = 2313.029`, `n_gauss = 24` (2026-06-01) | about 2 min |
| `.../C_corr_op_K_limber_FF/xi_C_corr_op_K_limber_FF.npz` | 2, 3, 4, 5, 6 | same geometry, `vertex_types: [F]` (2026-06-01); a control re-run on 2026-08-27 reproduced it bit for bit | about 3 h |
| `.../C_corr_op_K_limber_FK_cut15360_permfix/xi_...npz` | 2, 3, 4, 6 and every FK number | the fold of the permutation-aware vertex at `ell_max = 15360` (2026-08-26); **re-run 2026-09-02, byte-identical** | 2 min 50 s |
| `.../cutoff_ladder/table_{tree,bihalofit}_cut{960..15360}_r4_xi.npz` | 17, the convergence numbers | `rebuild/sweep_cutoff.sh`, ten folds (2026-08-26) | minutes each, plus the band builds |
| `.../multiz/multiz_ff_real_all5.npz` | 5 | five single-threaded FF sweeps, one per source plane (2026-08-27) | about 3.1 h per plane |
| `.../multiz/order0_5planes_figure_lambdas_xi.npz` | a cross-check on 5 | exact Order-0 at the figure's own five planes (2026-09-02) | 85 s |
| `analyses/analysis3/outputs/multiz_kappa_2pcf_5z.npz` | 5 | Order-0 interpolated in z from the vendored talk sweep, FK folded at the five planes, FF from the file above | 1034 s for the FK planes |
| `callables/kappa3_vertex/equal_time_limber_cut15360_permaware/table_permclosed.npz` | every FK result | `rebuild/assemble.py` over the octave-band pieces (2026-08-26) | the band builds are hours; `rebuild/products/pieces/` holds them (untracked, see its note) |
| `callables/kappa3_vertex/equal_time_limber/outputs/zeta_bands_cut15360_gmax85.npz` | 7 | `ell_band_decomp.py` with eight HIGH windows to 15360 (2026-08-27) | about 5 min |
| `analyses/mc_fk_complete/_markers_pooled.npz` | the FK markers of 6 | four independent 24-seed blocks, `sigma_lambda -> 0` extrapolation (2026-08-28) | hours |
| `callables/C_propagator/corr_op/corr_op_table_*.npz` | 8, 9 and every sweep | canoes `corr_op` table build (2026-05-30) | hours |

Two inputs are deliberately not regenerable here: `rebuild/products/pieces/` (50 MB of
octave-band builds, untracked, listed with md5 in
`reorg_2026-09/dfe_manifest_2026-09-02.txt`) and the arXiv PDFs under `docs/papers/`.

## 5. The numbers

`reproduce/check_numbers.py` recomputes the quantitative statements and prints, for each,
what the manuscript says next to what the products give. Run it after changing any
product. The table as of 2026-09-02 is `reorg_2026-09/phase3/numbers_after_moves.md`;
every row in that historical record agrees with the manuscript then present. The
active revision also labels harmonic residuals and finite-cutoff sensitivity as
diagnostics rather than physical accuracy claims.

The following four discrepancies were recorded against the historical products.
They are retained for traceability and do not describe the revised source-plane
setup or establish the accuracy of the replacement products:

* the nonlinear-propagation term overtakes Order-0 at 180 arcmin, where the text says
  "beyond about 200";
* the total Order-2 correction to `C_kk` runs 1.30 to 2.47 per cent over the quoted
  multipole range, where the caption says "one to two percent";
* the extrapolated `ell_max -> infinity` value is 17 per cent above the quoted one, where
  the text says "of order ten percent";
* `t_final = 2313.029` is labelled `z_s = 5` but is `z = 4.70` on the closed-form
  background used by `D_callable.py` (which gives `lambda(5) = 2317.96`). Every run and
  every figure shares that value, so all ratios are consistent; someone computing
  `lambda(z = 5)` independently would not reproduce the absolute amplitudes to better
  than about 5 per cent.

## 6. Re-running a sweep

Sweeps write where their config says, and `scripts/run_FF_single.py` clears the caches and
replaces the output. A run folder carrying a `PRODUCTION` marker is protected: the runner
refuses to overwrite it unless `SFT_WICK_FORCE_OVERWRITE=1` is set. A variant config that
inherited a production output path destroyed the deployed FK sweep once, on 2026-08-25.

To fold an alternative vertex table without touching anything deployed:

```bash
cd sachs_sft/callables/kappa3_vertex/rebuild
$SFTW run_fk_variant.py <table>.npz \
    --callable ../equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py \
    --config ../../../sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/config_L2.yaml \
    --tag mytag --n-jobs 6 --work-dir /tmp/myvariant --out /tmp/myvariant/xi.npz
```

`--n-jobs` is a memory control, not a speed knob: sft-wick forks with loky and every
worker reloads the inputs. Never run two folds at once. Never fan out the canoes vertex
builds; that cost this machine two out-of-memory reboots in August. `rebuild/README.md`
has the chunking and queueing rules.

## 7. Building the manuscript

```bash
latexmk -g -pdf -halt-on-error main.tex
```

from the repository root. The 2026-09-02 build had 32 pages, 17 figures, and no
undefined references. Check the new build independently after revision.

## Revision repository packaging (2026-09-29)

The revision commit retains all selected arrays, comparison and convergence
probes, run configurations, source snapshots, input manifests, and validation
records. Local rollback manuscript copies, rendered previews, joblib caches,
lock files and new execution logs are ignored; they remain on disk. The
rollback archive manifest remains tracked. These ignored files are not needed
to redraw the figures from the accepted numerical products.

Recorded absolute paths are provenance for the original workstation. Replaying
all numerical runs on a different machine requires rebasing those paths and
installing the pinned external environments; the guards intentionally reject
changed input or runtime hashes. Saved `*.source.py` files are immutable run
snapshots, not additional entry points. Figure regeneration and claim checks
were verified locally, not in a clean machine installation.
