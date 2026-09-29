# Sequential FK replay plan

Prepared on 2026-09-28. The user selected a true source redshift of 5 for the
updated main figures, with matching Order-0 and FF recomputation. No fold was
started while preparing this plan.

The current prepared batch retains the historical K inputs for traceability.
The separate September input-table correction is awaiting an explicit input
decision. These held-fixed historical inputs are not established as the final
scientific inputs for deployment.

## Framework and source convention

- sft-wick 0.6.1, clean commit
  `6945815ad9e6efa2b9e768891d6bf92b1214670d`.
- `framework_provenance.json` records the independently verified framework and
  canoes loader state. The runner checks both before execution.
- `source_map.json` supplies source endpoints and the explicit
  `aligned_response.py:response_R` module. It matches the coordinate mapping
  used by the saved C table.
- For the main source, lambda changes from `2313.0288751857356` to
  `2317.9695958203943` Mpc. The former is the historical endpoint, whose old
  z = 5 label is inaccurate.
- The five multi-z endpoints also change slightly. They are all read from the
  reviewed source map, without recomputing a different z-to-lambda conversion.
- True-redshift runs retain the one-Mpc propagator grid spacing and set its
  upper extent to the ceiling of the largest requested endpoint. Thus the
  response is never tabulated to the historical 2360 Mpc endpoint beyond the
  matched mapping's finite affine range. The main and ladder grids end at
  2318 Mpc, and the multi-z grid ends at 2298 Mpc.

`replay_manifest.json` is the complete plan with absolute input paths, SHA256
hashes, old and new endpoints, exact callable identities and figure dependencies.
Every required input was present when the manifest was generated. It does not
certify the scientific convergence of the input tables.

## Historical provenance and target list

| Target | Exact K input | Callable | Historical source endpoints | Reproduction figure IDs |
|---|---|---|---|---|
| `main` | `equal_time_limber_cut15360_permaware/table_permclosed.npz` | permutation aware | production `config_L2.yaml` | 2, 3, 4, 6 |
| `multiz` | the same production table | permutation aware | `analysis3/outputs/multiz_kappa_2pcf_5z.npz`, `lam` array | 5 |
| `tree_960` to `tree_15360` | each `rebuild/products/table_tree_cut<N>_r4.npz` | historical sorting callable | the corresponding recorded ladder YAML | 17 |
| `bihalofit_960` to `bihalofit_15360` | each `rebuild/products/table_bihalofit_cut<N>_r4.npz` | historical sorting callable | the corresponding recorded ladder YAML | 17 |

The ladder cutoffs are 960, 1920, 3840, 7680 and 15360. In total there are 12
sequential fold targets. Each retains its recorded 24-point Gauss rule and
40-angle geometry. The main target retains all six historical observable
component pairs. Multi-z and all ten ladder targets request only the convergence
pair `(a,b) = (0,0)`, which is the only observable they contribute to the figures
and quoted claims. The stochastic field still has three components, so this
selection changes only requested output pairs. Internal component sums and
input tensors are unchanged. `component_pairs`, `historical_component_pairs`
and `n_components` record the choice explicitly in the manifest.

Four workers are the default, with six available explicitly.
BLAS and OpenMP thread counts are limited to one within each worker process.

The main and ladder inputs are not interchangeable. The main production table
contains the permutation-closed geometry and extra rows used with its
permutation-aware callable. The ladder was originally made with its own r4
tables and the sorting callable. Keeping this distinction makes the framework
change traceable. A separate callable/table comparison is required before
replacing the historical ladder inputs with a new common input model.

The multi-z recorded YAML still contains the single baseline source endpoint.
`regenerate_multiz.py` overrode it at runtime through
`compute_multiz_kappa_2pcf.py`. Therefore the actual deployed NPZ supplies the
historical endpoint array. The old compute script's default FK configuration
points to the superseded cut1000 run and is deliberately not called.

## Separate table correction

The September `pieces_nphi512/table_permclosed_np512.npz` rebuild includes a
pair-phase correction and different angular quadrature. None of the targets
silently adopts it. The first framework replay preserves the historical table,
and a subsequent controlled table comparison can establish whether that input
correction should also enter the revised paper.

## Commands

Run from the reproduction-package root with its sft-wick interpreter. Planning
and preparation do not run any folds.

Inspect a historical-input plan:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/replay_suite.py plan \
  --source-mode historical --targets main multiz ladder
```

Prepare all true-redshift targets in a new directory:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/replay_suite.py prepare \
  --source-mode true-redshift \
  --source-map sachs_sft/analyses/r1_sft061/source_map.json \
  --targets all --n-jobs 4 \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24
```

After the current controlled fold has finished and code review is complete,
execute that prepared batch:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -u \
  sachs_sft/analyses/r1_sft061/replay_suite.py run \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24
```

The component-selection update has already prepared the distinct batch
`true_redshift_fk_gl24_historical_inputs` with four workers. Its exact execution
command, after the input choice is resolved, is:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -u \
  sachs_sft/analyses/r1_sft061/replay_suite.py run \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24_historical_inputs
```

This batch preserves the old K tables. Adoption of the corrected tables needs
an explicitly named new input mode and a fresh prepared batch, with complete
provenance for each main, multi-z and cutoff input.

`--targets` also accepts individual names, such as `main`, `multiz`,
`tree_960` or `bihalofit_15360`. To preserve an interrupted target, prepare a
fresh batch with only that target and give it a new directory name. A completed
target is skipped only after its result hash is checked. Partial runs are never
overwritten. Do not start a second fold while this batch is active.

## Output safety

- Reuses `run_fk_variant.materialise_variant` and `run_fk_variant.run`, which
  invoke the established `run_FF_single.py` numerical workflow.
- Each batch must be a new directory. Each target gets private expansion and
  propagator caches, a rewritten K callable and a config.
- Every generated output/cache path is validated to remain inside the batch.
- The verbatim historical YAML is copied to `historical_config.yaml`. The new
  YAML and callable hashes enter `prepared.json`.
- The run checks all input hashes, framework version and interpreter before
  starting. `started.json` and `completed.json` document execution.
- A shared advisory lock serializes suite batches. An additional process check
  rejects overlap with the established fold drivers, which do not use the lock.
- The runner neither replaces production data nor deploys figures.

## Required companion products and figure integration

1. The true-z = 5 main result must be paired with freshly computed Order-0 and
   FF at the same endpoint and response convention. Existing historical arrays
   cannot provide the denominator or Full = O0 + FF + FK at the new source.
2. The multi-z result needs matching Order-0 and FF at all five updated
   endpoints. The `ff` field inside `multiz_kappa_2pcf_5z.npz` is a historical
   placeholder. The old plot actually loads `multiz_ff_real_all5.npz`, whose
   source endpoints now also need alignment.
3. Reproduction IDs 2, 3 and 4 need the new main FK, O0 and FF products. Their
   angular transforms and E/B decomposition must be recomputed from those
   arrays rather than rescaling old spectra.
4. Reproduction ID 5 needs a new combined multi-z file with matched components.
   Its plot must not silently reload the old fixed-path FF file.
5. Reproduction ID 17 needs all ten new ladder folds and the new source-matched
   Order-0 denominator. Its current generator hard-codes the historical Order-0
   denominator and must be adapted before publication.
   The multi-z and ladder products now contain only the convergence pair. The
   existing `load_observables` helper expects shear pairs as well. Their new
   adapters must extract `(0,0)` directly and must not invent zero-valued shear
   fields to satisfy that helper.
6. Reproduction ID 6 includes Monte Carlo markers. Their source plane and
   driving inputs must match the new theoretical curves before interpreting
   agreement or disagreement. The framework update alone does not update those
   markers.
7. Update derived numerical claims only after these matched products are
   verified. Preserve previous figure artifacts, then use `reproduce/deploy.py`
   to copy reviewed figures into the manuscript.

## Preparation checks completed

Both historical and true-redshift modes prepared all 12 targets successfully
in separate temporary directories. Checks confirmed the intended endpoint
arrays, original K-table/callable distinction, explicit matched response in
true-redshift mode, bounded workers, isolated write paths and unstarted outputs.
Existing-directory preparation and a missing true-redshift source map were
rejected. These are preparation checks only. No new suite fold or FF job was
started and no scientific result is claimed from them.

The updated component-selection preparation was also checked across all 12
targets. It retains all 40 angles and three internal components, keeps six
requested pairs for the main target, and requests only `(0,0)` for the five
multi-z source planes and ten ladder targets. The previous manifest is preserved
as a timestamped `replay_manifest_archive_*.json` file. Previously prepared
products and configs are unchanged.
