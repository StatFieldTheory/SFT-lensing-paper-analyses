# Matched FF companion plan

Prepared on 2026-09-28. This plan creates new source-matched FF products for the
latest-framework FK figures. It does not start folds, modify production arrays
or change the covariance input table.

## Scope and checkpoints

`companion_suite.py` prepares one independent target per source plane. Each
target has a private config, expansion cache, propagator cache and result file.
The reviewed `replay_suite.run_batch` supplies the shared fold lock, active-driver
check, dependency hashes and completed-result checkpoints. FK and companion
batches therefore use the same exclusion mechanism. Completed targets survive
an interruption and are skipped only after their result hash is checked. An
interrupted target is preserved and must be repeated in a fresh directory.

| Target | Requested source | Required output pairs | Purpose |
|---|---|---|---|
| `ff_main` | z = 5 | all six historical pairs | main real-space, angular-spectrum and E/B figures |
| `ff_z1p0` | z = 1 | (0,0) | first multi-z convergence panel |
| `ff_z1p7` | z = 1.7 | (0,0) | second multi-z convergence panel |
| `ff_z2p5` | z = 2.5 | (0,0) | third multi-z convergence panel |
| `ff_z3p2` | z = 3.2 | (0,0) | fourth multi-z convergence panel |
| `ff_z4p0` | z = 4 | (0,0) | fifth multi-z convergence panel |

The five multi-z redshifts are read from the actual figure product's `z` array.
Only convergence is plotted there, so its five new FF targets omit unused
component-pair tasks while keeping the same per-pair calculation. The main
source retains all six pairs. Full targets retain all 40 historical angular
directions, including the tail needed by the angular transforms.

The six true source endpoints and response module come from `source_map.json`.
FF retains its historical propagator grid spacing of 2 Mpc. Its upper extent is
the next 2 Mpc grid point enclosing the requested source, within the matched
response module's finite affine range. The recorded local F tensor, C callable,
24-point outer Gauss rule and `vertex_types: [F]` are otherwise preserved.
The old K declaration remains in the FF configuration, with its import-time
dependency hashed, but no K diagram contributes to the F-only selection.

## Historical controls require the actual FF endpoints

The five real FF sidecars contain the endpoints 1822.721, 2094.894, 2216.793,
2266.582 and 2297.289 Mpc. These were rounded to 0.001 Mpc when the August FF
runs were dispatched. They differ slightly from the full-precision endpoints
stored in the historical FK figure product. The runner reads each sidecar's
`lam` field and verifies that its FF array is exactly the corresponding slice
of `multiz_ff_real_all5.npz`.

The `historical` mode retains those actual endpoints and `D_callable.py` for a
controlled framework comparison. The `true-redshift` mode uses the reviewed
new endpoints and `aligned_response.py`. Consequently a historical-control
comparison isolates the executable change, while the final true-redshift
comparison also includes the explicitly requested geometry alignment.

## Why FF needs a fresh check across the framework update

The installed clean source is sft-wick 0.6.1 at
`6945815ad9e6efa2b9e768891d6bf92b1214670d`. The paper's reproduction pin predates
the intervening 0.4, 0.5 and 0.6 releases. The version change cannot be treated
as a K-only update merely because the first identified defect concerned FK.

The installed `CHANGELOG.md` records changes to batched off-diagonal covariance
access, propagator index checks and Gauss-Legendre integration-domain splitting.
Version 0.6.1 specifically changes the treatment of two swept external times.
The companion configuration uses `integrate_over: all`, full component C
(`diag_C: false`), an isotropic R with diagonal response bookkeeping, and a
supplied vectorized C callable. This makes a controlled FF replay appropriate.
It does not establish that every listed fix changes this particular FF result.

Several release changes concern different paths. The production configuration
uses `c_closed_form_only: true`, so the package's internal covariance table
rebuilds do not replace the supplied C table. Matrix-R-only fixes concern a
different response configuration. The supplied C callable does not declare
`has_diagonal_kink`, so fixes conditional on that declaration must not be
assumed to apply. The historical-endpoint probe measures the net effect instead
of attributing a change to one release note without evidence.

## Two-angle verification before a full FF batch

The proposed probe keeps the first and last historical angular directions,
approximately 0.5 and 5000 arcmin, with all six main-source component pairs.
They sample the small-angle signal and large-separation tail with only two
angular tasks per component.

1. Run the historical main source with the current framework at GL24 and
   compare the signed values with the same rows of the archived main FF file.
   This measures the executable change at fixed inputs.
2. Run the aligned main source at GL12 and GL24 in separate directories. Compare
   signed values, including the tail, rather than dividing by near-zero entries.
   A GL12/GL24 comparison can reveal inadequate outer quadrature but is not a
   proof that GL24 is converged.
3. Inspect run time and memory from these probes. If their difference leaves a
   material numerical question, run a targeted GL48 probe before selecting the
   full-production quadrature. Do not start a six-plane high-order batch merely
   to estimate convergence.
4. The historical FF helper documented loky deadlocks and used one worker per
   source. The default here is one worker. A separate two-angle run with four
   workers can check the current release's behavior and value agreement before
   selecting parallelism for the full batch. This is a proposal, not a claim
   that the old deadlock persists or is fixed.

GL12 is refused for full figure runs. The structural result check verifies
finite values, source endpoint, selected perturbative order, complete component
pairs and the exact selected angular directions. Scientific convergence remains
a separate comparison.

## Runtime estimates and limits

The August documentation reports about 3.1 hours for a single-threaded,
40-angle, six-component FF source run. This is a historical measurement, not a
measurement of the new release. If costs per angular task remained comparable,
the two-angle all-component GL24 probe would take roughly 9 minutes, and each
convergence-only multi-z source would take about half an hour. One full main
source plus five convergence-only sources would then total roughly 5.5 to
6 hours with one worker and sequential folds.

These estimates ignore changes in diagram evaluation, domain partitioning,
cache behavior and startup costs. The current FK replay already makes it
unsafe to assume historical timing is retained. Use the measured two-angle FF
runtime to replace the estimates before dispatching the full batch. More workers
also multiply memory use, and no two folds may run simultaneously.

## Commands

Run from the reproduction-package root. The following command only prepares a
historical two-angle control:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/companion_suite.py prepare \
  --source-mode historical --targets main --profile probe \
  --n-gauss 24 --n-jobs 1 \
  --work-dir sachs_sft/analyses/r1_sft061/ff_probe_historical_gl24
```

Prepare a matched two-angle probe:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/companion_suite.py prepare \
  --source-mode true-redshift \
  --source-map sachs_sft/analyses/r1_sft061/source_map.json \
  --targets main --profile probe --n-gauss 24 --n-jobs 1 \
  --work-dir sachs_sft/analyses/r1_sft061/ff_probe_true_gl24
```

Use `--n-gauss 12` or `--n-gauss 48` with a new directory for the other probe
settings. After the current fold has finished and review is complete, the
explicit execution command for any prepared batch is:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -u \
  sachs_sft/analyses/r1_sft061/companion_suite.py run \
  --work-dir sachs_sft/analyses/r1_sft061/ff_probe_true_gl24
```

After probe assessment, prepare the full six-source FF batch:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/companion_suite.py prepare \
  --source-mode true-redshift \
  --source-map sachs_sft/analyses/r1_sft061/source_map.json \
  --targets all --profile full --n-gauss 24 --n-jobs 1 \
  --work-dir sachs_sft/analyses/r1_sft061/ff_true_all_gl24
```

Individual source selectors are `main`, `z1p0`, `z1p7`, `z2p5`, `z3p2` and
`z4p0`. `multiz` selects only the five lower-redshift sources. The run resumes
completed planes by hash. To repeat a failed plane, prepare only that selector
in a fresh directory.

## Optional Order-0 diagnostics and covariance limitation

`--channels o0-diagnostic` prepares separately labelled Order-0 snapshots from
the ordinary production O0 configuration, using the same source map. They are
diagnostics of the held-fixed 21-source-node C table and are not a replacement
for the independent, matched Order-0 reference being established by the
comparison task. The default channel is FF, so requesting the companion batch
does not silently dispatch redundant Order-0 work.

The C table remains unchanged for every companion run. Aligning source
coordinates and increasing outer Gauss nodes do not establish convergence of
that 21-node input interpolation. A claim of complete numerical convergence
would require its own input-table study. No covariance-table correction is
silently mixed into this framework and source-alignment update.

## Preparation checks completed

Both historical two-angle and true-redshift full-grid modes prepared all six
FF sources plus six optional O0 diagnostic sources in temporary directories.
Checks confirmed one source per checkpoint, matched response selection,
preserved FF time spacing, correct F-only selection, six main-source pairs,
one multi-z pair, isolated outputs and unstarted results. GL12 full-run
preparation was correctly refused. No fold or Monte Carlo job was started.
