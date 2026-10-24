# z_s = 3.0 T-001 ladder replay (zs3p0_ladder_t001_20261023)

Scope: candidate z_s=3 finite products, no scientific acceptance.

`driver.py` replays four T-001 cutoff-ladder rungs at z_s = 3.0 with the
inputs behind the accepted z_s = 5 ladder: sft-wick 5198fb9 from the frozen
`framework_export/src`, the C23 covariance, the PCHIP loader, the
`pieces_t001` tables and GL36. Only the source endpoint changes:
`sweep.t_final_grid = [lambda_s]`, `propagators.t_max = ceil(lambda_s)` and the
contract `query_interval`. Two z_s = 5 gate folds must first reproduce the
accepted products.

| Target | Table | Guard |
|---|---|---|
| `tree_main_15360` | `pieces_t001/table_permclosed_np512_t001_refined.npz` | `main_guard` (2142 x 16), 11 ladder directions |
| `tree_auto_7680` | `pieces_t001/ladder/table_tree_auto_cut7680.npz` | `continuation_guard` (34 x 16) |
| `bihalofit_auto_7680` | `pieces_t001/ladder/table_bihalofit_auto_cut7680.npz` | `continuation_guard` |
| `bihalofit_auto_15360` | `pieces_t001/ladder/table_bihalofit_auto_cut15360.npz` | `continuation_guard` |

Optional controls `tree_auto_15360`, `bihalofit_given_7680` and
`bihalofit_given_15360` can be prepared and run in separate batches. They are
never assembled.

## Driver hash binding

The source map (`producer.driver`), every `prepared.json` (`plan.driver`),
every contract (driver pin) and every gate report (through its batch) record
the sha256 of `driver.py`. Any edit to the driver after plan step 3
invalidates the source map, all prepared batches and the gate report, and
`run`, `gate` and `assemble` then refuse. The driver must be final before
step 3 (`source-map`).

## Served endpoints

Two endpoints only. The reviewed `r1_sft061/source_map.json` serves z_s = 5.0
(the gate batch, target tag `zs5`). A source map written by this driver
serves z_s = 3.0 (production, tag `zs3p0`); its `target_tag` must equal the
tag derived from z_s. Any other z or map refuses. There is no tag override.

## Subcommands

- `source-map --z 3.0 --out FILE`: the same code must reproduce the six
  reviewed `source_map.json` endpoints bitwise, `max_affine =
  2353.929874188644`, `lambda_s = 2255.3607197066854` and `z(lambda_s) = 3.0`
  within 5e-13. The file is checked with `replay_suite.read_source_map`
  before it appears under its name.
- `prepare`: writes `config_main.yaml` (canonical JSON), `callable_main.py`,
  `main_contract.json` (canonical JSON; its sha256 is bound into the loader)
  per target, then `prepared.json`. It takes no lock and runs no fold.
  Every materialised config, contract and loader is diffed field by field
  against the accepted record of the same rung, and preparation refuses any
  difference outside paths, the endpoint (z_s != 5), the 11-direction main
  subset (`tree_main_15360`), an explicit worker count or a non-36 GL order.
  The loader is imported once in the fold environment to prove that its
  guard installs.
- `run --target ID [--gate-report FILE]`: one fold. A z_s = 3 batch first
  needs a PASS gate report (default `<driver dir>/gates_zs5_gl36/gate_report.json`),
  checked as for `assemble` before the lock is taken. Then, under `fcntl`
  `LOCK_EX|LOCK_NB` on `r1_sft061/.replay_suite.lock`: the `ps` check, a full
  rebuild of the plan from current inputs, a memory admission check and the
  loader self-check. It writes `started.json`, calls
  `run_fk_variant.run(config, result, gl, '0,2')`, validates the 11-row
  result and writes `completed.json`, or `failed.json` on any exception. An
  existing `started.json`, result or cache refuses.
- `gate`: compares the z_s = 5 gate folds with the accepted products at 1e-12
  relative and reports bitwise status. G1 uses
  `t001_accepted/cutoff_ladder/bihalofit_auto/table_bihalofit_cut15360_r4_xi.npz`.
  G2 uses `t001_accepted/products/main_fk.npz` pair (0,0) at main indices
  [0, 5..14] and the accepted tree 15360 rung. It writes `gate_report.json`
  and exits 1 on failure. An existing report refuses.
- `assemble --out DIR [--gate-report FILE]`: takes only the z_s = 3.0
  production batch. The gate report must be PASS, sit in its gate batch, equal
  a recomputation from that batch's completed folds, and share the driver,
  runtime, framework, canoes state and shared pins with the production batch.
  Every production rung that has a gate must use the gate's table, template,
  GL and worker count. It writes
  `tree_auto/table_tree_cut{7680,15360}_r4_xi.npz`,
  `bihalofit_auto/table_bihalofit_cut{7680,15360}_r4_xi.npz` and
  `assembly_manifest.json`, including a plausibility table against the
  accepted multiz z = 2.5 and 3.2 planes with no acceptance threshold. The
  tree 15360 member is the `tree_main_15360` result.

## Worker counts

The gate batch and the production batch must use the same `--n-jobs` value;
`run` and `assemble` refuse a GL36 production rung whose worker count differs
from its gate. With the plan's 8 workers:

```bash
$SFT -B $NEW/driver.py prepare --source-map $NEW/source_map_zs3p0.json --z 3.0 --gl 36 --n-jobs 8 \
  --targets tree_main_15360 tree_auto_7680 bihalofit_auto_7680 bihalofit_auto_15360 --work-dir $NEW/folds_zs3p0_gl36
$SFT -B $NEW/driver.py prepare --source-map $R1/source_map.json --z 5.0 --gl 36 --n-jobs 8 \
  --targets tree_main_15360 bihalofit_auto_15360 --work-dir $NEW/gates_zs5_gl36
```

`--n-jobs accepted` takes `sweep.n_jobs` from the accepted run of each rung
(4 for the main fold, 12 for `tree_auto_7680` and the BiHalofit auto rungs);
then pass `accepted` to both prepare commands. The accepted BiHalofit rungs
ran with 12 workers, so at 8 workers G1 also tests whether the worker split
changes values. GL24 and control batches have no gate rung to match; `run`
still requires a PASS gate report with the same driver, runtime and inputs.

The memory admission requires `vm_stat` free + inactive memory of at least the
measured T-001 peak plus 12 GiB: 37.7 GiB at 4 workers, 59.2 GiB at 8 and
75.1 GiB at 12. `--min-available-gib` replaces the rule and is recorded.
`completed.json` records `max_single_process_rss_bytes`
(`getrusage(RUSAGE_CHILDREN).ru_maxrss`), the largest single reaped process.
The T-001 memory figures above are summed process-group RSS, so the two are
not comparable.

## Failures and fresh names

Records are published complete or not at all (temporary file, fsync, then an
exclusive `renamex_np(RENAME_EXCL)`), and nothing is replaced. A failure leaves
these states:

- failed `prepare`: the `--work-dir` holds `prepare_failed.json` and no
  `prepared.json`; it never runs. Prepare again with a new `--work-dir`.
- failed or interrupted `run`: `started.json` (and `failed.json`) stay; that
  target refuses to run again. Prepare a new batch.
- FAIL `gate`: `gate_report.json` records FAIL and is never replaced. Do not
  tune tolerances or workers without an author decision.
- failed `assemble`: `--out` is not created; partial files stay in a hidden
  sibling `.<out>.<pid>.<token>.staging` directory, and a retry with the same
  `--out` works.

## Test-only flag

`prepare --allow-stage-test` relaxes one check: that the driver sits at
`r1_sft061/zs3p0_ladder_t001_20261023`. It exists so that a copy in a scratch
directory can prepare test batches. Such batches are marked
`stage_test: true`, and `run`, `gate` and `assemble` refuse them, as they
refuse any driver outside its planned directory. Never use it for real runs.

## Tests

`test_driver.py` (pytest) runs in a pytest tmp home with `run_fk_variant`
replaced by a stub, so no fold starts; the lock path, `ps` and `vm_stat` are
replaced too. It covers source map, prepare, the accepted-record diff, gate
PASS/FAIL and forged reports, the run preconditions, stubbed folds and
assembly.

```bash
$SFT -B -m pytest -p no:cacheprovider --import-mode=importlib --basetemp=<scratch dir> test_driver.py
```

## Environment

```bash
PKG=/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses
R1=$PKG/sachs_sft/analyses/r1_sft061
NEW=$R1/zs3p0_ladder_t001_20261023
FW=$R1/t001_candidate_inputs/main_fk_gl36_optimized_v6/framework_export/src
SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
export CANOES_ROOT=$HOME/projects/angular_statistics/canoes
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$FW:$CANOES_ROOT/src
```

Every subcommand refuses another interpreter, a different `PYTHONPATH` or
`CANOES_ROOT`, or bytecode writing. The fold environment is the caller's
environment without `SFT_WICK_*` and `PYTHON*` variables, plus the pinned
`PYTHONPATH`, `CANOES_ROOT`, `PYTHONDONTWRITEBYTECODE=1`, `PYTHONNOUSERSITE=1`,
`MPLCONFIGDIR` and single-thread BLAS variables. Run `run` outside the Claude
sandbox, because `ps` is blocked there.

## Limits

T-001 acceptance covers finite amplitude and sampled sign at z_s = 5. The
16-shell radial grid is concentrated near z = 5, so radial accuracy at
z_s = 3 is not established. C23 source-grid accuracy is not granted. There is
no GL24/GL36 comparison at z_s = 3 unless the optional GL24 batch is run.
Quote every FK value with its l_max. `reproduce/active_products.json` and all
historical records stay unchanged.
