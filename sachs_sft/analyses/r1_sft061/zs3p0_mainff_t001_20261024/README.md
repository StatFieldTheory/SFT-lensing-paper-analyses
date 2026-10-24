# z_s = 3.0 Order-0, FF and main-FK replay (zs3p0_mainff_t001_20261024)

Scope: candidate z_s=3 finite products, no scientific acceptance.

`driver.py` is version 2 of `../zs3p0_ladder_t001_20261023/driver.py`. It keeps
v1's guards (fresh directories, atomic exclusive publish, the shared `fcntl`
lock, the `ps` check, pin re-hashing, framework and canoes binding, contract
binding, driver-hash binding, typed refusals, memory admission) and reuses v1's
z_s = 3.0 source map and its `tree_main_15360` fold. v1 products are read and
pinned, never modified. It adds three pipelines for the NHFP figure:

| Pipeline | Gate stage | Production stage | Runner |
|---|---|---|---|
| Order-0 | `order0_parity`: z = 3.2 recorded argv, z = 5 full channels | `order0_production`: z = 3.0 production, five-multipole probe, full-channel source | `source_folded_reference.py`, PyCCL interpreter |
| FF (C23, Gaussian driving fields) | `ff_canary`: z_s = 5, main directions 0 and 6 | `ff_production`: z_s = 3.0, 40 directions | `run_fk_variant.run` -> `run_FF_single.py` |
| Main FK (tree, l_max = 15360) | `fk_regression`: z = 3.2, main directions 0, 5, 10, 15, 25, 35 | `fk_production`: z_s = 3.0, 40 directions | `run_fk_variant.run` -> `run_FF_single.py` |

All folds use pair (0,0), orders `0,2`, the frozen sft-wick 5198fb9 export,
canoes 53850f3 (clean), `t001_candidate_inputs/corr_op_C_callable.py` with the
C23 handoff table `d487c51b`, `aligned_response.response_R`, `F_tensor.npy` and
the zero-kappa2 seed.

## Served endpoints

Three only. The v1 map `../zs3p0_ladder_t001_20261023/source_map_zs3p0.json`
(pinned, written by the pinned v1 driver) serves z_s = 3.0, `lambda_s =
2255.3607197066854` Mpc. The reviewed `../source_map.json` serves the z = 3.2
regression (`2266.590169124188`) and the z_s = 5 canary (`2317.9695958203943`).
Every endpoint is re-derived with `aligned_response` and must match bitwise.

## Configurations and accepted records

Each fold config is a copy of an accepted config. Only these fields may change,
and `prepare` refuses any other difference:

| Target | Base config | Accepted records diffed field by field | Allowed differences |
|---|---|---|---|
| `fk_regress_zs3p2_gl36` | main v6 `config_main.yaml` | main fold; multiz v10 admission | owned paths, endpoint, 6-direction subset, pair (0,0); against multiz also workers (12 there) and guard scope |
| `fk_main_zs3p0_gl36` | main v6 `config_main.yaml` | main fold | owned paths, endpoint, pair (0,0) |
| `ff_canary_zs5_gl24` | C23 main release (`ff_maincap66_full_main_preparation_20261005`) | release; final_main_v11 FF admission (identical) | owned paths, C module, 2-direction subset, pair (0,0) |
| `ff_zs3p0_gl24` | C23 z3p2 release (`ff_singlepair24_selected_releases_20261005/full_ff_z3p2_gl24`) | z3p2 release; main release | owned paths, C module, endpoint (and pair (0,0) against main) |

The endpoint fields are `sweep.t_final_grid = [lambda_s]` and
`propagators.t_max = dt * ceil(lambda_s / dt)`: 2256 for z_s = 3.0 (FK dt = 1,
FF dt = 2), 2267 for the FK regression, 2318 for the FF canary. With
`c_closed_form_only`, `t_max` enters only the cache key and the R probes. The
z_s = 3.0 FK config equals v1's `tree_main_15360` config apart from owned
paths and the 40-direction grid.

The FK loader is the accepted construction: the PCHIP callable with
`TABLE_PATH` set to `pieces_t001/table_permclosed_np512_t001_refined.npz` and a
`main_guard` install bound to the canonical contract's path and sha256. It must
equal the accepted main loader after contract binding, and its body must equal
the multiz loader body (which used `continuation_guard`). The contract has the
accepted main keys; its pins must agree with every common pin of the accepted
main and multiz contracts.

The accepted main FF ran GL24 with `sweep.n_jobs = 4` in its config, executed
by the C23 row-partition runner with 24 workers (each shard `n_jobs = 1`). That
runner's plan and closure are deleted. This driver runs one `run_FF_single.py`
process with the config's 4 workers; the C23 records show 4-, 12- and 24-worker
partitions bitwise equal, and the canary tests it here. The C23 runner built C
with `make_C_fn_lambda_project(table, 'ridge', output_h_power=0)` on the handoff
table `d487c51b`; `corr_op_C_callable.py` builds the same function on the same
table, which is why it replaces the deleted `corr_op_C_candidate.py`.

The main FK runs with the accepted main fold's 4 workers (recorded as
`accepted:main_fk_gl36`). `--n-jobs N` at `prepare` overrides it, is recorded as
`cli`, and must then be used for both the gate and the production batch.

## Order-0 jobs

Every argv is the recorded z = 3.2 multiz argv of
`analysis1/r1_aligned/source_folded_multiz_plan.json` with explicit changes.
Before a run the argv is checked against the recorded settings of an accepted
product; after a run the output's recorded settings must equal them apart from
`out`, `z_source` and `source_map`.

| Job | Changes to the recorded argv | Settings reference |
|---|---|---|
| `o0_parity_zs3p2` | none (only `--out`) | `source_folded_multiz_z3p2.json` |
| `o0_parity_full_zs5` | `--z-source 5 --source-map <reviewed>`, n_chi 4096, k 2048, block 4, channels `cl_pp,cl_pw,cl_ww` | `source_folded_full_boundary_n4096_k2048.json` |
| `o0_prod_zs3p0` | `--z-source 3 --source-map <v1 map>` | `source_folded_multiz_z3p2.json` |
| `o0_probe5_zs3p0` | z and map as above, ells `2,10,100,500,800`, n_chi 4096, k 2048, block 4 | `source_folded_multiz_z3p2_probe.json` |
| `o0_full_zs3p0` | z and map as above, n_chi 4096, k 2048, block 4, channels `cl_pp,cl_pw,cl_ww` | `source_folded_full_boundary_n4096_k2048.json` |

The two full-channel jobs go beyond `order0_plan.txt`. The NHFP `[order0]`
product holds six pairs built from `kk`, `xip`, `xim` and `kgt`
(`prepare_plot_inputs.order0`), and `nhfp_fig_evidence.window_check` reads
`xi_plus` from it. A `cl_pp`-only run gives `kk` alone. `o0_full_zs3p0` uses
the protocol of the accepted z = 5 `[order0]` source, and `o0_parity_full_zs5`
checks that protocol against that source under the current canoes.

The environment is the recorded one (`OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, `CANOES_SUPPRESS_METAL_WARNING=1`) plus
`PYTHONDONTWRITEBYTECODE=1`, without `PYTHON*` or `SFT_WICK_*` variables.

## Gates

| Gate | Criterion |
|---|---|
| `order0_parity` | theta, taper and ell grids bitwise; max abs difference over peak <= 1e-5 for `kk` (z = 3.2) and for `kk`, `xip`, `xim`, `kgt` (z = 5); the accepted z = 3.2 `kk` must equal `multiz.npz` `o0[3]` bitwise |
| `ff_canary` | `t001_accepted/products/main_ff.npz` pair (0,0) at main indices 0 and 6, relative difference <= 1e-12 |
| `fk_regression` | `t001_accepted/products/multiz.npz` `fk[3]` at main indices 0, 5, 10, 15, 25, 35, relative difference <= 1e-12 |

Reports record bitwise status, the maximum differences and the workers used.
A production `run` refuses without its pipeline's PASS report, which must sit
in its gate batch, equal a recomputation from that batch, share the driver,
runtimes, framework, canoes state and shared pins, and (folds) match the
production target in GL, workers, guard, C module and table.

After the FK production, `consistency` compares the z_s = 3.0 main FK at main
indices [0, 5..14] with v1's `tree_main_15360` result (and checks that v1's
assembled l_max = 15360 member equals it), relative difference <= 1e-12.

## Export

`export` re-verifies the three gates, the consistency report and all
production results, then publishes into a new directory (default
`export_nhfp_zs3p0`):

| File | Layout |
|---|---|
| `main_fk_zs3p0.npz` | `t001_accepted/products/main_fk.npz`: members a, b, metadata_json, order, t_final, value, x, y; x, y float64 (40, 3); 40 pair (0,0) order-2 rows in main-grid order |
| `main_ff_zs3p0.npz` | same, from the FF fold |
| `order0_zs3p0.npz` | `figure_products/order0.npz`: the `prepare_plot_inputs.order0` transformation of `o0_full_zs3p0` (240 rows, order 0, pairs 00, 01 = -kgt, 02 = 0, 11 = (xip+xim)/2, 12 = 0, 22 = (xip-xim)/2) |
| `order0_comparison_zs3p0.npz` | byte copy of the `o0_full_zs3p0` output, the analogue of `[order0_comparison]` (`ell`, `ell_taper`) |
| `export_manifest.json` | sha256 of every file and member, every source product, gate and consistency evidence, checks and diagnostics |

The transformation is checked first: applied to the accepted z = 5 source it
must reproduce `figure_products/order0.npz` bitwise. Export refuses unless the
FK and FF grids equal the accepted z_s = 5 main grids, the Order-0 layout
equals the accepted one, the NHFP separation grids of the three files are
identical and equal the z_s = 5 grid, every row has `t_final = lambda_s`, the
v1 ladder grids equal the accepted z_s = 5 ladder grids and lie on the main grid
at indices [0, 5..14], and the exported main FK equals the v1 tree l_max =
15360 rung there within 1e-12 (the `nhfp_fig_evidence.load_products` checks).
The Order-0 theta grids must equal `multiz.npz` gamma and the tapers the z = 5
taper. Three diagnostics are recorded without refusing: the fine probe against
the production spectrum (flag above 5e-3), whether `kk` lies between the
multiz z = 2.5 and 3.2 Order-0 curves, and the full-channel against the
production `kk` (flag above 1e-3 of the peak).

The NHFP generator reads products through `reproduce/active_products.json`,
which is unchanged. Read the exported files by explicit path and sha256.

## Driver hash binding and test flag

Every `prepared.json`, contract (driver pin), gate report, consistency report
and export manifest binds the sha256 of `driver.py`. Copy it unchanged; an
edit afterwards invalidates all batches. `prepare --allow-stage-test` relaxes
only the check that the driver sits at `r1_sft061/zs3p0_mainff_t001_20261024`;
such batches are marked `stage_test` and `run`, `gate`, `consistency` and
`export` refuse them. Never use it for real runs.

## Memory admission

`vm_stat` free + inactive memory must cover the recorded peak plus 12 GiB:
37.7 GiB for FK at 4 workers (main v6 sampled peak), 28 GiB for FF at 4
workers (the reviewed C23 GL24 four-worker cap, a bound rather than a
measurement), 13.1 GiB for an Order-0 job (z = 5 full-channel run). Other
worker counts need `--min-available-gib`.

## Failures and fresh names

Records are published complete or not at all and nothing is replaced. A failed
`prepare` leaves `prepare_failed.json` and no `prepared.json`. A failed or
interrupted `run` leaves `started.json` (and `failed.json`); that target
refuses to run again, so prepare a new batch with `--work-dir`. A FAIL gate or
consistency report is never replaced; do not change tolerances or workers
without an author decision. A failed `export` leaves `--out` absent and the
partial files in a hidden `.<out>.<pid>.<token>.staging` sibling.

## Commands

Run outside the Claude sandbox (`ps` is blocked there), one job at a time, with
nothing else folding. `set -o pipefail` keeps the exit status through `tee`.

```bash
PKG=/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses
R1=$PKG/sachs_sft/analyses/r1_sft061
NEW=$R1/zs3p0_mainff_t001_20261024
FW=$R1/t001_candidate_inputs/main_fk_gl36_optimized_v6/framework_export/src
SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
export CANOES_ROOT=$HOME/projects/angular_statistics/canoes
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$FW:$CANOES_ROOT/src
set -o pipefail
D="$SFT -B $NEW/driver.py"
mkdir $NEW/logs

# prepare (no fold), about 1 min in total
for s in order0_parity order0_production ff_canary ff_production fk_regression fk_production; do
  $D prepare --stage $s || break
done

# gates
$D run --work-dir $NEW/order0_parity --target o0_parity_zs3p2 2>&1 | tee $NEW/logs/o0_parity_zs3p2.log
$D run --work-dir $NEW/order0_parity --target o0_parity_full_zs5 2>&1 | tee $NEW/logs/o0_parity_full_zs5.log
$D gate --work-dir $NEW/order0_parity
$D run --work-dir $NEW/ff_canary_zs5_gl24 --target ff_canary_zs5_gl24 2>&1 | tee $NEW/logs/ff_canary.log
$D gate --work-dir $NEW/ff_canary_zs5_gl24
$D run --work-dir $NEW/fk_regression_zs3p2_gl36 --target fk_regress_zs3p2_gl36 2>&1 | tee $NEW/logs/fk_regression.log
$D gate --work-dir $NEW/fk_regression_zs3p2_gl36

# production (each refuses without its PASS gate)
for t in o0_prod_zs3p0 o0_probe5_zs3p0 o0_full_zs3p0; do
  $D run --work-dir $NEW/order0_production_zs3p0 --target $t 2>&1 | tee $NEW/logs/$t.log || break
done
$D run --work-dir $NEW/ff_production_zs3p0_gl24 --target ff_zs3p0_gl24 2>&1 | tee $NEW/logs/ff_zs3p0.log
$D run --work-dir $NEW/fk_production_zs3p0_gl36 --target fk_main_zs3p0_gl36 2>&1 | tee $NEW/logs/fk_main_zs3p0.log
$D consistency --work-dir $NEW/fk_production_zs3p0_gl36

# export
$D export
```

## Tests

`test_driver.py` runs in a pytest tmp home with `run_fk_variant` and
`launch_order0` replaced by stubs that write copies of accepted products, so no
fold or Order-0 job starts; the lock path, `ps` and `vm_stat` are replaced too.
It prepares all six stages (real loader, C-wrapper and PyCCL self-checks),
covers the endpoint, diff, gate, forgery, binding, refusal and consistency
paths, and runs a full export whose files are read with a copy of
`nhfp_fig_evidence.load_channels`.

```bash
$SFT -B -m pytest -p no:cacheprovider --import-mode=importlib --basetemp=<scratch dir> test_driver.py
```

## Limits

T-001 acceptance covers finite amplitude and sampled sign at z_s = 5. The
16-shell radial grid of the K table is concentrated near z = 5, so radial
accuracy at z_s = 3 is not established. The z_s = 3.0 endpoint lies between
the C23 nodes 2250 and 2280 Mpc; C23 source-grid accuracy is not granted. FF
and FK carry the C23 l <= 5000 tapered window. The exported Order-0 is the
windowed SFT projection; the NHFP baseline is the PyCCL all-multipole sum,
which this driver does not compute. Quote every FK value with its l_max.
`reproduce/active_products.json` and all historical records stay unchanged.
