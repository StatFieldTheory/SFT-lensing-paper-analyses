# Factored rebuild of the kappa3 vertex

Everything the note reports comes from these scripts. The design point is
that the expensive build is done ONCE and every variant is derived from it
exactly, rather than rebuilt.

## Why a factored build is exact

Three properties, all tested rather than assumed:

* **Rows are independent.** `zeta[i]` depends only on `cosine_triples[i]`.
  Computing a row in a six-row table and in a two-row table gives identical
  doubles, for both the LOW Wigner-3j branch and the HIGH Limber branch
  (`np.array_equal` true). So one build at the finest triple set yields
  every coarser subset by row selection, and chunking a build over rows is
  exact.
* **LOW is independent of the HIGH cutoff and of the bispectrum model.** It
  is always tree level, because its FFTlog machinery needs a separable
  bispectrum. One LOW build serves every variant on the same grid.
* **The HIGH octave windows are disjoint in maximum leg.** Their sum is the
  same integral, so partial sums give every intermediate cutoff.

Validation: assembled on the deployed 1671 rows with the deployed settings,
the rebuild reproduces the deployed table to a median relative difference
of 3e-16 across all six channels.

## Memory

`compute_kappa3_sigma3_high` holds a working array proportional to
(rows x n_ell x n_ell x n_phi), which is about 8 GB at 1828 rows with the
default quadrature. Eighteen of those at once cost a 103 GB machine an
out-of-memory reboot on 2026-08-25. Two guards:

* `build_band.py --chunk-rows 192` bounds the per-process peak to about
  3.2 GB;
* `run_queue.py` caps concurrency and refuses to start a job when free
  memory is below a floor, re-checked before every launch because other
  sessions share the machine.

Folds are a separate risk: sft-wick forks with loky and every worker
reloads the inputs, so `run_fk_variant.py --n-jobs` is a memory control,
not a speed knob. Never run two folds at once.

## Pipeline

```sh
CAN=/Users/zzhang/projects/angular_statistics/canoes
P=../../products
export PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1

# 1. triple sets: production mesh plus a densified collapsed family
python ../make_dense_triples.py --refine 4 --out $P/triples_full_r4.npz

# 2. the LOW branch, once
$CAN/.venv/bin/python -u build_low.py --triples-npz $P/triples_full_r4.npz \
    --out $P/pieces/low_r4.npz

# 3. the HIGH windows, queued
python run_queue.py --jobs-json jobs_bands.json --max-jobs 5 --min-free-gb 30

# 4. assemble any (rows, cutoff, model) combination
$CAN/.venv/bin/python assemble.py --low $P/pieces/low_r4.npz \
    --bands "$P/pieces/band_tree_r4_"*.npz --cutoff 3840 \
    --subset-npz $P/triples_full_r1.npz --out $P/table_tree_cut3840_r1.npz

# 5. fold
sh sweep_cutoff.sh tree 960 1920 3840 7680 15360
```

## Files

| file | what it does |
|---|---|
| `_common.py` | shared loaders, octave windows, grid fingerprints |
| `build_low.py` | the LOW exact-Wigner-3j branch, once per grid |
| `build_band.py` | one HIGH multipole window, chunked over rows |
| `assemble.py` | LOW + selected windows + row subset -> deployed-format table |
| `run_queue.py` | concurrency and free-memory guarded job runner |
| `observables.py` | the four plotted observables from a folded sweep |
| `analyse_cutoff.py` | partial sums against cutoff, and the k-space mapping |
| `probe_weight_identity.py` | is the published curve the interpolation weight |
| `probe_permutation.py` | do the spin channels depend on the leg ordering |
| `probe_radial_interp.py` | cost of linear interpolation between the shells |
| `probe_eb_structure.py` | which vertex channel carries E minus B, and where it stops being small |
| `table_nonlinear.py` | tree against BiHalofit by cutoff; emits the note's LaTeX table |
| `fig_redshift.py` | FK against source redshift, and the cost of the lambda floor |
| `cl_ratio.py` | FK relative to Order-0 in harmonic space |
| `fig_*.py`, `plotstyle.py` | the note's figures |
| `sweep_cutoff.sh` | assemble and fold a series of cutoffs, one fold at a time |

## Traps

* `chi_shells` in a table is comoving distance in **Mpc**, not Mpc/h. The
  bispectrum is in h units, so the multipole-to-wavenumber map needs
  `chi_h = chi_Mpc * h`. Getting this wrong misplaces every mode by 1/h.
* The callable deduplicates sorted triples after rounding to 8 decimals, so
  collapsed configurations closer than `1 - cos g = 5e-9`, that is
  `g = 0.344 arcmin`, cannot be represented however finely the table samples.
* The production config's `output.path` and both `cache_path`s are absolute
  and point into the deployed run directory, and the runner unlinks the
  output before running. `run_fk_variant.py` redirects all three and refuses
  to run if the output resolves into production.
* P(k) stops at k = 206 h/Mpc and canoes refuses to extrapolate, so the
  innermost usable shell is `chi_h = 2 ell_max / k_max`.
