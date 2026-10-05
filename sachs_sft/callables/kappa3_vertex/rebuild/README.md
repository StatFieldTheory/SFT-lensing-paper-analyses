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
| `assemble_flip_odd_high.py` | the same sum with the spin-2 sign convention of every band checked; `--flip-odd-high` reverses the odd HIGH channels of old bands |
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

## Spin-2 sign convention (2026-10)

canoes before commit 810b133 built the HIGH branch with `response_psi =
-response_phi`. In the `exp(2 i phi_l)` basis that branch uses, `Psi0(l) =
+exp(2 i phi_l) Phi00(l)`, so every HIGH band built before that commit has
`zeta_TTP`, `zeta_PPP` and `dmod` (one or three spin-2 legs) with the wrong sign;
`zeta_TTT`, `zeta_TPP` and `bmod` are right, and the LOW (exact) branch is not
affected. A band's convention is in `cosmo_meta["potential_convention"]`:
`Psi0=-A(a)` is the old one, `Psi0=+A(a)` the corrected one. Every table in
`products/` and in `r1_sft061/` assembled before 2026-10-03 carries the old
sign in its odd channels, through its HIGH part.

`assemble.py` cannot tell the two kinds apart; `assemble_flip_odd_high.py`
reads the label of every band and stops on a mixture or an unknown label.
`pieces_nphi512/table_permclosed_np512_spin2fix.npz` is the stored production
table with the odd HIGH channels reversed (`--flip-odd-high`, no new
quadrature): the control that isolates the sign fix (paper task T-001). Bands
rebuilt with canoes 810b133 or later are assembled with the same script
without the flag.

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
* Two linear P(k) tables exist in canoes `examples/data`: `PCAMBz0.txt`
  (sigma8 0.808988, no header), which every piece built before 2026-10-03
  used, and `PCAMB_pyccl_stf_fid_z0.txt` (sigma8 0.810000, cosmology in its
  header), which the C table, FF, Order 0 and the PyCCL reference use. Both
  share one k grid. `build_band.py --pk-table` selects the table for the
  bispectrum model, `pk_matter_today` and the BiHalofit sigma8, and records
  its name and sha256 in the build metadata; a piece without that record was
  built on `PCAMBz0.txt`. `merge_bands.py` and `assemble.py` refuse bands on
  different tables. The LOW piece is reused across choices, so an assembled
  table records `pk_table_low` and `pk_table_high` separately.
