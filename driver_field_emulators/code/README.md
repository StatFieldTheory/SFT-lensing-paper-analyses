# code/ - the nonlinear-bispectrum FK rebuild

Status of the phased plan in `../notes/implementation_roadmap.md`.

## Environment

One environment runs everything: the canoes virtualenv, with canoes on
PYTHONPATH and the emulator code on the path as well.

```bash
CAN=/Users/zzhang/projects/angular_statistics/canoes
cd /Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code
PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q
```

Environment facts established 2026-08-25:

* canoes lives at `/Users/zzhang/projects/angular_statistics/canoes` (src
  layout). It is importable from its own `.venv` only with
  `PYTHONPATH=$CAN/src`; the editable installs in the PyCCL and sft-wick
  conda envs still point at the deleted `/Users/zzhang/projects/canoes`.
* That `.venv`'s `pip` script has a broken shebang from an earlier rename
  (`projects/AngStats`); use `.venv/bin/python -m pip` instead.
* BiHalofit comes from fastnc 1.2.3, installed with
  `python -m pip install --no-deps fastnc` (the full dependency set pulls in
  astropy, pandas and mpi4py, none of which the BiHalofit path needs). Its
  `halofit.py` is loaded directly as a module file, bypassing the package
  `__init__`, which does import the heavy modules.
* The three hardcoded `/Users/zzhang/projects/canoes` paths in the callable
  folder now resolve through `CANOES_ROOT` / the importable package / known
  locations instead.

## b_model package

`b_model` serves `b_delta_fn(k1, k2, k3, z)` in exactly the units the
canoes `equal_time_limber` HIGH branch consumes: k in h/Mpc, B in
(Mpc/h)^6, matter (no Poisson factor), equal time, one scalar z per call.

Models: `tree` (canoes-identical reference), `bihalofit` (gravity-only
nonlinear, blended to tree below the low-k seam), `bihalofit_baryon`
(times the Takahashi baryon ratio, only in fastnc builds that ship it).

Pinned inputs:

| quantity | value | provenance |
|---|---|---|
| sigma8 of `PCAMBz0.txt` (3-point table) | 0.808988 | direct top-hat integration; fastnc's independent integration agrees |
| sigma8 of `PCAMB_pyccl_stf_fid_z0.txt` (2-point side) | 0.810000 | same method |
| Omega_m, h, n_s | 0.3160919980475834, 0.6711, 0.97 | canoes defaults = STF fiducial |

Measured behaviour (see the test suite for the assertions):

* tree branch reproduces the canoes internal assembly to 6.6e-13 on its own
  quadrature geometry, and to 3.7e-6 through the full vertex, the residual
  being entirely the growth-factor interpolation route (canoes interpolates
  D on its lightcone z-grid, b_model on the growth ODE's a-grid). Pass
  `growth_fn` to remove even that.
* BiHalofit does NOT reduce to tree in the linear limit: it sits 2-4% high
  for all k <= 1e-3. The blend seam therefore sits at the crossover
  (5e-3 to 3e-2 h/Mpc), not deeper. In production the blend never
  activates, because the HIGH window forces k_max > 0.038 h/Mpc.
* squeezed boost B_nl/B_tree at (k_h, k_h, 1e-3): 1.4 / 3.2 / 7.8 / 16.4 at
  k_h = 0.2 / 0.5 / 1 / 2 h/Mpc for z = 0.1, softening to 1.08 at z = 5.7.
  These are larger than the P_NL/P_L proxy used in the earlier notes, since
  the 1-halo term contributes directly.
* cost is negligible: 7.9 s for the full production quadrature grid
  (96 x 96 x 64 nodes) over all 16 shells, 0.3 GB peak.

## canoes patch

Branch `feat/nonlinear-bdelta-vertex` in the canoes repository adds an
optional `b_delta_fn` to `compute_kappa3_sigma3_high` and
`compute_kappa3_mod_sigma3_high`. When supplied it replaces the internal
tree assembly AND bypasses the hardwired `growth**4`, because a nonlinear
bispectrum does not obey that separable scaling. `b_delta_fn=None`
reproduces the previous code path exactly. The helper `_eval_b_delta_fn`
fails loud on wrong shape or non-finite values, and the table metadata
records `b_delta_source`.

The LOW exact-Wigner-3j branch keeps tree level: its FFTlog machinery
factorises the bispectrum into per-leg one-dimensional integrals, which a
general nonlinear model does not admit.

## Build

```bash
CAN=/Users/zzhang/projects/angular_statistics/canoes
cd .../callables/kappa3_vertex/equal_time_limber
PYTHONPATH=$CAN/src $CAN/.venv/bin/python build_equal_time_limber_table.py \
    --smoke --b-model bihalofit
```

`--b-model` selects `internal_tree` (default, unchanged behaviour), `tree`,
`bihalofit` or `bihalofit_baryon`. The build records the choice in the
metadata and no longer carries the stale `+1.219e-4` baseline stamp: the
deployed tree baseline is +4.2877e-6, read from the production FK 2PCF.

## Diagnostics

`sweep_ell_high_max.py` separates quadrature resolution from genuine
multipole truncation in the HIGH branch, which a naive cutoff sweep at
fixed `n_ell` conflates.

## Redshift dependence (2026-08-25/26)

Four artifacts behind `notes/theory_redshift_dependence.md`, all run from
this directory with the canoes venv (`PYTHONPATH=$CAN/src:.`):

* `verify_redshift_scalings.wl` (wolframscript): 15 symbolic checks, the
  power-law soft/hard factors, the projection-slice identity, the
  per-shell exponent assembly, and the low-z fold closed forms with
  chi_s exponents verified on both sides of n = -1.
* `redshift_pershell.py`: the per-shell certification ([A2] brute force,
  the decisive bookkeeping test), the factorized-form accuracy vs shell
  redshift [A], the nonlinear excess [B] and the per-shell cutoff
  undercount [C]. Writes `products/redshift/pershell.npz`.
* `redshift_fold_trends.py`: the fold at lambda_f(z_s) on the dense
  tables (logPCHIP in lambda), controls T1 (production comparison at
  MATCHED gammas; the dense grid interleaves half-step rows, so nearest
  matching across grids pairs gammas 11% apart and fakes mid-angle
  inflation), trends T2/T3 and the fold-weighted diagnostics T4. Writes
  `products/redshift/fold_trends.npz`.
* `probe_lambda_interp.py`: the interpolation truth probe at mid-gap
  shells (linear/truth up to 1.17, logPCHIP/truth within 1.1%), the
  basis for the +11%-at-arcminutes production fold systematic and the
  +3.84e-6 grid-converged tree baseline.
