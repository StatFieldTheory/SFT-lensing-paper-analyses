# 2PCF / C_corr_op_K_limber_FK_cut15360_permfix

> **Result arrays archived on 2026-10-02.** The arrays of this run are superseded and were moved to
> `_archive/cleanup_2026-10-02/` (index: `reproduce/archived_inputs.json`). The configuration stays here
> because the referee-revision drivers in `sachs_sft/analyses/r1_sft061/` use it as their template. The
> manuscript's products are the ones selected by `reproduce/active_products.json`. The text below
> describes the run as it stood and is kept as a record.

**The FK two-point sweep of the manuscript** (figures 2, 3, 4, 6 and every FK number, ell_max = 15360).

* `xi_C_corr_op_K_limber_FK_cut15360_permfix.npz`: byte-identical to
  `driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz`, produced 2026-08-26
  09:48 by `run_fk_variant.py table_permclosed_cut15360.npz --callable perm_aware_kappa3_callable.py
  --tag permfix --n-jobs 6` (fold time about 2.5 min, sft-wick env). Order-2 kappa-kappa at 0.5' is
  `+1.9480e-5`, 2.31% of Order-0.
* `config_L2.yaml`: the same run expressed with paths relative to this folder (vertex =
  `../../../callables/kappa3_vertex/equal_time_limber_cut15360_permaware/`); `config_record_2026-08-26.yaml`
  is the verbatim config of the August run (absolute paths of that machine).
* `PRODUCTION`: marker read by `scripts/run_FF_single.py`, which refuses to overwrite the result here
  unless `SFT_WICK_FORCE_OVERWRITE=1`.

Re-run (writes into this folder; the runner clears the two `.cache_*` directories first):

    cd sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix
    SFT_WICK_FORCE_OVERWRITE=1 SFT_WICK_CONFIG_YAML=$PWD/config_L2.yaml SFT_WICK_EXPAND_ORDERS=0,2 \
      SFT_WICK_SWEEP_N_GAUSS=24 /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
      ../../../scripts/run_FF_single.py ../../../scripts/inputs/kappa2_grid_lcut0_low_empty.npz \
      ./xi_C_corr_op_K_limber_FK_cut15360_permfix.npz

The superseded cut1000 sweep of the June draft is `../C_corr_op_K_limber_FK/`.
