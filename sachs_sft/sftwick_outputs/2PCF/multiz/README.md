# 2PCF / multiz

> **Result arrays archived on 2026-10-02.** The arrays of this run are superseded and were moved to
> `_archive/cleanup_2026-10-02/` (index: `reproduce/archived_inputs.json`). The configuration stays here
> because the referee-revision drivers in `sachs_sft/analyses/r1_sft061/` use it as their template. The
> manuscript's products are the ones selected by `reproduce/active_products.json`. The text below
> describes the run as it stood and is kept as a record.

Source-redshift sweeps behind figure 5 (`multiz_kappa_xi_cl.pdf`) and the redshift-trend notes.

* `order0_multiz_xi.npz` + `source_distances.npz`: exact Order-0 at seven source planes
  (z_s = 0.5, 1, 1.7, 2.5, 3.2, 4, 5; lambda = 1330.7 to 2318.0 Mpc), 2026-08-26, config
  `configs/config_order0_7planes_record_2026-08-26.yaml`. Note the figure's own planes (via the
  analysis3 `multiz_sweep.py` mapping) sit at lambda = 1822.721, 2094.894, 2216.793, 2266.582,
  2297.289 Mpc, so the z = 1.7 plane differs by 0.14 Mpc from this file's 2095.038.
* `multiz_ff_real_all5.npz` (and the per-plane files, plus the z = 5 control which reproduced the
  June production FF bit for bit): genuine FF sweeps at the five figure planes, one single-threaded
  process each, about 3.1 h per plane, 2026-08-27 (`analyses/revision_2026-08/run_ff_one_lambda.py`;
  logs `callables/kappa3_vertex/rebuild/logs/ff_z*.log`). `_ff_caches/` are their joblib caches (gitignored).
* The FK planes of the figure were folded with the perm-aware cut15360 vertex at the five figure planes
  (config `configs/config_fk_cut15360_permfix_5planes_record_2026-08-26.yaml`, 1034 s) directly into
  `analyses/analysis3/outputs/multiz_kappa_2pcf_5z.npz` by `analyses/revision_2026-08/regenerate_multiz.py --reuse`.
* The Order-0 planes of that npz are PCHIP-in-z interpolants of the 20-plane talk-era sweep vendored at
  `analyses/analysis3/inputs/multiz_components_talk_2026-06-10.npz`.
