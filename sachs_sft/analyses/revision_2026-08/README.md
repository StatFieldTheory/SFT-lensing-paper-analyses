# revision_2026-08

The drivers that produced the 2026-08 revision figures by importing the paper's own generators and
rebinding their inputs (moved from `driver_field_emulators/code/figures_corrected/` on 2026-09-02):

| driver | what it made |
|---|---|
| `regenerate.py` | figures 2, 3, 4 with the FK input rebound to the cut15360 perm-aware sweep |
| `regenerate_multiz.py --reuse` | the FK planes of `analysis3/outputs/multiz_kappa_2pcf_5z.npz` (figure 5) |
| `run_ff_one_lambda.py`, `run_ff_multiz_serial.py` | the five real FF planes of figure 5 (`sftwick_outputs/2PCF/multiz/`) |
| `make_val_figure.py` | figure 6 with the corrected FK line and the `mc_fk_complete` markers |
| `regenerate_zeta_slices.py` | figure 7 at the converged cutoff (`equal_time_limber/outputs/zeta_bands_cut15360*.npz`) |
| `regenerate_mc.py`, `extend_ff_markers*.py` | superseded Monte-Carlo routes, kept for the record |

Since the generators now default to the corrected inputs (Phase 2 of the 2026-09 reorganisation),
`regenerate.py` is no longer needed to reproduce figures 2 to 4; it is kept as the record of how the
deployed files were made. Outputs go to `outputs/` here.
