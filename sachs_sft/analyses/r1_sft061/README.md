# Referee R1: sft-wick 0.6.1 numerical alignment

The accepted numerical products are now selected by `reproduce/active_products.json`.
See `closeout_validation.json` and the root `REPRODUCE.md`; the original
execution narrative below is historical. Agent run-state notes remain local.

The user requested alignment with the latest sft-wick after the Order-0
comparison was started. On 2026-09-28, a live `git ls-remote origin HEAD`
returned `6945815ad9e6efa2b9e768891d6bf92b1214670d`, the same commit as the
installed editable source at `/Users/zzhang/projects/SFT/sft-wick`.
The installed distribution reports version 0.6.1.

## Controlled first replay

`fk_same_inputs_gl24` keeps the manuscript's existing permutation-aware K
callable, deployed ell_max = 15360 input table, affine endpoint, 40-angle grid,
and 24-point Gauss rule. Only the framework version and fresh isolated caches
change relative to the historical product. It is driven by the existing
`run_fk_variant.py` with two workers. No production output is overwritten.

The release's internal covariance-table changes are bypassed by the production
configuration's `c_closed_form_only: true`. The release also fixes integration
domain splitting for two swept external times. Its effect on this particular
configuration must be checked, rather than inferred from its effect on package
demonstrations.

The historical reproduction pin is `53a3631` (`v0.3.0-8`), so the relevant
upgrade includes the intervening 0.6.0 changes. In particular, commit `9c000b4`
corrects callable nonlocal-vertex leg ordering in numerical FK evaluation.
`tests/test_nonlocal_leg_order.py` explicitly covers an equal-time FK example.
This is directly relevant to the component-dependent K callable used here.
The original August product does not record its executable commit, so the
replay comparison measures the total historical-to-current change, rather
than proving that one commit alone explains it.

## Separate input-table issue

`reproduce/provenance/vertex_rebuild_20260906.md` (relative to the repository
root) records a rebuilt K input with a
pair-phase correction and improved angular quadrature. That is a separate
change from the framework upgrade. The controlled replay above deliberately
keeps the old table so the two effects can be distinguished. The corrected
table and its acceptance evidence must be checked before adoption. Historical
figure products and data remain available for comparison.

## Original work plan (completed; retained for provenance)

- Compare the controlled replay with the historical FK arrays and record
  signed, component-level changes.
- Check the relevant integration-path behavior and numerical convergence.
- Establish the exact list of affected figures and quoted numbers.
- Recompute required products sequentially, review, and deploy through the
  reproduction tools after archiving previous artifacts.
- Update the manuscript and `referee-r1/response.md` after the numerical
  comparison and its limitations have been established.
