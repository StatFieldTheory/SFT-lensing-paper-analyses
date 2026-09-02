# Which FK kappa-kappa number is which

Written 2026-09-02 during the reproducibility reorganisation. The FK contribution to
`xi_kappa` at `gamma = 0.5'` appears in this repository as five different numbers,
because the pipeline was corrected four times. They are all "correct" for their own
configuration and none of them is a typo. This page is the key; every note that quotes
one of them links here.

| value at 0.5' | ell_max | vertex table | callable | lambda interpolation | where it is quoted |
|---|---|---|---|---|---|
| `+3.086e-5` | 1000 | pre-h^4-fix | sorting | linear | run READMEs, `equal_time_limber/README.md`, the production callable docstring, `gen_sftwick_configs.py`, `corr_op/README.md` |
| `+1.219e-4` | 1000 | post-h^4, pre-Jacobian | sorting | linear | `build_equal_time_limber_table.py` docstring |
| `+4.29e-6` | 1000 | deployed June table | sorting | linear | `notes/assumptions_and_limitations.md`, `notes/input_ready_b_to_zeta_spec.md`, `notes/implementation_roadmap.md`, `b_model/tree.py`, the August `code/README.md` |
| `+3.84e-6` | 1000 | deployed June table | sorting | log-PCHIP | `SUMMARY.md` section 4d, `HANDOVER.md`, `notes/theory_redshift_dependence.md`, `notes/finding_grid_artifact_measured.md` |
| **`+1.9480e-5`** | **15360** | **permutation-closed** | **permutation-aware** | **linear** | **the manuscript** (2.31% of Order-0), every figure since the 2026-08 revision |

Two corrections were applied to the table between the first two rows: the `h^6 -> h^4`
unit conversion and the `(1+z)^-4` collapse Jacobian, both landed 2026-06-09/10 on the
canoes side. The third row is the value the deployed June table actually produces once
both fixes are in. The fourth row is the same table folded with log-PCHIP interpolation
across the sparse mid-path lambda shells instead of the production linear interpolation,
which is a property of the fold, not of the table. The last row is the manuscript's: the
multipole cutoff raised from 1000 to the converged 15360 and the leg-permutation defect
of the lookup fixed.

**Rules that follow.**

* Never quote an FK amplitude without its `ell_max`. The cutoff changes it by a factor
  4.5 between 1000 and 15360, which is larger than any other correction in the table.
* The ratios between variants at a fixed cutoff are unaffected by the lambda
  interpolation, so the August notes' comparative statements stand as written even where
  their absolute numbers carry the `+4.29e-6` normalisation.
* The manuscript's value is reproduced by
  `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/`, whose README gives the
  command. Verified 2026-09-02: a fresh fold reproduces the deployed npz byte for byte.

Sources: `notes/finding_grid_artifact_measured.md` (the sampling defect),
`notes/finding_ell_max_resolved.md` (the cutoff), `notes/theory_redshift_dependence.md`
section 2.2 (the lambda interpolation), `note/main.pdf` (all four together).
