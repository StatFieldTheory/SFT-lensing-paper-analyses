# driver_field_emulators

Can an existing emulator supply the three-point input statistics of the SFT
lensing driving fields (Ricci focusing Phi_00, Weyl shear Psi_0), and how does
one convert between an emulator's output and what the pipeline consumes?

Survey date: 2026-08-25 (multi-agent literature sweep, per-candidate source
verification, completeness check; coverage 2018-2026).

> **Next session: read `HANDOVER.md`.** It carries the open question
> (does the emulator give a more realistic prediction), the commands, and
> the traps. Its status update of 2026-08-26 marks B4 (the redshift
> trend) measured at ell <= 1000 and adds two rebuild requirements: the
> cutoff fixed in k, and the lambda-shell grid or its interpolation
> fixed.
>
> **Read `SUMMARY.md` first.** It is the conclusion of the work: what was
> found, what it means for the draft, and what should happen next.

## Bottom line

1. **No emulator of the absolute gravity-only matter bispectrum
   B(k1,k2,k3,z) exists** (as of 2026-08). The community backbone remains the
   BiHalofit fitting formula (Takahashi et al. 2020, arXiv:1911.07886,
   ~10-15% at k < 3-10 h/Mpc, calibrated z = 0-3).
2. **The only genuine 3-point emulator in the Euclid ecosystem** is the
   baryonic matter-bispectrum BOOST emulator (Burger et al.,
   arXiv:2506.18974, A&A 705 A170): R = B_hydro/B_gravity-only, shipped as
   `baccoemu.Matter_bispectrum()` (`pip install baccoemu`); k in
   [0.017, 17.7] h/Mpc, a in [0.24, 1.01] (z <= ~3.2). The gravity-only part
   is not emulated; the paper multiplies the boost onto BiHalofit.
3. **No Weyl/potential/Ricci-field 3-point emulator exists**, and none is
   needed: the driving fields are per-multipole LINEAR transfers of the one
   scalar potential, so any matter-bispectrum source converts exactly
   (see notes/conversion_interface.md).
4. Lensing-level 3-point emulators exist but are observable-level surrogates
   (aperture-mass skewness <Map^3> for KiDS/DES/HSC, integrated shear 3PCF),
   all trained on BiHalofit-based theory codes; useful as methodology
   references (CosmoPower/cosmoemuJAX training patterns), not as inputs for
   this pipeline. The one simulation-trained lensing-bispectrum surrogate
   (MassiveNuS GP, Coulton et al. 2019) was never publicly released.
5. **Recommended input model for SFT lensing**:
   B = BiHalofit(k, z) x baccoemu_boost(k, z) inside the training box, with
   tree-level SPT at low k (exact there) and tree/BiHalofit fallback outside
   the box. The equal-time restriction of emulators is NOT a structural
   obstacle: the kappa3 vertex is already equal-shell collapsed, and the
   surviving cross-shell structure is restored by per-leg growth transfers
   (see notes/cross_shell_vs_equal_time.md).

## Notes index

- `notes/input_ready_b_to_zeta_spec.md` - THE implementation spec: the
  complete source-verified map B_model(k1,k2,k3;z) -> configuration-space
  six-channel zeta table -> sft-wick coupling_fn, with quadrature, kernels,
  units, insertion diff, query set, grid requirements, and contract
  checklist. Start here for coding.
- `CHANGES_OUTSIDE_THIS_FOLDER.md` - the complete audit of what this
  session touched beyond this folder, and how to revert each item. The
  draft (sections, figures, main.tex) was never modified.
- `notes/theory_redshift_dependence.md` - the redshift dependence made
  explicit: the exact per-shell growth-times-geometry factorization
  certified on the real table to 0.5%, the measured z_s trends (tree
  FK/O0 flat to -16%/+5% in z_s; nonlinear FK/O0 FALLING with z_s, 2.7%
  at z_s = 0.5 vs 0.8% at z_s = 5), the four measured distortions of the
  deployed multi-z trend, and a new fold systematic: the production
  lambda-shell interpolation inflates xi_FK at z_s = 5 by +11% at
  arcminutes, decaying through zero near 2 degrees (grid-converged tree
  baseline +3.84e-6, not +4.28e-6; all nonlinear/cutoff percentages are
  at ell <= 1000).
- `notes/theory_shape_amplitude_factorization.md` - the theory: the
  squeezed-limit factorization of the collapsed vertex into (projected
  linear correlation shape) x (response-weighted hard variance), with
  every coefficient machine-verified (wolframscript + independent sympy)
  and the theorem itself confirmed by direct quadrature to ~1%, including
  the sign-change and 17/7 - n/2 coefficient predictions.
- `notes/finding_three_effect_decomposition.md` - measured at the
  observable level: the ell_high_max cutoff changes xi_FK by x1.97 per
  doubling, MORE than the tree-to-BiHalofit upgrade (x1.66), while the
  cosine grid governs the wide-angle shape. Read this before quoting any
  FK amplitude.
- `notes/finding_grid_artifact_measured.md` - the deployed cosine grid
  freezes the FK vertex at the (1,1,1) corner from 0.5' to about 600',
  which is what produces the reported wide-angle FK enhancement and the
  crossover with Order-0 near 2 degrees. The 0.5 arcmin amplitude
  survives the GRID fix at the 0.15% level (the separate lambda-shell
  interpolation systematic of `theory_redshift_dependence.md` still
  applies to it).
- `notes/finding_ell_max_resolved.md` - the cutoff question settled: the
  vertex integral converges, and the deployed ell_high_max = 1000 sits
  below the integrand's peak, undercounting by factors of 1.5 to 57.
- `notes/finding_cl_and_credibility.md` - the same story in harmonic
  space, where the ell dependence REVERSES (the deployed FK is concentrated
  at ell < 10 and exceeds Order-0 there by 18x; the corrected FK rises with
  ell and stays at 1 to 2%), plus the assessment of how much the
  underconverged loop undermines, which is the amplitude and not the shape.
- `notes/finding_lambda_min_cut.md` - the LOS floor: 0.1% at z_s = 5,
  but 8 to 13% of the kernel weight at z_s = 0.5, so it can fake the
  redshift trend of the multi-redshift figure.
- `notes/finding_mc_crosscheck.md` - the independent diagram assembly
  agrees to 0.1% on BOTH the flat deployed curve and the decaying dense
  one, which is what proves it validates the fold and not the input.
- `notes/finding_ell_max_uv_sensitivity.md` (SUPERSEDED) - measured 2026-08-25: the
  HIGH-branch cutoff sets the vertex value rather than converging, because
  the collapsed family carries a UV-divergent coincident-pair moment.
  Affects the deployed table, so read before quoting any FK amplitude.
- `notes/finding_measurement_addenda.md` - three findings that came out of
  checks run for other purposes: the deployed build's own quadrature is off
  by ~10% at wide angles (banding is more accurate), what the nonlinear
  bispectrum does at the vertex before folding (20 to 46x at the near shell,
  with sign changes), and an upstream constructor bug in fastnc 1.2.3.
- `notes/cosmological_relevance.md` - is a one to two percent correction
  significant? It sits at the Stage-IV threshold, but being flat in angle
  it is largely degenerate with sigma_8, so what decides its importance
  is the redshift dependence, not the size; that dependence is now
  measured (`theory_redshift_dependence.md`): flat at tree level,
  falling with z_s for nonlinear input.
- `notes/assumptions_and_limitations.md` - the four-layer register of every
  assumption and model limitation (physics / B_model ingredients /
  numerics / software-provenance), each with its guarding test, plus the
  known-unknowns list. Independently verified 2026-08-25.
- `notes/emulator_candidates.md` - the verified candidate landscape: every
  genuine 3-point emulator found, with domains, output forms, code status,
  caveats; closest relatives and confirmed negative results.
- `notes/conversion_interface.md` - the exact conversion formalism and the
  two concrete integration points in the pipeline (NPZ table-level and
  canoes build-level), with unit/measure bookkeeping.
- `notes/cross_shell_vs_equal_time.md` - why single-redshift emulators match
  the pipeline's equal-shell collapse; power-conservation analysis; lifts to
  cross-shell if ever needed.
- `notes/coverage_ell_z_mapping.md` - (ell, z) -> k mapping tables; why the
  soft leg at the (z = 5, ell = 10) corner is tree-exact while the coincident
  hard pair is not; hybrid evaluation scheme.
- `notes/implementation_roadmap.md` - phased plan for the nonlinear-B FK
  rebuild (the "redo the 3pt leakage into the 2pt function" project).

## Where each kind of record lives

* `notes/` - physics and measurement findings, one file per finding
* `code/README.md` - implementation-level records: the sigma8 pins
  (0.808988 for the 3-point table, 0.810000 for the 2-point side), the
  measured BiHalofit limits (its 2 to 4% offset from tree below
  k = 1e-3 h/Mpc, the squeezed boost table), the environment facts and the
  canoes patch
* `CHANGES_OUTSIDE_THIS_FOLDER.md` - everything touched beyond this folder
* `HANDOVER.md` - the open questions and how to run them
* `PROMPT_normalisation_study.md` - a self-contained brief for a dedicated
  session on the one gap no existing check covers: the absolute
  normalisation of the three-point transfer chain

## Folder layout

- `papers/` - key papers (arXiv PDFs).
- `notes/` - analysis notes (above).
- `code/` - future integration code (baccoemu/BiHalofit wrapper for the
  equal_time_limber vertex rebuild).
