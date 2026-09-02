# Handover

Written 2026-08-25. Three tasks, in order. Task A must be done FIRST and
BLIND: do not read the rest of this folder's findings before finishing it,
and do not tell the reviewer what to expect.

Background material, once Task A is complete: `SUMMARY.md`, then `notes/`.

## Status update, 2026-08-26 (read after Task A, like the rest)

A redshift-dependence round ran after this handover was written
(`notes/theory_redshift_dependence.md`; scripts
`code/{verify_redshift_scalings.wl, redshift_pershell.py,
redshift_fold_trends.py, probe_lambda_interp.py}`; data
`products/redshift/`). It changes the task list as follows.

* **B4 is MEASURED at ell <= 1000** by a cheaper route than prescribed
  below: the existing dense tables were folded at lambda_f(z_s) directly
  (their shells extend past lambda_f(5), so no per-plane rebuilds were
  needed), against the exact 20-plane O0. Result: the tree FK/O0
  fraction is flat in z_s to -16%/+5% (0.40-0.50% at 0.5 arcmin at every
  z_s in 0.5-5), and with BiHalofit input the fraction FALLS with z_s
  (~2.7% at z_s = 0.5 against ~0.8% at z_s = 5). On the cosmology
  argument of `notes/cosmological_relevance.md` this is the favourable
  outcome: at tree level sigma_8 absorbs the term, and the nonlinear
  z_s dependence points DOWN with z_s, not up. What remains of B4 is
  only the repeat at a converged cutoff, which folds into B1/B2.
* **A NEW fold systematic was found and must shape B1**: the production
  fold interpolates the vertex LINEARLY in lambda across the mid-path
  shell gaps, where the z_s = 5 fold carries half its weight, inflating
  xi_FK by +11% at arcminutes (decaying through zero near 2 degrees,
  -10% at 145 arcmin). A truth probe at mid-gap shells certifies
  log-PCHIP instead (within ~1%). Grid-converged tree cut1000 baseline:
  +3.84e-6, not +4.28e-6. Every absolute FK number in the older notes
  carries this; all variant RATIOS are unaffected. So the B1 rebuild has
  a SECOND grid dimension to fix: the lambda shells (densify mid-path,
  or make the callable interpolate in log), not just the cosine rows.
* **The cutoff must be fixed in k, not in ell** (new B2 requirement):
  at fixed ell the k-cap 1000/chi varies with the shell, so the
  undercount is z_s-dependent (fold-weighted 1.39 at z_s = 0.5 against
  2.14 at z_s = 5), and any fixed-ell cutoff redshift-biases a multi-z
  figure by ~1.5x across the range.
* **B3's low-z_s empirical test is DONE**: folding the floor-397 against
  the floor-50 table gives -6.6% at z_s = 0.5 falling to -0.1% at
  z_s = 5 (0.5 arcmin), confirming the kernel-weight estimate.
* Two convention traps quantified: t_final = 2313.029 sits 0.2% from
  lambda(z_s = 5) = 2318.0 in lambda but the fold is sensitive at
  ~1.1%/Mpc, so that 0.2% is a 5.5% amplitude difference (and 2313.029
  corresponds to z_s = 4.70 on the closed-form background). And the
  dense table's cosine rows interleave HALF-STEPS between the production
  40-gamma grid: nearest-match pairing across the two grids silently
  compares gammas 11% apart; always compare at the coinciding gammas.

---

# Task A: independent verification, blind

**Purpose.** A previous session reached a conclusion about the FK figures
in the draft. That conclusion needs confirming by someone who has not been
told what it is, working forward from the paper's own figure-production
path rather than backward from an answer.

**How to run it.** Give an independent reviewer (a fresh agent or a
colleague) the brief below, verbatim, and nothing else from this folder.
Do not paste any note, table or number from `notes/`. Do not mention the
words grid, interpolation, corner, or artifact.

## Reviewer brief (hand over exactly this)

> The paper `STF_lensing` reports an Order-2 "FK" contribution to the
> lensing two-point function, shown in `figures/analysis3_NLO_FFFK.pdf`
> and `figures/analysis3_cl_full_vs_O0.pdf`, and quoted in
> `sections/insights.tex` and `sections/conclusion.tex`.
>
> Your task is to audit the production path of those figures, end to end,
> and report whether the plotted FK curve is a faithful evaluation of the
> quantity the paper defines.
>
> Start from `SFT-lensing-paper-analyses/figure_manifest.md`, which lists
> each figure's generator and inputs. Work forward: generator script, the
> `xi_*.npz` it reads, the sft-wick run that produced it, the config, the
> callables that config names, and the tables those callables load. For
> the three-point vertex the relevant definition is in
> `sections/appendix.tex`, appendix `append: driving-field spectra`.
>
> Questions to answer with evidence:
>
> 1. What quantity does the FK fold actually request from the kappa3
>    vertex, as a function of the plotted separation? Determine this from
>    the code, not from the documentation.
> 2. How does the fold obtain that quantity from the tabulated data?
>    Characterise the mechanism and establish its accuracy.
> 3. Over the angular range the figures plot, how well does the pipeline
>    resolve the requested quantity's dependence on separation? Support
>    the answer with a measurement.
> 4. Are the numerical parameters of the vertex build (any multipole
>    cutoffs, quadrature rules, radial grids) adequate for the plotted
>    range? Test rather than assume.
> 5. Anything else in that path that would change the plotted curve.
>
> Report what you measure, with the commands used. If you find nothing
> wrong, say so; a clean audit is a useful result. Do not modify
> `sections/`, `figures/` or `main.tex`.
>
> Environment: canoes is at
> `/Users/zzhang/projects/angular_statistics/canoes`, importable from its
> own `.venv` with `PYTHONPATH=$CAN/src`. sft-wick runs in the `sft-wick`
> conda env. Warning: the production run config names an ABSOLUTE output
> path and `run_FF_single.py` unlinks it before running, so a naive re-run
> destroys the deployed result; copy the config and redirect its `output`
> and `cache_path` entries before running anything.

## After the reviewer reports

Compare their findings with `notes/finding_grid_artifact_measured.md`,
`notes/finding_ell_max_resolved.md` and `notes/finding_lambda_min_cut.md`.
Record agreements and disagreements. Only proceed to Task B if the
independent audit confirms the substance; if it does not, resolve the
discrepancy first.

---

# Task B: redo the three analyses to revision standard

Assume Task A confirmed the problem. The existing measurements were
exploratory: 80 cosine rows on the collapsed family only, one component
pair, one source redshift, cutoffs chosen for speed. A referee would ask
for more. Redesign each, then produce the note in Task C.

Four sub-tasks. B4 is listed last but decides whether the term matters
at all; it is MEASURED at ell <= 1000 (status update at the top), so
what remains of it is the converged-cutoff repeat on the B1/B2 grid.

## B1. Dense-grid rebuild, production scale

What exists is a collapsed-family-only table with 80 rows, valid for the
FK two-point fold and nothing else. What a revision needs:

* the full production triple set (1671 rows) with the collapsed family
  densified, so the same table serves the order-1 three-point vertex and
  the multi-redshift figure. Densify in log(1 - cos gamma), NOT uniformly
  in cos: the diagnosis is not that the grid is coarse in angle but that
  the queries crowd against the c = 1 endpoint, twenty-six of the forty
  sitting within 1% of it, so a uniformly finer cosine grid would not fix
  the problem (`notes/finding_grid_artifact_measured.md`)
* a resolution study rather than one grid: build at two or three
  densities and show the answer has stopped moving, since "we refined the
  grid and the answer changed" invites "how do you know it has converged
  now"
* the SECOND grid dimension (status update above): densify the lambda
  shells across the mid-path gaps (130-450 Mpc between 397 and 1823 Mpc,
  where the fold carries half its weight), or switch the callable's
  lambda interpolation to log-PCHIP, and show the fold is stable under
  both. `code/probe_lambda_interp.py` is the acceptance test: linear
  interpolation across the current gaps inflates the z_s = 5 fold by
  +11% at arcminutes
* all four observables, not just kk. The draft's figures show
  xi_kappa, xi_+, xi_- and xi_kappa_gamma, and the selection rule that
  sends FK into the first two is a claim the corrected table should be
  shown to preserve
* the old and new tables compared at fixed everything else, so the figure
  isolates the grid

Cost note: the LOW branch is independent of the cutoff and always
tree-level, so `--cache-low` lets one LOW build serve every variant. That
was the difference between 235 s and 103 s per variant at 80 rows; at 1671
rows the LOW build is the expensive part and the cache matters more.

## B2. ell_max convergence

What exists: octave-band slopes, and folded results at 1000, 2000 and
8192 for tree, 1000 for BiHalofit. What a revision needs:

* the nonlinear model carried to the same cutoffs as tree, which has not
  been done and is the main gap
* the amplitude as a FUNCTION of the cutoff, for both models, rather than
  a claimed converged number. The bands are additive, so one sweep gives
  every cutoff
* the k-space statement made explicit: which physical wavenumbers
  dominate at each cutoff, plotted against BiHalofit's calibration range
  and the scale where baryons matter. This is what turns "underconverged"
  into a quantified model-error statement
* the practical ceiling documented: the P(k) table stops at k = 206 h/Mpc
  and the Limber branch needs k = 2 ell_max / chi_h, so the innermost
  shell sets the maximum reachable cutoff. Either extend the table or
  state the limit
* a quadrature control, since banding and a single run at fixed n_ell do
  not agree at wide angles (see `notes/finding_measurement_addenda.md`)
* the cutoff parameterised in k, not in ell (status update above): a
  fixed-ell cutoff undercounts the far shells 3.4x more than the near
  ones (per-shell U runs 1.19 to 4.06 across z = 0.1-5.7), which makes
  the undercount z_s-dependent and redshift-biases any multi-z figure
  built at fixed ell

## B3. lambda_min cut

What exists: a kernel-weight estimate across source redshifts, and one
empirical test at z_s = 5 with shells extended to 50 Mpc. What a revision
needs:

* the empirical test repeated at low source redshift: DONE 2026-08-26,
  -6.6% at z_s = 0.5 falling to -0.1% at z_s = 5 at 0.5 arcmin (up to
  -15% at 145 arcmin, z_s = 0.5), by folding the floor-397 against the
  floor-50 table at each lambda_f
  (`notes/theory_redshift_dependence.md` section 4 item 3)
* the two effects that both grow toward low z_s separated: this cut, and
  whatever the grid does to slices covering different angular ranges
* the structural floor stated: chi_min = 2 ell_max / (k_max h), which
  couples the cutoff choice to the smallest reachable shell

## B4. The redshift trend, and why it comes FIRST

STATUS 2026-08-26: MEASURED at ell <= 1000 (see the status update at the
top and `notes/theory_redshift_dependence.md` section 3; data in
`products/redshift/fold_trends.npz`). The route taken was cheaper than
the one prescribed below: the dense tables' shells extend past
lambda_f(z_s = 5), so ONE table folds at every lambda_f(z_s), no
per-plane rebuilds. Headline: tree FK/O0 flat in z_s to -16%/+5%;
BiHalofit FK/O0 falls from ~2.7% (z_s = 0.5) to ~0.8% (z_s = 5). What
remains is the repeat at a converged, k-fixed cutoff on the B1 grid,
plus the open question of whether the deployed multi-z figure should be
replaced (a draft decision, out of scope here). The per-plane t_final
table below remains correct and useful for that repeat.

Priority note: on the cosmology argument in
`notes/cosmological_relevance.md` this is the most important of the four,
and it should be run before the amplitude work. The corrected FK is flat
at one to two percent of Order-0 across angle, and a flat multiplicative
offset is largely absorbed by sigma_8. What CANNOT be absorbed in a
tomographic analysis is a term that grows with source redshift. So the
draft's claim, that the leakage "strengthens steadily with source
redshift" (conclusion.tex:37, and the figure `multiz_kappa_xi_cl.pdf`), is
what decides whether this term belongs in a cosmological analysis at all.

It is also the one claim this session did not test, and two measured
effects both bias it in the same direction, making the trend look stronger
than it is:

* the cosine-grid freezing, whose severity depends on the angular range
  each source plane covers
* the lambda_min floor, which removes 8 to 13% of the kernel weight at
  z_s = 0.5 against under 1.5% at z_s = 5
  (`notes/finding_lambda_min_cut.md`)

**Do not reuse the existing multi-redshift machinery for FK.** The
figure's FK slice is PCHIP-interpolated in z from a cached grid under
`talk/`, built on the same coarse cosine table, so it inherits both
problems. Compute each source plane directly. (The 2026-08-26
measurement did exactly this split: it reused the cache's EXACT O0
planes and lambda_f list, which are clean, and folded the dense tables
itself for FK; the cache's FK slice was used only as the
"deployed trend" comparison curve.)

The table to use already exists:
`products/collapsed_dense_tree_cut1000_lowshells.npz` carries both fixes,
the densified collapsed family and shells extended down to lambda = 50 Mpc.
Rebuild it at the converged cutoff before the production run.

Each source plane needs its own `t_final`, the affine parameter at the
source, with dlambda/dchi = a^2 in physical Mpc:

| z_s | chi [Mpc] | t_final |
|---|---|---|
| 0.5 | 1959.1 | 1330.73 |
| 1.0 | 3413.8 | 1822.72 |
| 1.7 | 4855.8 | 2094.89 |
| 2.5 | 5989.9 | 2216.79 |
| 3.2 | 6716.0 | 2266.58 |
| 4.0 | 7356.1 | 2297.29 |
| 5.0 | 7971.2 | 2317.96 |

Cross-check: the deployed config uses t_final = 2313.0289 for z_s = 5,
which matches the table to 0.2% in lambda; but the fold is sensitive to
lambda_f at ~1.1% per Mpc, so that 0.2% is a 5.5% amplitude difference
(2313.029 corresponds to z_s = 4.70 on the closed-form background).
Whichever convention you pick, use it on BOTH sides of every ratio.
Derive the values you use from the pipeline rather than from this table,
and beware the unit trap: pass chi in physical Mpc, since
`lambda_of_chi` defaults to `distance_unit="Mpc"`. Passing Mpc/h
silently gives a number about 39% too large.

What to report: FK / Order-0 at each source redshift, on the corrected
table, with the trend shown before and after each of the two fixes so the
reader can see how much of the reported trend was real. At ell <= 1000
this is now answered: the tree trend flattens (sigma_8 absorbs it), and
the nonlinear trend runs OPPOSITE to the draft's claim, falling with
z_s. The converged-cutoff repeat decides whether those statements
survive quantitatively; the direction of the nonlinear one is expected
to (deepening the cutoff adds low-redshift small-scale power, which
strengthens exactly the low-z_s end).

## Cross-cutting

Every rebuilt table should be folded with `code/run_fk_variant.py`, which
isolates the run from the deployed artifacts, and cross-checked with
`code/mc_crosscheck_dense.py`. Note what that cross-check does and does
not establish: it agrees to 0.1% with the production fold on any table,
so it validates the assembly, not the input.

---

# Task C: the note, LaTeX and PDF

Write a self-contained note in this folder, `note/` alongside `notes/`,
building to PDF with `latexmk -pdf`. It is the document a referee report
would be built from, so: plain language, short sentences, no rhetoric,
every number traceable to a script in `code/`.

Suggested structure, about six to ten pages:

1. what was checked and how
2. the vertex sampling problem, with the evidence internal to the deployed
   table so a reader need not trust a second table
3. the cutoff analysis, with the amplitude-versus-cutoff curve and the
   k-space model-error statement
4. the lambda_min analysis across source redshift
5. corrected results: the figures below, before and after
6. what changes in the draft's text, quoted line by line, separating
   confirmed statements from those that do not survive
7. what remains uncertain

## Figure inventory: everything that must be covered

From `SFT-lensing-paper-analyses/figure_manifest.md`. Six of the sixteen
paper figures depend on the kappa3 vertex and need regenerating; one more
is worth re-examining; the rest are unaffected.

| # | figure | why it depends on the vertex |
|---|---|---|
| 2 | `analysis3_NLO_FFFK.pdf` | reads the FK `xi_*.npz` directly |
| 3 | `analysis3_cl_full_vs_O0.pdf` | same three `xi_*.npz` |
| 4 | `cl_EB_polarization.pdf` | reads the FF and FK `xi_*.npz` |
| 5 | `multiz_kappa_xi_cl.pdf` | two-stage; the compute step reads the equal_time_limber callable, and its cached FK slice comes from the same table. This is the figure B4 addresses, and it carries BOTH biases: the grid and the lambda_min floor |
| 6 | `appendix_mc_workflow.pdf` | reads the FK `xi_*.npz` and, through `driver_stats.py`, the kappa3 callable |
| 7 | `zeta_driving_field_slices_draft.pdf` | plots the vertex itself, via `ell_band_decomp_results.npz` |

Unaffected, for completeness: #1 `analysis1_O0_vs_pyccl.pdf` (Order-0
only), #8 to #13 (schematics and two-point operator slices), #14
(coordinate maps), #15 and #16 (Feynman diagrams).

Two cautions from the manifest. Figure 7's generator writes DIRECTLY into
the repository `figures/` directory with a hardcoded path, so running it
overwrites the paper's file; copy it aside first. Figure 5 depends on an
external cache under `talk/`, so its FK slice must be regenerated rather
than reused.

The note should show, for each of the six, the current and corrected
version side by side, or state explicitly why a given one does not change.

---

# Task D: the original open question

Once the above is settled, the question this programme started from:

> Does the emulator, that is a nonlinear bispectrum, give a more
> realistic prediction?

The comparison to run is tree versus BiHalofit at a CONVERGED cutoff on
the corrected grid; only BiHalofit at ell_max = 1000 has been done, which
is far below convergence. Build it with

```bash
CAN=/Users/zzhang/projects/angular_statistics/canoes
P=.../driver_field_emulators/products
PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1 \
  $CAN/.venv/bin/python -u build_equal_time_limber_table.py \
    --b-model bihalofit --ell-bands --ell-high-max 8192 \
    --triples-npz $P/collapsed_triples_dense.npz \
    --cache-low $P/_low_cache_collapsed_dense.npz \
    --out $P/collapsed_dense_bihalofit_conv8192.npz
```

The tree equivalent took 617 s. `baccoemu` is installed nowhere yet; the
`b_model` dispatcher has the slot for the Euclid baryon boost.

Answer in four separate parts, because they do not move together:

* **shape**: settled and cutoff-independent, and the nonlinear model does
  not change it
* **amplitude**: rises, and rises more at a higher cutoff, since halofit's
  extra power sits in the modes the cutoff was excluding
* **redshift dependence** (measured 2026-08-26 at ell <= 1000, the piece
  the emulator genuinely changes): the nonlinear enhancement of the fold
  falls from 7.2x at z_s = 0.5 to 1.7x at z_s = 5, because the
  nonlinearity lives in the z <~ 1 shells and even the z_s = 5 fold
  takes half its weight below shell z = 0.88. So the emulator's main
  physical effect is to make FK a LOW-source-redshift phenomenon
  (`notes/theory_redshift_dependence.md` sections 1.4 and 3.3)
* **trustworthiness**: a converged nonlinear FK is dominated by k of
  several to tens of h/Mpc, outside BiHalofit's calibration and where
  baryons are not a small correction. "More realistic input, less
  controlled dominant modes" is the honest phrasing

One limitation to state explicitly when answering. BiHalofit is not
independent of perturbation theory: its three-halo term is built on the
tree-level F2 kernel with a fitted effective spectrum, and only its
one-halo term carries genuinely new, simulation-calibrated information
(see `notes/emulator_candidates.md`, section 1b). So swapping tree for
BiHalofit measures how much the FK term grows once small-scale nonlinear
power is included. It does NOT provide an independent check that the
tree-level calculation was right. The only independent checks available
are a direct measurement of zeta on N-body lightcones, which nobody has
done, and the Monte-Carlo cross-check, which validates the fold rather
than the input.

What the comparison DOES support is a bounded statement, and it is worth
making carefully because it holds at the deployed cutoff and weakens at a
converged one. At ell_max = 1000 the hard legs reach at most k = 3.4 h/Mpc
(nearest shell), inside BiHalofit calibration, so the tree-to-BiHalofit
factor of 1.66 bounds the bispectrum-model contribution to the FK error
budget. At ell_max = 8192 the same shell reaches k = 28 h/Mpc, outside it,
so the bound stops applying exactly as the integral converges. Section 1c
of `notes/emulator_candidates.md` has the numbers and lists what does
support the model instead, together with the one gap: the collapsed slice
that FK uses has no independent normalisation check.

---

# Task E: the validation gap, if anyone wants to close it

Section 1d of `notes/emulator_candidates.md` argues that neither BiHalofit
nor baccoemu can validate the model, since both share the gravity backbone
and the whole projection chain, and lists what could:

1. **A convergence-bispectrum comparison.** Predict B_kappa(l1,l2,l3) with
   the same machinery and compare against a direct measurement on public
   ray-traced lightcones (MassiveNuS maps at Columbia Lensing, or the
   Takahashi et al. 2017 full-sky set). This independently tests the
   Poisson factor, the screen-Hessian normalisation, the Limber projection
   and the measure. It cannot reach the collapsed corner, so it does not
   close the gap, but nothing currently available gets closer.
2. **A direct measurement of zeta on N-body lightcones**, which does close
   it: build Phi_00 and Psi_0 from a simulated density field and measure
   the collapsed three-point cumulant. Nobody has done this.
   `/Users/zzhang/Workspace/stf-transfer` already computes driving-field
   two-point spectra from gevolution using the same transfer catalogue, so
   it is the natural starting point.

Neither is required for Tasks A to D. Both are what a referee would ask
for if they pressed on where the absolute normalisation of the collapsed
vertex comes from.

Route 1 is written up as a self-contained brief in
`PROMPT_normalisation_study.md`, ready to hand to a dedicated session. Its
key design point: the decisive test is analytic, not numerical. Assemble
the convergence bispectrum from the pipeline's OWN per-leg factors and
check whether it reduces to the standard Limber form W^3/chi^4; any
discrepant power of (1+z) or factor of h IS the normalisation error, and
that step localises it. The comparison against the published squeezed
measurements is the empirical anchor that follows. References and PDFs are
in `notes/route1_convergence_bispectrum_references.md` and `papers/`.

# Environment and traps

* canoes: `/Users/zzhang/projects/angular_statistics/canoes`, importable
  only from its own `.venv` with `PYTHONPATH=$CAN/src`. Its `pip` script
  has a broken shebang; use `.venv/bin/python -m pip`.
* The production config's `output.path` and both `cache_path`s are
  absolute and point into the deployed run directory, and
  `run_FF_single.py` unlinks the output before running. Use
  `code/run_fk_variant.py`, which redirects them and refuses to run if the
  output resolves into production. This destroyed the deployed FK result
  once; it was recovered from `_PROD_REF_2026-06-10`.
* P(k) stops at k = 206 h/Mpc and canoes refuses to extrapolate; fastnc
  does extrapolate, so the two branches differ above that.
* The draft was not modified in this session and should not be modified
  in the next one without an explicit instruction. Task C's note is the
  deliverable; edits to `sections/` are a separate decision.
