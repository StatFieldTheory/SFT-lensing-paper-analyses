# Prompt for a dedicated session: the normalisation of the driving-field
# three-point vertex

Hand this file to a fresh session. It is self-contained. Deliverable is a
short report, not a code project.

---

## The question

The SFT lensing pipeline builds a three-point cumulant of the driving
fields, `zeta_XYZ(gamma, chi)`, from the matter bispectrum, through a
chain of per-leg factors: the Poisson relation, a screen-space Hessian
projection, a Sachs photon weight, a Limber radial collapse, and a
measure conversion between comoving distance and affine parameter. That
table then feeds the FK contribution to the lensing two-point function.

Every ingredient of that chain has been checked for internal consistency.
What has NOT been checked is its ABSOLUTE NORMALISATION against anything
external. Specifically:

* the Order-0 comparison against PyCCL validates the TWO-point transfer
  chain, and never touches the three-point vertex
* the fastnc closure test validates the three-point machinery on OPEN
  triangles, and cannot reach the collapsed configurations the FK term
  uses
* the Monte-Carlo cross-check reads the same table as the production fold,
  so it validates the diagram assembly and not the input

Your task is to establish whether the normalisation is right, and if not,
by how much and where the error enters.

## The method, in two steps

**Step 1, analytic, and it is the decisive one.** The pipeline's own
factors should, when assembled into a convergence bispectrum, reproduce
the standard Limber expression

    B_kappa(l1, l2, l3) = INT dchi  W(chi)^3 / chi^4
                          * B_delta(l1/chi, l2/chi, l3/chi; z(chi))

with the usual lensing efficiency
W(chi) = (3/2) Omega_m H0^2 (1+z) chi (chi_s - chi) / chi_s.

Do NOT assume this and fit to it. Assemble it from the pipeline's pieces
independently: the per-leg Poisson amplitude, the per-leg screen-Hessian
multiplier, the per-leg Sachs weight, the 1/chi^4 Limber factor, the
affine lensing kernel, and the comoving-to-affine measure conversion. Then
compare, symbolically or numerically, against the standard form. Any
discrepancy in the prefactor, in a power of (1+z), or in a factor of h IS
the normalisation error, and this step localises it to a specific factor.

Two features of the pipeline are worth verifying rather than trusting:
under Limber the screen-Hessian multiplier is said to cancel the Poisson
1/k^2 per leg, leaving a scalar; and the collapse of two radial integrals
is said to carry a Jacobian (dlambda/dchi)^2 = (1+z)^-4. Check both.

**Step 2, empirical anchor.** Predict the convergence bispectrum with the
pipeline and compare against published measurements on ray-traced
lightcones, which depend on no fitting formula.

The comparison already exists in print and, usefully, covers squeezed
configurations, which is the shape family the FK vertex integrates over:

* Takahashi et al. 2020, arXiv:1911.07886 (the BiHalofit paper), section 6
  and Figure 13. Section 6.2 shows the convergence bispectrum at source
  redshift z_s = 1 in equilateral, flattened, and two squeezed
  configurations, B_kappa(l, l, 84) and B_kappa(l, l, 510), binned
  dlog10 l = 0.13, measured from 1000 ray-traced maps of 5 x 5 deg^2
  (Sato et al. 2009; Kayo, Takada, Jain 2013, arXiv:1207.6322). Section
  6.1 does the same for CMB lensing against Namikawa et al. 2019,
  arXiv:1812.10635, in the Takahashi et al. 2017 full-sky simulations,
  arXiv:1706.01472, with a squeezed panel B_kappa(l, l, 50).

PDFs of all of these are in `papers/`. The reference note is
`notes/route1_convergence_bispectrum_references.md`.

Reproduce the z_s = 1 cosmic-shear squeezed panels with the pipeline's own
projection, using the same matter bispectrum the published curves use
(BiHalofit), so that only the projection differs. Agreement validates the
transfer chain; disagreement localises the error, and the analytic step
should already tell you which factor.

**Step 3, if the analytic factorisation applies.** The squeezed-limit
factorisation in `notes/theory_shape_amplitude_factorization.md` predicts
the vertex in closed form inside its validity window. Comparing the
pipeline's numerical vertex against that prediction tests the same
normalisation without going through a projection to the convergence
bispectrum at all, so a disagreement here versus in Step 2 would separate
an error in the vertex itself from an error in the projection to
B_kappa. Report both if you run both.

## What counts as an answer

Report the ratio of the pipeline's B_kappa to the published measurement,
as a function of multipole, for at least the two squeezed configurations
and one equilateral one. State whether that ratio is consistent with unity
within the quoted simulation accuracy (5% up to l ~ 4000 for the Kayo
maps), and if it is not, identify the factor responsible.

Do not tune anything to make it agree. If the pipeline is off by a clean
factor, say so and name the factor. If it is off by something
scale-dependent, that is more informative still.

## Scope limit worth stating in the report

A measured convergence bispectrum lives on resolvable triangles with a
finite soft leg. It cannot reach the coincident limit the FK vertex
actually samples, so this study can validate the transfer chain and the
projection at finite soft leg, and cannot by itself certify the collapsed
corner. Say how close to the corner the published measurements let you
get, and what residual gap remains.

## Where things are

* the formalism: `sections/appendix.tex`, appendix
  `append: driving-field spectra`, and `sections/cosmology.tex`,
  subsection `subsec: driving fields`. Read-only; do not edit `sections/`,
  `figures/` or `main.tex`.
* the implementation: `canoes` at
  `/Users/zzhang/projects/angular_statistics/canoes`, importable from its
  own `.venv` with `PYTHONPATH=$CAN/src`. The Limber branch is
  `src/canoes/sachs/kappa3.py`, function `compute_kappa3_sigma3_high`.
  Its per-leg response product is `_kappa3_limber_response_product_h`.
* a matter bispectrum in the right units is already wrapped:
  `driver_field_emulators/code/b_model` serves
  `b_delta_fn(k1, k2, k3, z)` in h-units, models `tree` and `bihalofit`.
  Run its tests first to confirm the environment:
  `PYTHONPATH=$CAN/src:. $CAN/.venv/bin/python -m pytest b_model/tests -q`
* the extracted chain, factor by factor with file:line, is in
  `driver_field_emulators/notes/input_ready_b_to_zeta_spec.md`, sections 2
  and 5. Treat it as a map, not as evidence: it was written by reading the
  code, not by testing it.
* `driver_field_emulators/notes/theory_shape_amplitude_factorization.md`
  derives an analytic factorisation of the collapsed vertex in the
  squeezed limit, into a projected linear correlation shape times a
  response-weighted hard variance, with coefficients machine-verified and
  the result confirmed by direct quadrature. This gives you a THIRD handle
  on the normalisation, and possibly the sharpest one: the factorised form
  has an explicit closed-form amplitude, so the pipeline's numerical
  vertex can be checked against it directly, factor by factor, without
  any simulation. Use it as a cross-check on Step 1, and note where its
  validity window ends.

## Traps

* `canoes`'s `pip` script has a broken shebang; use
  `.venv/bin/python -m pip`.
* The P(k) table stops at k = 206 h/Mpc and canoes refuses to extrapolate.
* If you run anything through sft-wick, use
  `driver_field_emulators/code/run_fk_variant.py`. The production config
  names an absolute output path and the runner unlinks it before running,
  which destroys the deployed result.
* Units: canoes works internally in Mpc/h and converts with h^4 at the
  end. The published measurements are in the usual convergence-bispectrum
  conventions. Getting this wrong will look like a normalisation error, so
  pin the convention explicitly before concluding anything.

## Deliverable

A short report, two to four pages, in
`driver_field_emulators/notes/`. Plain language. Structure:

1. what was compared and how
2. the analytic assembly and whether it reduces to the standard Limber
   form, with any discrepant factor named
3. the numerical comparison against the published measurements, with the
   ratio plotted or tabulated
4. the verdict on the normalisation, and the residual gap at the collapsed
   corner
5. commands to reproduce

If the answer turns out to be "the normalisation is correct to within the
simulation accuracy", that is a complete and valuable result. Say it
plainly and stop.
