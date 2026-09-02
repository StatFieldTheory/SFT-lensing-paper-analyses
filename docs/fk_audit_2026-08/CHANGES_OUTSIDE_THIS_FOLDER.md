# Changes made outside driver_field_emulators/ during 2026-08-25

The session's deliverables live in this folder. The draft was never
touched: `git status` on `sections/`, `figures/` and `main.tex` is clean,
and no figure was regenerated. Three pipeline files and one separate
repository were modified, because building a table at all required it.
They are listed here so nothing is hidden and everything can be reverted.

## 1. sachs_sft callable folder (gitignored, so no VCS history to fall back on)

`SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/`

* `equal_time_limber_kappa3_callable.py`, `build_equal_time_limber_table.py`,
  `_local_cosmo_pk.py`: replaced the hardcoded `/Users/zzhang/projects/canoes`
  (a path that no longer exists, so these files could not import canoes at
  all) with a resolver that tries `CANOES_ROOT`, then the importable
  package, then known locations. This is a repair, not a feature.
* `build_equal_time_limber_table.py` additionally gained three optional
  flags, all defaulting to previous behaviour: `--b-model` (default
  `internal_tree`, the unchanged code path), `--triples-npz`, `--cache-low`.
  Its stale `fk_kk_baseline` docstring stamp (`+1.219e-4`, a pre-Jacobian
  value) was replaced with the measured `+4.2877e-6` and an explanation.

To revert: the flags and the resolver are additive and the default path is
byte-compatible, verified by the regression test in
`code/b_model/tests/test_vertex_integration.py`.

## 2. canoes repository (separate repo, isolated on a branch)

`/Users/zzhang/projects/angular_statistics/canoes`, branch
`feat/nonlinear-bdelta-vertex`, uncommitted. One file,
`src/canoes/sachs/kappa3.py`:

* optional `b_delta_fn` on `compute_kappa3_sigma3_high` and
  `compute_kappa3_mod_sigma3_high`; when absent, the previous code path
  runs unchanged
* `_eval_b_delta_fn` helper that fails loud on wrong shape or non-finite
  values
* `b_delta_source` recorded in the HIGH metadata, and
  `low_b_delta_source` / `high_b_delta_source` carried through
  `kappa3_combine_low_high`, which previously dropped the HIGH metadata

To revert: `git checkout master -- src/canoes/sachs/kappa3.py` or delete
the branch.

## 3. One production artifact was deleted and restored

`sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz`
was deleted by the first variant run, because the variant config inherited
the production config's absolute `output.path` and `run_FF_single.py`
unlinks that path before running. Restored from the repository's own
`xi_C_corr_op_K_limber_FK_PROD_REF_2026-06-10.npz`, and verified: the kk
order-2 value at 0.5 arcmin reads `+4.287723e-06`, the documented baseline,
with the original June timestamp preserved. The two `.cache_*` directories
in that folder were cleared and rebuilt by the same run; they are caches.

The cause is fixed in `code/run_fk_variant.py`, which now redirects the
output and both cache paths into the variant folder and refuses to run if
the output resolves inside the production directory.

## 4. Environment

* `fastnc` 1.2.3 installed into the canoes virtualenv with `--no-deps`
  (needed for BiHalofit; the full dependency set pulls in astropy, pandas
  and mpi4py). Removable with
  `.venv/bin/python -m pip uninstall fastnc`.
* Nothing was installed into the PyCCL or sft-wick conda environments.

## 5. One paper figure generator edited (2026-08-26)

`SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/plot_cl_EB_polarization.py`

The EB panel's annotation read `$\Delta C_\ell^{EB}=0$ (exact)`, which reads
as a measured null. It is a parity theorem that the vertex callable enforces
structurally: every coupling-tensor entry with an odd number of Im slots is
never assigned. Changed to `by parity`, with eight lines of comment recording
why. No numerical behaviour changes.

To revert: `git checkout -- <that file>`.

Every other figure regeneration this session left the paper's generators
untouched and only rebound their input paths, via
`driver_field_emulators/code/figures_corrected/regenerate.py`.

## 6. Paper figures replaced, old versions archived (2026-08-26)

`figures/analysis3_NLO_FFFK.pdf`, `figures/analysis3_cl_full_vs_O0.pdf`,
`figures/cl_EB_polarization.pdf` and `figures/multiz_kappa_xi_cl.pdf` were
replaced with versions built at `ell_max = 15360` on a permutation-closed
vertex table. The originals were moved with `git mv` to
`figures/archived_2026-08-26/` with a `_preRevision` suffix, so the history
is intact. For the multi-z figure only the FK slice was recomputed; the O0
and FF slices are bit-identical to the June result (neither reads the
three-point vertex), reused via `regenerate_multiz.py --reuse`. The refreshed
data npz is `sachs_sft/analyses/analysis3/outputs/multiz_kappa_2pcf_5z.npz`
(June original preserved at
`driver_field_emulators/products/multiz_kappa_2pcf_5z_preRevision.npz`).

`main.tex` gained an `\edited{...}` highlight macro and `xcolor`. Note that
REVTeX already defines `\revised`, so that name could not be used. Setting
`\revhlfalse` produces a clean copy.

`sections/insights.tex` gained one sentence in the C_ell figure caption,
inside `\edited{}`.

## 7. Second paper generator edited: FK Monte-Carlo markers removed (2026-08-26)

`SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/fig_xi_channels.py`

An exact second-order expectation analysis of `simulate_fk_vr` showed the FK
Monte-Carlo measures a placement share of the vertex (~0.24 at 0.5', crossing
zero near 10') plus a finite-sigma_lambda term sourced by the lambda-grid
graininess of the tabulated tensors; the apparent MC/analytic agreement at the
adopted sigma_lambda = 8 was those two terms summing through unity at the tuned
value. The FK markers were therefore removed from the validation figure
(_compute, _plot, argparse, and the superseded convergence docstring, which is
preserved in git history). The FK channel is validated instead by
`driver_field_emulators/code/rebuild/fk_kernel_crosscheck.py`, a direct
quadrature contraction of the tabulated cumulant with the FK kernel (agreement
~1% at 0.5'-5', ~7% at 17', vs the converged fold). FF markers keep the
original 15-point grid; the figure was rebuilt by
`code/figures_corrected/make_val_figure.py` and deployed to
`figures/appendix_mc_workflow.pdf` (June original archived in
`figures/archived_2026-08-26/`).

To revert: `git checkout -- <that file>`.

## 8. Third paper generator edited: multi-z FF slices are now real (2026-08-27)

`SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/plot_multiz_kappa_2x5.py`

The multi-z FF slices had never been computed. The talk-era cache filled them
with the exact z=5 FF scaled by the PER-GAMMA ratio fk(z)/fk(5)
(talk/scripts/build_multiz_from_cache.py step 5, docstring "ff(z):
NEGLIGIBLE"). At large separation fk sits on its +-1e-8 numerical floor, so
that ratio is noise over noise: the tails carried 1e-7-level dips, went
negative at z_s=1 (impossible, the FF moment is a connected part plus
<kappa>^2 > 0), and through the truncated-range transform produced a spurious
low-ell excess of order the Order-0 signal.

The generator now reads `driver_field_emulators/products/multiz_ff_real_all5.npz`,
five GENUINE per-redshift FF sweeps produced by
`code/figures_corrected/run_ff_one_lambda.py`: one single-threaded process per
source distance (sweep.n_jobs=1, so loky is out of the picture entirely),
about 3.2 h each, run concurrently with private expansion and propagator
caches. It fails loudly if that file is missing or its grids do not match.

Two independent checks on the result:

* each slice's large-separation plateau against the deterministic
  <kappa>^2(z_s) from `code/rebuild/mean_kappa_z.py` (the expectation of the
  simulate_ff_crn recursions, no Monte-Carlo): ratios 1.022, 1.005, 1.005,
  0.999, 1.008 at z_s = 1, 1.7, 2.5, 3.2, 4. Two unrelated computations
  agreeing to 0.1-2.2% on the statement that the FF two-point function tends
  to a constant at large separation.
* a z_s=5 CONTROL slice, run the same way with fresh private caches, against
  the June production FF: BIT-IDENTICAL at all 40 separations (max relative
  difference exactly 0.0, `np.array_equal` True). The production result was a
  real computation, it is deterministically reproducible three months on, this
  path is the same path, and no stale cache is involved.

An interim version of this figure used the z=5 FF scaled by the computed
plateau ratio <kappa>^2(z)/<kappa>^2(5); that approximation was accurate to
1.3% at large separation and 4.6% at 0.5' but up to 22% near 10', and has
been replaced by the real sweeps. The figure's O0, FF and FK are now all
genuine computations with no placeholder or scaled component.

To revert: `git checkout -- <that file>`.

## 9. C_ell generator: ratio strips converted to percent (2026-08-26)

`SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/plot_analysis3_cl_decomposition.py`

The (full-O0)/O0 ratio strips under each panel now plot percentages with
y-range 0.1-100% (previously fractional, 1e-3-1e2). Cosmetic units change
only; figure regenerated through the corrected-input driver and deployed.

## 10. Fig. 16 (zeta slices) rebuilt at the converged cutoff (2026-08-27)

`figures/zeta_driving_field_slices_draft.pdf` was the last product still drawn
at ell_max = 1000, while every other number and figure in the revision uses
15360. That is not just a normalisation: at shell 8 the normalised zeta_TTT at
gamma = 5' moves from -0.87 to -0.27 and the peak |zeta| grows by 13.7x, so
the printed panels showed a cumulant falling about three times more slowly in
angle than the one the paper's numbers come from.

Rebuilt by `driver_field_emulators/code/figures_corrected/regenerate_zeta_slices.py`,
which imports the production `ell_band_decomp.py` and `plot_zeta_figures.py`
and rebinds their constants rather than forking them: the exact Wigner-3j LOW
ladder stays at ell <= 60 and four flat-sky Born-Limber HIGH windows are
appended, reaching (8000, 15360]. Cost: ~2.5 min for LOW, ~20 s per HIGH
window.

The angular axis is clipped at 85'. Justification, recorded here rather than
in the caption: the single wide (60,15360] window and the sum of the eight
disjoint windows agree to 3-6% out to 85' but disagree by factors of several
to several tens beyond 117', where |zeta| is already below 2% of its peak --
the cumulant there is the residual of a near-cancellation between the low- and
high-multipole parts. Verified against the production cut15360 vertex table:
the plotted `full_*` agrees to 0.3% over gamma <= 21'.

The caption's blanket "change sign with angular scale" was also wrong at the
converged cutoff: within the plotted range zeta_TTP never changes sign,
zeta_TTT and zeta_Bmod change sign at 59.3' only for z >~ 4, and zeta_Dmod
from z ~ 1. It now says "change sign near 60' in the higher source shells".

Old version archived at `figures/archived_2026-08-26/zeta_driving_field_slices_draft_preRevision.pdf`.
No production file was edited; to revert, restore the archived PDF.

## 11. FK series de-emphasised outside its converged range (2026-08-27)

Five files: `analyses/analysis3/_plot_style.py`,
`analyses/mc_sachs_2pt/_plot_style.py` (a separate copy),
`analyses/analysis3/plot_analysis3_nlo_decomposition.py`,
`analyses/analysis3/plot_analysis3_cl_decomposition.py`,
`analyses/analysis3/plot_cl_EB_polarization.py`,
`analyses/mc_sachs_2pt/fig_xi_channels.py`.

The FK multipole sum converges only for gamma <~ 1 degree. The separation
enters zeta only through P_l3(cos gamma), which is unity for every term at
coincidence but oscillates with period 2*pi/gamma once l3*gamma >~ 1; beyond a
degree the high-l terms alternate in sign and largely cancel against the low-l
branch, and raising the cutoff shifts the residual instead of settling it
(measured: the ratio between the ell_max = 7680 and 15360 folds is 1.08-1.14
for gamma <= 56' but 0.69, 0.28, 2.59, 0.73 at 114', 294', 597', 1535').
Two independent evaluations of the same integral differ by a factor 5.2 at
114' while agreeing to 0.4-6% below 20'.

Rather than clip the axes -- Order-0 and FF are reliable over the full range
and carry a real result there (FF overtakes Order-0 near 190') -- the FK
series alone is drawn faint outside its converged range: gamma > 60' in
Figs. 11 and 17, ell < 50 in Figs. 12 and 14. Nothing is hidden and nothing is
over-claimed.

Implementation: both `_plot_style.py` copies gained an optional
`faint_outside=(lo, hi)` argument on `plot_signed_line` and (in the analysis3
copy) `plot_signed_markers`. It defaults to None, which reproduces the
previous single-style call exactly, so no other figure changes. The four
generators pass it for the FK series only.

To revert: `git checkout -- <those files>` and rebuild.

## 9. FK Monte-Carlo markers restored (2026-08-28)

Reverses part of item 7, which removed them. Item 7's reasoning was correct and
still stands: `sachs_mc_core.simulate_fk_vr` measures a placement share of the
vertex (0.211 at 1', crossing zero between 8' and 12'), and it must never be
used for FK markers again. The markers are restored from a DIFFERENT estimator,
`driver_field_emulators/code/mc_fk_complete/`, which is placement-complete (all
three legs of the deformation reach the vertex; proved as a symbolic identity in
generic symbols, `review/proof_placement_identity.wl`) and variance-controlled.

Two files edited outside this folder:

* `SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/fig_xi_channels.py`
  -- `_plot` now draws FK markers when the caller supplies `g_fk`/`fk_mc`/`fk_se`;
  `_compute` still does not produce them, so the behaviour of the generator run
  on its own is unchanged. The "NO FK MARKERS" docstring block is rewritten to
  record both the withdrawal and the restoration.
* `sections/insights.tex`, `sections/conclusion.tex`, `sections/appendix.tex`
  -- the Validation paragraph, the conclusion's numerical-side sentence, and the
  Workflow-validation appendix plus its figure caption. See
  `code/mc_fk_complete/PAPER_PROPOSAL.md` for the before/after of each block.

Figure rebuilt by

    code/figures_corrected/make_val_figure.py \
        --fk-xi products/table_permclosed_cut15360_permfix_xi.npz \
        --fk-markers code/mc_fk_complete/_markers_pooled.npz

(`--fk-markers` is new; omitting it reproduces the FF-markers-only version.)
Deployed to `figures/appendix_mc_workflow.pdf`; the previous version is archived
as `figures/archived_2026-08-28/appendix_mc_workflow_preFKmarkers.pdf`.

Marker values: four independent 24-seed blocks, `n_lambda = 1000`,
`n_real = 24000`, paired-seed linear `sigma_lambda -> 0` extrapolation from
`sigma_lambda = 8` and `4`; error bars are the block-to-block scatter, since
independent blocks scatter by 1.75x the within-block seed error for reasons not
yet understood. Ratios to the fold: 1.012 +/- 0.005 (1'), 1.012 +/- 0.008
(2.6'), 1.008 +/- 0.014 (6.7'), 0.999 +/- 0.033 (17.3').

To revert: `git checkout -- <the .tex files> figures/appendix_mc_workflow.pdf`
and `git checkout -- .../fig_xi_channels.py`.

## 10. FK quoted over 2'-12' instead of at 0.5' (2026-08-28)

`0.5'` was the paper's quoted separation and is also the endpoint of the gamma
grid. It was used for cosmology by exactly one cosmic-shear analysis, KiDS-1000
v1 (Asgari et al. 2021), which measured `xi_+` over `[0.5', 300']` with no lower
cut; the KiDS team retracted it (Li et al. 2023 moved to `theta_min = 2'` on
baryon-feedback grounds, a 0.7-0.8 sigma shift in `S_8`), and KiDS-Legacy
(Wright et al. 2025) keeps `2'`, reporting that a re-test at `0.5'` fails their
B-mode null test. DES Y3 cuts `xi_+` at ~2.5', DES Y6 tighter, HSC Y3 at 7.1',
UNIONS at 12'.

FK/Order-0 in `xi_kappa` over 2'-12' is 1.3-1.7% at `z_s=5` and 0.9-1.1% at
`z_s=1`, from `analysis3/outputs/multiz_kappa_2pcf_5z.npz` and the converged
`cut15360` fold. (Beware: `products/table_permclosed_cut1000_permfixmultiz_xi.npz`
is an `ell_max=1000` product and is a factor ~5 low; it is not what the deployed
multi-z figure plots.)

NOT a claim that the tree-level bispectrum is valid at 2'-12'. It is not valid
at any separation the paper quotes: tree-level SPT for the bispectrum breaks at
`k ~ 0.1 h/Mpc`, a few degrees in this geometry. The argument for 2'-12' is that
it is the range the data are used over and that it avoids quoting at a grid
endpoint. The existing "underestimates the nonlinear small-scale clustering"
caveat carries the validity statement and was not weakened.

One generator edited outside this folder:
`SFT-lensing-paper-analyses/sachs_sft/analyses/analysis3/plot_analysis3_nlo_decomposition.py`
gained `QUOTE_RANGE_ARCMIN = (2.0, 12.0)` and one `axvspan` per panel, shading
the quoted range. Rebuilt with

    code/figures_corrected/regenerate.py --figure nlo \
        --fk-npz ../../products/table_permclosed_cut15360_permfix_xi.npz

and deployed to `figures/analysis3_NLO_FFFK.pdf`; previous version archived as
`figures/archived_2026-08-28/analysis3_NLO_FFFK_preShading.pdf`.

`0.5'` deliberately KEPT in the cutoff subsection: it is the most demanding case
for convergence, so demonstrating it there is stronger than at 12'. A clause now
says so.

To revert: `git checkout -- sections/ figures/analysis3_NLO_FFFK.pdf` and
`git checkout -- <the generator>`.
