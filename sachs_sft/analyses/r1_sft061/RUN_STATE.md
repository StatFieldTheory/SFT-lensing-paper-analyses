# Revision execution state

## Completion on 2026-09-29

The numerical revision, eight affected figures, 34 dependent prose replacements,
and the separately approved one-sentence abstract change are complete. All 12
products are selected in `reproduce/active_products.json`. The final number and
figure checks passed, and the manuscript compiled to 33 pages. Rendered pages
were inspected. See `referee-r1/response_letter_r1.tex` in the manuscript root for Comment
14 and the validation limits, and `final_numbers.md` and `final_figure_check.txt`
in the revision run folder for numerical evidence. This completion record
predates the repository closeout commit.
Earlier pending-state entries below are historical execution records. The
2026-09-29 approval supersedes the earlier abstract restriction for one sentence.


Updated 2026-09-28, approximately 22:05 UTC. This records execution,
not scientific acceptance or manuscript deployment.

## Author decisions

- Actual source redshift 5, with corresponding Order-0 and FF recomputation.
- Current sft-wick 0.6.1 and corrected pair-phase, n_phi=512 K inputs.
- Minimal consistent figure and text revision. Keep the abstract unchanged.
- Preserve historical products. Deploy only through `reproduce/deploy.py`.

## Active and queued work

**Scientific hold discovered after this queue was dispatched:** the restricted
52-row cutoff contract reproduces the historical geometry-sorting callable,
which can misassign spin-dependent cumulants. The scalar convergence FK
contraction also contains Ricci-Weyl modulus terms, so this is relevant to
the cutoff ladder. It does not affect the permutation-aware main or multiz
FK, FF, or the updated MC. The future tree sorting batch is blocked by
archiving its active plan as `prepared_blocked_sorting_callable.json`, before
any target started. The future BiHalofit sorting batch directory is reserved
with `SCIENTIFIC_HOLD.md`, making its queued preparation fail the fresh-folder
check. Both coordinator sessions are therefore expected to stop at these
holds. Do not restore the old plan or clear the reserved directory to resume.
The valid band integrations continue and will be reused where applicable.
A permutation-aware replacement contract and fresh batches are being designed.

1. Main corrected FK completed in 1235 seconds, batch
   `true_redshift_fk_gl24_corrected_main`, four workers, GL24.
   Result SHA256 starts `892e1eb5e4f2`. At 0.5 arcmin and ell_max=15360,
   convergence FK is `1.61873463520648e-5`, or 1.91235% of accepted Order-0.
2. Serial framework continuation: session 97610,
   `revision_fold_queue.log`. Currently runs `true_redshift_ff_gl24_all`,
   followed by corrected multiz FK and corrected
   tree52 ladder, in that order. Never start a manual overlapping fold.
3. Serial BiHalofit build: session 34512,
   `corrected_kappa3_ladder52/queue.log`. Four bands complete, fifth active.
   Eight bands total. No second vertex build may run concurrently.
   Root session 65800 waits for all eight successful bands, assembles and
   validates the five BiHalofit cutoffs, then waits for every tree fold and
   the shared lock before starting the BiHalofit folds. Its log is
   `bihalofit_followup_queue.log`. Do not launch a duplicate continuation.
4. FF baseline `mc_update/ff_n4000` complete, all four angles.
   The selected N8000 check completed after recovery of the last two seeds.
   Its previous agent-owned execution had stopped after six seeds. The
   recovery preserved those seeds and is recorded in `recovery.json`.
   The first dense reference stopped at angle 21.8774 arcmin because a
   pointwise relative tolerance was inappropriate at covariance zero
   crossings. A reviewed component-peak diagnostic passed all 40 angles.
   Both reference grids for all previous 16 angles remained bit-identical.
   That partial directory and log are archived. The root-owned dense
   reference completed all 40 angles. The half-F coupling check also
   completed, with a 0.1948% difference in the normalized signal at 1 arcmin.
5. Root-controlled FK MC queue session 4207 now runs the first of four
   independent blocks, followed by pooling, source-grid quadrature
   and cache assembly serially. Log: `mc_update/fk_mc_queue.log`.

Each queued process has a predecessor timeout and fails rather than replacing
an existing artifact. Check its log and session before any resume.

## Completed gates

- Main z5 Order-0/PyCCL comparison and all five lower-source Order-0 products.
- True-z5 FF endpoint quadrature probe: GL24 to GL48 changes convergence
  by at most 0.1975% at the two tested endpoints. Full batch uses GL24.
- Restricted tree geometry: all tensor values at the contracted queries
  agree bitwise with the full geometry, including scalar and batch calls.
- General and Python reviews passed for restricted fold guards, assembly,
  MC wrappers, coupling probe, product routing and plot adapters.
- The FF covariance diagnostic correction passed independent general and
  Python review. Cache compatibility for its exact recorded runner-hash
  transition passed both reviews as well.
- `deployment_archive_2026-09-28` preserves 52 pre-deployment files,
  including all 17 live PDFs, affected generator outputs, numerical records,
  and the current manuscript sources and PDF. The manifest records hashes.
- Assembler's local type-name warning was corrected to `row_key`.

## Pending

- Finish numerical queues. Assemble and validate BiHalofit tables, prepare
  their permutation-aware replacement batch, and run it after the framework
  queue releases the shared lock. The original 52-row sorting batches are
  not final scientific products.
- Assemble multiz and ladder products with `assemble_figure_products.py`.
- Audit measured FK polarization residuals against the exact scalar FK BB
  null. Do not interpret the old or new residual as physical FK B power.
- Create the complete hash-pinned active product manifest only after all
  final products exist. It is not currently present.
- Archive previous generated/deployed figures and number records. Render,
  inspect, deploy, update affected prose and the referee response, compile.
- Flag the abstract's comparable FK E/B claim without editing it.

## Resume on 2026-09-29

The earlier session list above is historical. Main FK, all FF source planes,
FK multiz, four FK MC blocks, the pooled markers, refined quadrature, and
MC figure cache are complete. The multiz products have been assembled.
The original sorting-callable cutoff plans remain on scientific hold.

The replacement uses exactly 34 ordered input rows with the original
permutation-aware callable. The full and restricted tensors are bitwise equal
for 2728 scalar and 2728 batch comparisons. Unsupported direction, mixed-angle,
and nonfinite-time queries are rejected. Evidence: `permaware_ladder/guard_validation.json`.
General and Python review approved the new preparation, merge, and guard files.

Current root-owned jobs:

- Session 76042: serial tree folds, `permaware_ladder/folds_tree/run.log`.
- Session 57580: serial eight-band BiHalofit build of the 22 missing ordered
  rows, `permaware_ladder/queue.log`. The other 12 required rows are reused
  exactly from validated existing pieces. Do not run another vertex build.

After both jobs complete, merge the BiHalofit pieces, prepare and run its
five folds serially, and assemble the ladder. The assembler requires bitwise
agreement of the highest tree cutoff with the 11 corresponding full-main
convergence values, and verifies physical and quadrature setup equality.
`active_products_candidate.json` pins the eleven completed products but is
not active and lacks the final cutoff directory. Do not activate it before
the remaining gates pass.

The two plot-comment corrections on this resume acknowledge that PCHIP is
applied separately to each total curve. The quadrature matrix is shared, but
the complete finite-grid interpolation is not a linear operator. Detailed
harmonic amplitude claims are omitted from the minimal prose revision.

### Sequential finish coordinator

Tree session 76042 exited successfully. The real 11-angle highest-cutoff
comparison is bitwise identical to the full main fold, recorded in
`permaware_ladder/main_fold_validation.json`.
Root session 17835 now waits for the successful eighth vertex band, then runs
`pieces.py merge-bihalofit`, `suite.py bihalofit`, the serial BiHalofit fold
batch, and the guarded final ladder assembly. Its log is
`permaware_ladder/bihalofit_finish.log`. Do not manually duplicate these steps
while that coordinator is active. It stops on any failure and will not replace
existing outputs. Its vertex wait is bounded to one hour.

### Final reviews and author-approved abstract change

The code-reviewer approved the final assembler and product-resolution checks,
including all 11 completed candidate-product hashes. The python-reviewer approved
`check_numbers.py` and the final CL/MC plot changes, with targeted Ruff
E9,F63,F7,F82 passing. The assembler rejects mismatched runtime/config quadrature
orders. No active cutoff numerical claim is accepted by these code reviews alone.

The staged manuscript replacements are in the parent manuscript's ignored
`referee-r1/numerical-text-replacements.json`. They have unique, nonoverlapping
old strings. The final scientific review corrected an adjacent 3PCF scope issue:
the direct K projection is not the whole leading Order-1 connected 3PCF, because
F contributes at that order too. Master formulas are unchanged.

The author explicitly approved replacing the single abstract E/B sentence with:
"The two-point function, for instance, receives corrections from three-point
cumulants of Ricci focusing and Weyl shearing through nonlinear Sachs evolution."
This edit is implemented in main.tex and recorded in response.md. It supersedes
the earlier abstract hold for this one sentence only. Vector/tensor sources can
have linear B, so the scalar FK BB condition is confined to the worked example.

### Final-fold recovery

All eight missing-row bands and their merge completed successfully. The first
finish coordinator correctly stopped before writing the fold batch because its
preparation step used the PyCCL interpreter, whose installed distribution
metadata says sft-wick 0.1.0 despite the shared source being at the pinned commit.
No fold ran under that interpreter. Root resumed with the designated sft-wick
interpreter in session 17679. It prepares and runs the five BiHalofit folds,
then assembles the final ladder. The existing completed merge is preserved.
The same `bihalofit_finish.log` records both attempts.
