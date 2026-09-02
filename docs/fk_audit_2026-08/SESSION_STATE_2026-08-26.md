# Session state, 2026-08-26

Written before a context compaction. Everything needed to resume.

## Deliverables that exist

| what | where |
|---|---|
| audit note, 28 pages, 12 sections | `note/main.pdf` (build: `latexmk -pdf main.tex`) |
| email to editor/referee, ~335 words | `revision_plan/letter_email.txt` (plain text, ready to paste) |
| long letter (superseded by the email) | `revision_plan/letter.tex` |
| revision plan, 8 pages, 12 numbered text changes | `revision_plan/plan.pdf` |
| corrected callable + 7 tests | `code/callable_fixed/` |
| factored rebuild machinery | `code/rebuild/` (see its README) |
| xAct beam-width derivation, 15/15 pass | `../SFT-lensing-paper-analyses/mathematica/beam_width_regulator.wl` |

## The final numbers (ell_max = 15360, permutation-aware callable)

Table: `products/table_permclosed_cut15360.npz`
Fold:  `products/table_permclosed_cut15360_permfix_xi.npz`

| observable | peak | at | % of Order-0 |
|---|---|---|---|
| xi_kappa | +1.948e-05 | 0.5' | 2.31 |
| xi_plus | +1.959e-05 | 0.5' | 2.33 |
| xi_minus | +5.687e-07 | 5.3' | 0.93 |
| xi_kappa_gamma_t | -1.628e-06 | 2.6' | 1.63 |

* no crossover with Order-0 at any cutoff
* B/E = 0.47 (max|dC_EE| = 3.69e-06, max|dC_BB| = 1.73e-06)
* |E-B|/|E+B| reaches 1 at gamma = 31'
* gamma^4 law survives only below 1': slope 3.94 at 0.5', 2.51 at 1.3', <1 above 2'
* DO NOT quote the -9.7% for xi_kappa_gamma at 0.5': both O0 and FK vanish as
  gamma -> 0 with slopes 1.55 and 0.44, so it is a quotient of two small
  numbers. Quote the plateau, about -1.2% over 3'-5', easing to -0.7% by 70'.

## Draft changes made so far (all highlighted)

`main.tex` gained a highlight macro. NOTE: REVTeX already defines `\revised`,
so the macro is `\edited{...}`. Set `\revhlfalse` for a clean copy.

1. **Figure 2** `analysis3_NLO_FFFK.pdf` replaced. Old -> `figures/archived_2026-08-26/analysis3_NLO_FFFK_preRevision.pdf` (git mv).
2. **Figure 3** `analysis3_cl_full_vs_O0.pdf` replaced, old archived likewise.
3. **Figure 3 caption** in `sections/insights.tex`: added, inside `\edited{}`,
   that structure below ell ~ 20 is finite-range and apodisation dominated and
   quantitative statements are restricted to ell >= 50.
4. **Figure 4** `cl_EB_polarization.pdf` replaced, old archived.
5. **Generator edit** (the only one): `analyses/analysis3/plot_cl_EB_polarization.py`,
   the EB annotation "$\Delta C_\ell^{EB}=0$ (exact)" -> "by parity", with eight
   lines of comment explaining that it is a theorem the code enforces
   structurally, not a measurement. MUST be recorded in
   `CHANGES_OUTSIDE_THIS_FOLDER.md`.

Draft builds: 31 pages, `latexmk -pdf main.tex` from the repo root.

## Figure 5: DONE (2026-08-26 afternoon)

The first regeneration attempt swept all three channels and stalled in the FF
loky sweep (killed at 1h47m, the known failure mode). `regenerate_multiz.py`
gained `--reuse`: O0 and FF are served bit-identical from the June npz (neither
touches the three-point vertex), only FK is recomputed (1034 s). FK/O0 at 0.5'
by z_s: 1.36 / 1.77 / 2.04 / 2.14 / 2.24 % (was flat ~0.53% pre-fix), monotone
toward the z=5 value 2.31%. Figure rebuilt with the paper's own
`plot_multiz_kappa_2x5.py`, deployed to `figures/multiz_kappa_xi_cl.pdf`; old
version archived as `figures/archived_2026-08-26/multiz_kappa_xi_cl_preRevision.pdf`
(git mv). Awaiting user approval.

NOTE the corrected figure still shows a flat Full-line tail at gamma > 300':
that is now the FF disconnected mean-field^2 piece (see the convention finding
below), not the FK artifact.

## Figure 6: CLOSED (2026-08-26, option A executed)

The FK Monte-Carlo discrepancy was fully explained by an exact second-order
expectation analysis of `simulate_fk_vr` (3-agent workflow audit; 21 probe
points matched to ~1.5 seed-SEM): the estimator measures a placement SHARE of
the vertex (~0.245 at 0.5', zero near 10', negative beyond) plus a
finite-sigma_lambda term from lambda-grid graininess of the tabulated tensors;
the old sigma_lambda=8 "agreement" was the two terms summing through unity at
the tuned value. User chose option A: FK markers removed; FK validated by
`code/rebuild/fk_kernel_crosscheck.py` (direct quadrature vs converged fold:
0.99-1.02 over 0.5'-5', 1.07 at 17.3'). Generator `fig_xi_channels.py` edited
(CHANGES_OUTSIDE_THIS_FOLDER.md #7); figure rebuilt with FF markers at the
original 15-point grid, deployed, June original archived.


## Convention check: RESOLVED, no action

The user asked whether the MC shares the paper's connected-cumulant
convention. Checked all four layers: vertex table, Q injection, FK estimator
(its disconnected subtraction is variance isolation for the FK diagram, not a
convention), FF channel. Everything is mutually consistent. The paper's
convention is: DRIVING-FIELD statistics are cumulants; OBSERVABLES (xi) are
raw moments <kappa kappa>. The FF flat tail at large gamma is the
second-order mean-field square <k>_Fon^2 (~5.8e-7) and is CORRECT under that
definition -- it is physics, not an artifact, and needs no paper change. (I
briefly misread it as a convention inconsistency; the user corrected me.)

## Plan status: ALL text changes applied (2026-08-26 afternoon)

Changes 1-8, 11, 12, new item 9 (subsection "Convergence in the multipole
cutoff" in insights.tex, label subsec: cutoff, + figures/fk_cutoff_convergence.pdf
via code/rebuild/fig_cutoff_paper.py) and item 10 (UV sentences in the
conclusion), plus the four positive statements A-E, all applied with
\edited{} highlighting; build green, 31 pp. Validation reframed everywhere
(appendix two paragraphs + caption, insights Validation paragraph, conclusion
sentence): FF by Monte-Carlo, FK by the deterministic kernel cross-check.
biblio.bib gained carrasco2012effective. Email numbers already converged
(2.3%, half, 31'). A final 3-agent review workflow was launched at the end;
apply its findings before shipping.


## Two decisions still open

* the headline number: quote 2.31% at the converged cutoff, or keep 0.5% with
  a caveat. The plan recommends the former.
* the email says a revised version "has been posted to arXiv today". That is
  written ahead of the fact; do not send until the revision is actually made.

## Machine

103 GB, shared with another session. An out-of-memory reboot happened once
today. Chunk canoes builds over rows and cap concurrency; see
`code/rebuild/README.md`. Folds: `run_fk_variant.py --n-jobs` is a memory
control, never run two folds at once.

## Final review: DONE (2026-08-26, end of day)

A 3-agent review workflow (numbers audit / stale-claim sweep / LaTeX-style)
checked all 24 edited blocks and every section in full. All findings fixed:
three stale unedited captions rewritten (C_ell figure caption in insights,
EB caption and multi-z caption in conclusion), the insights closing paragraph
("FF the only correction" -> FK leads all four), the dangling "angular
crossover" clause removed, intro's "conservative" replaced by the
smoothing-scale statement, two unqualified gamma^4 statements qualified,
the z=5-vs-figure-range citation recast (figure spans z<=4, 2.2%; z=5 = 2.3%
in the main analysis), the reduced-shear sentence pinned to "matched
multipole cutoff", "third cumulant" -> "three-point cumulant", and the
fragile \edited{\subsection...} restructured. Verified-correct and left
alone: FF-overtakes-O0 at ~200' (recomputed: ~190'), the EB caption's
"eight orders" residual (matches the regenerated figure's annotation).
Build green, 31 pp, no undefined refs. Not committed; Overleaf sync happens
on master push per the usual convention.
