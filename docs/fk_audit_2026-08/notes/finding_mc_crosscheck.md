# MC cross-check: it validates the fold, and it proves it cannot judge the input

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Run 2026-08-25 with `code/mc_crosscheck_dense.py`, which points
`analyses/mc_sachs_2pt/fk_analytic.py` at a chosen kappa3 table. That module
re-assembles the FK two-point diagram with its own quadrature: its own
causal simplex, its own response product, its own F-K index contraction.

## Result on the dense-grid table

| gamma ['] | independent assembly | production fold | ratio |
|---|---|---|---|
| 0.500 | 4.27573e-06 | 4.28149e-06 | 0.9987 |
| 3.307 | 4.05544e-06 | 4.06085e-06 | 0.9987 |
| 8.506 | 3.00002e-06 | 3.00373e-06 | 0.9988 |
| 21.877 | 8.72784e-07 | 8.73386e-07 | 0.9993 |
| 56.267 | 2.22925e-07 | 2.22997e-07 | 0.9997 |
| 144.713 | 1.15744e-08 | 1.15421e-08 | 1.0028 |
| 372.190 | -3.30800e-09 | -3.31202e-09 | 0.9988 |

Median ratio 0.9988, spread 0.0012, range [0.9971, 1.0038]. The agreement
holds across three orders of magnitude of decay and through the sign
change. So the production fold handles the densified table correctly: the
steep angular decay is not an artefact of the fold, it survives an
independent assembly of the same diagram.

## The control, and what it proves

The same independent assembly, pointed at the DEPLOYED table:

| gamma ['] | independent | production | ratio |
|---|---|---|---|
| 0.500 | 4.28195e-06 | 4.28772e-06 | 0.9987 |
| 8.506 | 4.28173e-06 | 4.28750e-06 | 0.9987 |
| 144.713 | 4.21902e-06 | 4.22471e-06 | 0.9987 |
| 372.190 | 3.88811e-06 | 3.89335e-06 | 0.9987 |

Median 0.9987, spread 0.0000. Equally good agreement, on a curve that is
FLAT instead of decaying.

That is the point. The cross-check reproduces whatever the table says, to
one part in a thousand, whether the table is frozen at its corner or
resolved. It validates the diagram assembly and the fold; it has no
opinion about the input statistics, because both sides read the same
input.

This settles, with a direct measurement, the reading recorded in
`impact_on_draft_discussion.md`: the appendix Monte-Carlo validation of the
FK channel does not certify the FK curve's angular shape. A Monte-Carlo
seeded from either table agrees just as well, at completely different
values.

## Scope

* Restricted to gamma < 700 arcmin; beyond that both sides sit in the
  cancellation floor where ratios are uninformative.
* The 0.13% offset is the same constant the original anchor reports
  (it quotes 0.992 against the older analysis-3 output), and comes from
  sft-wick's coarser internal discretisation, not from either table.
