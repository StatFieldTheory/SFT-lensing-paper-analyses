# Proposed revision: quote FK over 2'-12' rather than at 0.5'

*2026-08-28. Nothing applied. Four text blocks, one caption clause, one optional
figure change.*

## Why

`0.5'` was used for cosmology by exactly one cosmic-shear analysis, KiDS-1000
v1 (`asgari2021kids`), which measured `xi_+` over `[0.5', 300']` with no lower
cut. The KiDS team retracted it: Li et al. 2023 moved to `theta_min = 2'` on
baryon-feedback grounds, a 0.7-0.8 sigma shift in `S_8`, and KiDS-Legacy
(Wright et al. 2025) keeps `2'`. When KiDS-Legacy re-tested `0.5'-300'` the
cosmology moved only ~0.1 sigma but the **B-mode null test failed** their
`p > 0.01` requirement. Elsewhere: DES Y3 (`amon2022dark`) cuts `xi_+` at
~2.5', DES Y6 tighter, HSC Y3 at 7.1', UNIONS at 12'.

`0.5'` is also the endpoint of our own gamma grid, and at the FK kernel's
median distance it is a transverse separation of `0.12 Mpc/h`, `k ~ 26 h/Mpc`.

**A caveat that must not be over-read.** Quoting at 2'-12' is NOT a claim that
the tree-level bispectrum is valid there. It is not valid at any of these
separations: tree-level SPT for the bispectrum breaks at `k ~ 0.1 h/Mpc`, which
in our geometry is `theta` of a few degrees. The argument for 2'-12' is that it
is the range the data are actually used over, and that it avoids quoting at a
grid endpoint. The existing statement that the tree-level bispectrum
"underestimates the nonlinear small-scale clustering" continues to carry the
validity caveat and should not be weakened.

## The numbers

FK as a fraction of Order-0 in `xi_kappa`, from the data the deployed
multi-z figure plots (`analysis3/outputs/multiz_kappa_2pcf_5z.npz`) and the
converged `cut15360` fold for the main plane:

| `z_s` | 0.5' | 2' | 5' | 8' | 12' | **over 2'-12'** |
|---|---|---|---|---|---|---|
| 1.0 | 1.36 | 1.13 | 1.02 | 0.93 | 0.86 | **0.9-1.1%** |
| 1.7 | 1.77 | 1.41 | 1.28 | 1.16 | 1.07 | 1.1-1.4% |
| 2.5 | 2.04 | 1.58 | 1.43 | 1.30 | 1.19 | 1.2-1.6% |
| 3.2 | 2.14 | 1.64 | 1.49 | 1.35 | 1.23 | 1.2-1.6% |
| 4.0 | 2.24 | 1.70 | 1.55 | 1.40 | 1.28 | **1.3-1.7%** |
| 5.0 (main) | 2.31 | 1.74 | 1.58 | 1.43 | 1.31 | **1.3-1.7%** |

## 1. `sections/intro.tex` (l. 111)

### Current
```latex
acquires a correction of $2.3\%$ of the linear two-point signal at
$\gamma=0.5'$, growing with source redshift from $1.4\%$ at $z_s=1$
```
### Proposed
```latex
acquires a correction of $1.3$--$1.7\%$ of the linear two-point signal across
$2'$--$12'$, the range current cosmic-shear analyses
retain~\citep{asgari2021kids,amon2022dark}, growing with source redshift from
$0.9$--$1.1\%$ at $z_s=1$
```

## 2. `sections/conclusion.tex` (ll. 34-41)

Two changes: the numbers, and the phrase "the scales where cosmic shear is
measured best". Arcminute scales carry the most pairs, but they are precisely
the scales current analyses **discard**, so as written the clause invites the
objection it is trying to pre-empt.

### Current
```latex
\edited{On a source plane at
$z_{s}=5$ this leakage lifts the convergence and the parity-even shear by about
two percent at arcminute separations, the scales where cosmic shear is measured
best, and remains a small fraction of the linear signal at every separation
(Section~\ref{subsec: cutoff}); it strengthens with source redshift, from
$1.4\%$ at $z_{s}=1$ to $2.2\%$ at $z_{s}=4$ at $\gamma=0.5'$
(Fig.~\ref{fig: multiz kappa}), reaching $2.3\%$ on the $z_{s}=5$ plane of
the main analysis, and it leads the nonlinear-propagation term in
all four two-point combinations.
```
### Proposed
```latex
\edited{On a source plane at
$z_{s}=5$ this leakage lifts the convergence and the parity-even shear by
$1.3$--$1.7\%$ across $2'$--$12'$, the range current cosmic-shear analyses
retain after their small-scale
cuts~\citep{asgari2021kids,amon2022dark}, and remains a small fraction of
the linear signal at every separation
(Section~\ref{subsec: cutoff}); it strengthens with source redshift, from
$0.9$--$1.1\%$ at $z_{s}=1$ to $1.3$--$1.7\%$ at $z_{s}=4$
(Fig.~\ref{fig: multiz kappa}), and it leads the nonlinear-propagation term in
all four two-point combinations.
```
The `2.3\%` on the `z_s=5` plane is dropped because at 2'-12' the `z_s=4` and
`z_s=5` ranges coincide to two figures, so the clause no longer adds anything.

## 3. `sections/insights.tex` (ll. 213-216)

### Current
```latex
\edited{There, at small separation, FK reaches $2\times10^{-5}$, about
$2.3\%$ of Order-0 and an order of magnitude above the FF term; the ratio
declines gently with separation, to below one percent by a degree, the two
terms decaying together.}
```
### Proposed
```latex
\edited{There FK reaches $1.2\times10^{-5}$ at $2'$, $1.3$--$1.7\%$ of Order-0
across $2'$--$12'$ and an order of magnitude above the FF term; the ratio
declines gently with separation, to below one percent by a degree, the two
terms decaying together.}
```

## 4. `sections/insights.tex` (l. 242), one clause

`0.5'` should STAY in the cutoff subsection: it is the most demanding case for
convergence, so demonstrating it there is stronger than at 12'. But with the
headline moved, a reader meeting `0.5'` here will wonder why. One clause fixes
it.

### Current
```latex
$\xi_{\kappa}$ at $\gamma=0.5'$ is $0.48\%$, $0.97\%$, $1.53\%$,
```
### Proposed
```latex
$\xi_{\kappa}$ at $\gamma=0.5'$, the smallest separation computed and so the
most demanding for convergence, is $0.48\%$, $0.97\%$, $1.53\%$,
```

## 5. Optional: shade the quoted range on Fig. `fig: NLO 2pt FFFK`

A light band over `2'-12'` in the four panels would make the quoted range
visible rather than asserted, and would show at a glance that FK leads FF
throughout it. Generator
`analyses/analysis3/plot_analysis3_2pt_decomposition.py`; one `axvspan` per
panel plus a caption clause. Say if you want it and I will produce it for
review.

## Not changed

* `insights.tex` l. 66, `gamma in [0.5', 2000']` -- that is the computed range,
  not a quoted result.
* Fig. `fig: fk cutoff` caption -- the convergence figure is at `0.5'` for the
  reason given in item 4.
* The `xi_-` and `xi_kappa_gamma` numbers (`0.9%` near `5'`, `1.2%` over
  `3'-5'`) are already quoted inside the retained range.
