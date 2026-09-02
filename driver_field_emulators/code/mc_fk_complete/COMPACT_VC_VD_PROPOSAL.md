# Proposed compaction of Sections V C (end) and V D

*Nothing applied. Word counts exclude the figure environments.*

| | now | proposed |
|---|---|---|
| V C, spin-structure paragraph | 335 words, ~18 numeric tokens | **~135 words, ~9** |
| V D, cutoff subsection | 625 words, 6 paragraphs | **~380 words, 3 paragraphs** |

---

## V C: the fat is three restatements of the same split

The paragraph states the modulus-vs-other split **four times**: once as the spin
mechanism, once with real-space numbers, once as a summary sentence, and once
again in harmonic space. The last two add nothing the first two did not, and
the harmonic-space sentence is a near-verbatim duplicate of the caption of
Fig.~\ref{fig: NLO 2pt FFFK cl}, which already says "a slowly rising
correction, from about one percent to just above two percent, in
$C_\ell^{\kappa\kappa}$ and $C_\ell^{EE}+C_\ell^{BB}$ (top row), entering
$C_\ell^{EE}-C_\ell^{BB}$ and $C_\ell^{\kappa E}$ (bottom row) at a reduced but
nonzero level."

### Proposed

```latex
How large FK is depends on each observable's spin structure. The shear is a
spin-$2$ field, so a pair of shear legs can combine in two ways. Paired into a
modulus they make a net spin-$0$ quantity that stays finite as $\gamma\to0$ and
tracks the convergence; in their only other parity-even pairing the spins add
rather than cancel, and the result is washed out as $\gamma^{4}$. FK therefore
concentrates in the two 2PCFs built from the modulus, the convergence
$\xi_{\kappa}$ and the parity-even shear $\xi_{+}$, which are themselves nearly
equal at small angles: \edited{there it runs $1.3$--$1.7\%$ of Order-0 across
$2'$--$12'$, an order of magnitude above FF. In the shear difference $\xi_{-}$
and the convergence-shear cross $\xi_{\kappa\gamma_{t}}$ it is suppressed
rather than removed, to about a percent, the $\gamma^{4}$ law setting in only
below an arcminute; it stays above FF in all four. The same split carries over
to harmonic space over $50\lesssim\ell\lesssim1500$
(Fig.~\ref{fig: NLO 2pt FFFK cl}).}
```

### Cut

* `FK reaches $1.2\times10^{-5}$ at $2'$` -- an absolute value the figure shows.
* `the ratio declines gently with separation, to below one percent by a degree,
  the two terms decaying together` -- a trend the figure shows.
* `it reaches $0.9\%$ of Order-0 in $\xi_{-}$ near $5'$ and about $1.2\%$ in
  $\xi_{\kappa\gamma_{t}}$ over $3'$--$5'$` -- becomes "about a percent". Four
  numbers and two separations for a statement whose content is "about a percent".
* `Driving-field non-Gaussianity therefore does not contaminate every 2PCF
  equally, entering the modulus pair at around one and a half percent and the
  remaining pair at around one` -- a restatement of the two preceding sentences.
* The harmonic-space numbers -- kept in the caption, where they already are.
  Only the $\ell$ range survives in the text, because the caption lacks it.

---

## V D: six paragraphs to three

The subsection argues one thing: FK is a coincident-point moment, so the cutoff
is physics rather than bookkeeping; at tree level the sum converges and we show
it; with a nonlinear input it does not, and that is structural. Paragraphs 5, 6
and 7 all end on the same conclusion.

### Proposed structure

1. **Motivation and the tree-level convergence** (paras 2+3 merged, ~105 words).
   Unchanged in content; the five percentages become endpoints plus the four
   ratios, which are the actual evidence of convergence.
2. **Why the quoted range is limited** (para 4, ~95 words from 145). The
   $P_{\ell_3}(\cos\gamma)$ oscillation argument, trimmed.
3. **The nonlinear case and why it is structural** (paras 5+6+7 merged, ~180
   words from 360).

### Cut

* `The two are halves of one prescription, not alternatives.` -- a summary of
  the sentence before it.
* `That the two-point correction is finite at tree level is a property of the
  linear spectrum falling steeply enough, not a general feature.` -- duplicates
  the $n_{\rm eff}$ paragraph.
* `This ultraviolet sensitivity is a structural consequence of the formalism
  rather than an accident of the example.` -- a topic sentence that the very
  next sentence proves; folded in.
* `$n_{\rm eff}=-1.15$ at $k=1\,h/\mathrm{Mpc}$, the peak of the integrand` --
  the $0.21$ and $8.4\,h/\mathrm{Mpc}$ crossings carry the argument.
* `because the nonlinear Sachs terms do not commute with the smoothing` -- the
  reason smoothing changes the propagation; the fact is kept, the mechanism goes.

### KEPT deliberately, because they were added on purpose in earlier rounds

* the smoothing scale `$\ell_{\rm eff}\simeq2\chi/s_{b}$`;
* the closing caveat that a realistic evaluation would also have to carry the
  white-noise contribution to the small-scale variance, which the tree-level
  spectrum does not;
* `$\gamma=0.5'$, the smallest separation computed and so the most demanding
  for convergence`.

Say if any of the "cut" items should stay, or if V D should go further.
