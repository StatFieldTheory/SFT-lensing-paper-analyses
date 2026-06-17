# TASK 1 — the no-K (FF) MC-vs-analysis-3 discrepancy: VERDICT (2026-06-06)

Tool: `ff_channel_check.py` (F-on/F-off CRN, `sachs_mc_core.simulate_ff_crn`).
Figure: `figures/ff_channel_mc_vs_analysis3.pdf`.

## Verdict: the MC setup is CORRECT. No bug. The apparent discrepancy is two
## separate, fully-explained effects.

### (b) The large-gamma "discrepancy" is the DISCONNECTED <kappa>_2^2 piece.
analysis-3's FF flattens to a gamma-INDEPENDENT floor ~5.785e-7 for gamma > ~150',
while O0 decays and crosses zero (~200'). That floor is NOT a correlation: it is
the DISCONNECTED square of the mean second-order convergence. The F-vertex makes
`<kappa>_2 = int W F <s s> dlam != 0`, so as the two rays decorrelate
(`<f(n1)f(n2)> -> 0`), `<kappa(n1)kappa(n2)> -> <kappa>_2 <kappa>_2` = const.

CODE-LEVEL CONFIRMATION (2026-06-06 -- this was challenged as possibly "fitting an
explanation"; it is NOT a fit. analysis-3's FF is the full MOMENT, NOT connected-only):
1. DIAGRAM ENUMERATION (sft-wick `compute_moment([phi_a(x),phi_b(y)], order=2,
   vertex F)`): the FF order-2 moment has 6 diagrams = 5 connected + 1 DISCONNECTED.
   The disconnected one is literally
   `R(x,y0) R(y,y1) C(y0,y0) C(y1,y1)`  (x and y in separate components)
   = `<kappa(x)>_2 x <kappa(y)>_2` (each obs -> its own F-vertex with a C-tadpole).
   So sft-wick DOES include it (no connected-only filter; `is_connected` is used only
   for rendering).
2. ITO keeps it: the Ito prescription kills equal-point RESPONSE loops (R(x,x)=0) but
   the equal-point CORRELATION tadpole C(v,v) is always returned -> <kappa>_2 != 0.
3. VALUE MATCH with the EXACT analysis-3 inputs (NOT the MC's de-windowed Sigma2):
   `<kappa>_2 = int W(tau) F_0bc C_bc(tau,tau) dtau` with C = the corr_op propagator
   (the same `c_closed_form` analysis-3's FF uses) gives <kappa>_2 = -7.587e-4,
   `<kappa>_2^2 = 5.756e-7` vs the FF floor 5.785e-7 -> ratio 0.995 (0.5%). A
   connected-only FF would vanish as the cross-ray C -> 0, so a gamma-independent
   floor equal to <kappa>_2^2 can ONLY be the disconnected tadpole-squared diagram.

The MC (F-on/off CRN, anchor OFF) reproduces it independently:
* MC disconnected `<kappa>_2^2` (median, all gamma) = 6.11e-7  vs  analysis-3 floor
  5.785e-7  -> ratio 1.057 (5.7% high; the equal-time-density mean is slightly off).
* The MC FF MOMENT tracks analysis-3 FF across ALL gamma (0.92x at 0.5' rising to
  1.02x at large gamma).
* The MC FF CONNECTED piece (moment - disconnected) DECAYS like analysis-3's
  (FF - floor): they overlay from 1.2e-6 (0.5') down to ~1e-9 where both hit noise.

So FF "floors" while O0 "decays/crosses zero" simply because FF carries the
disconnected `<kappa>_2^2` monopole and O0 does not. For a CONNECTED correlation
function xi_kappa this `<kappa>_2^2` should be subtracted (it is the mean-
convergence-squared). This is ALSO the resolution of the original puzzle ("why does
the earlier MC differ so much from analysis-3 full at large gamma"): that MC was
F-OFF (pure O0, no disconnected piece, decays/crosses zero) while analysis-3 full
carries FF's disconnected floor. With F on, the MC floors too. NOT an MC bug; NOT a
sft-wick bug; it is a moment-vs-connected-cumulant definitional point.

### (a) The ~10% small-gamma offset is the equal-time finite-correlation deficit.
The MC FF CONNECTED is ~0.83x analysis-3 (FF-floor) at 0.5', rising to ~0.98x by
90'. This is the same finite-lambda-correlation tail the strict equal-time
(white-noise) collapse drops (the ANCHOR_C0 = 1/0.881 story for O0). IMPORTANT: FF
scales as Sigma2^2 (two F-legs, four f's -> Sigma2^2), while O0 scales as Sigma2^1.
So the equal-time deficit enters FF SQUARED, and a SINGLE ANCHOR_C0 (tuned to make
O0 ~ Sigma2^1 overlay) cannot simultaneously overlay O0 and FF. Anchor ON
overcorrects FF by ~1.135x (so MC FF moment / analysis-3 runs 1.15 .. 1.43 with
anchor ON; with anchor OFF it is 0.92 .. 1.02 -- closer, the honest equal-time
deficit). Conclusion: report FF with anchor OFF and the Sigma2^2 caveat; do not
force a single anchor across channels.

### (c) The CRN resolves the FF signal.
The raw single-arm MC `<kappa kappa>` has a ~6e-6 variance floor at large gamma
(cannot see ~5e-7). The F-on/off CRN difference has SE ~1e-7 (the O0 floor cancels
per realization), so the FF channel and its decomposition ARE resolvable.

## Numbers (anchor OFF, n_real=60k, n_lambda=600)
| gamma | MC FF mom | MC FF conn | MC disc | an3 FF | an3 - floor |
|------:|----------:|-----------:|--------:|-------:|------------:|
| 0.5'  | 1.85e-6   | 1.19e-6    | 6.61e-7 | 2.01e-6| 1.43e-6 |
| 5.3'  | 1.26e-6   | 6.08e-7    | 6.49e-7 | 1.27e-6| 6.91e-7 |
| 21.9' | 7.98e-7   | 1.72e-7    | 6.26e-7 | 7.67e-7| 1.88e-7 |
| 90'   | 6.12e-7   | 1.23e-8    | 6.00e-7 | 5.91e-7| 1.25e-8 |
| >300' | 5.92e-7   | ~0         | 5.92e-7 | 5.785e-7| ~0 |

## Implication for the paper figure
When overlaying the MC on analysis-3's "full" xi_kappa, either (i) compare the
MC (F-on) to analysis-3 full directly (both carry the <kappa>_2^2 floor), or
(ii) subtract the disconnected `<kappa>_2^2` from BOTH to show the physical
connected correlation. Mixing an F-off MC against analysis-3-full (the earlier
plot) is the apples/oranges that produced the large-gamma "discrepancy".
