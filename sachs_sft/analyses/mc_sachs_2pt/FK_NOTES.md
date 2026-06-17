# FK channel — findings (2026-06-06)

Two results this session: (1) the analytic FK anchor `fk_analytic.py` is DONE and
validates the FK expression against analysis-3; (2) the naive deformation-sampling
Path A FK MC is blow-up-free but variance-pathological — diagnosed, with the
variance-reduced fix specified.

## 1. Analytic FK anchor — DONE (`fk_analytic.py`)

`fk_analytic(cos_gamma, lam_f)` reproduces the FULL analysis-3 FK (3,3) matrix to a
flat, gamma-independent ratio **0.992**:

| pair | meaning | fk_analytic | analysis-3 | ratio |
|------|---------|-------------|-----------|-------|
| (0,0) | kk    | 5.417e-5 | 5.458e-5 | 0.9925 |
| (0,1) | k g+  | 4.015e-11 | 4.049e-11 | 0.9916 |
| (0,2) | k gx  | 0 | 0 | parity-zero ✓ |
| (1,1) | g+ g+ | 2.709e-5 | 2.729e-5 | 0.9925 |
| (1,2) | g+ gx | 0 | 0 | parity-zero ✓ |
| (2,2) | gx gx | 2.709e-5 | 2.729e-5 | 0.9925 |

The 0.8% is a STABLE constant under my grid refinement (n_lam 48..288 → 0.9917,
0.9926, 0.9913, 0.9918), so it is NOT my under-resolution; it is sft-wick's own
coarser internal discretization (config: dt=1.0, c_n_gauss=20, 16-shell zeta
lambda-interpolation). The selection rules (FK(0,2)=FK(1,2)=0) and spin-2 isotropy
(FK(1,1)=FK(2,2)) emerge correctly. MSR/combinatorial prefactor = 1 (calibrated on
the validated demo2 `_fk_spatial_integral`).

### The corrected FK expression (the earlier handoff schematic was wrong)
FK 2-pt diagram = 4 R-propagators, 0 C-propagators (4 psi legs → 4 phi; psi-psi=0).
Topology:

    obs a --R--> F-vertex ==R,R==> two K-legs ;   third K-leg --R--> obs b

So the external index a rides the F-VERTEX (F_{a m n}); only b rides ONE zeta leg
(K_{b m n}); F's two phi (m,n) contract the two CO-LOCATED zeta legs. Cosine-triple
= (1, cos gamma, cos gamma) (two legs co-located on ray n_a, one on n_b), NOT a
generic triple. Equal-time collapses demo2's 4-D time integral to 2-D; for
integrate_over='all' each external R is the LOS window W:

    xi_FK_ab = int_0^lam_f dtau int_0^tau dlam_v
        W(tau) W(lam_v) [D(lam_v)/D(tau)]^4  G_ab(lam_v; gamma),
    G_ab = (M + M^T)_ab,  M_ab = sum_{mn} F_{a m n} Z_{b m n}(1,cosg,cosg; lam_v),
    W(la) = int_la^lam_f [D(la)/D(t)]^2 dt.

(For kk: G_00 = -2 (zeta_TTT + zeta_Bmod).)

### What the anchor does / does NOT check
It cross-checks the sft-wick FK DIAGRAM ASSEMBLY (equal-time collapse, causal 2-D
measure, F-K index contraction, (1,cosg,cosg) geometry, prefactor) from a
from-scratch quadrature. It does NOT independently validate the zeta table's
internal UNITS — both fk_analytic and analysis-3 read the SAME zeta table, so a
units error there cancels. So the node skewness ~5 is NOT a zeta units inflation;
it is a physically-large tidal skewness, and the OBSERVABLE FK is small only
because it is a perturbative cross-term.

## 2. Path A (perturbative-in-zeta FK MC) — blow-up-free but VARIANCE-PATHOLOGICAL

`simulate_fk_pathA` in `sachs_mc_core.py`: evolves the Gaussian arm `s_g` (full
F-vertex) and a LINEAR field `s_d` driven by the skew increment
`df = 0.5 Q(z x z - V)` through the F-vertex LINEARISED about s_g; FK =
`<kappa_g x kappa_d + kappa_d x kappa_g>` cross-ray. Linearity in s_d removes the
Riccati blow-up (n_good = 40000, no caustics — the structural win).

BUT the estimator does not resolve FK. Diagnostics at gamma=10', sigma_lambda=30:

* fk_analytic target: 5.448e-5.
* Path A: 3.5e-3 .. 5.8e-3 (≈ 70-100x too big), with huge SE ~6e-4.
* F-OFF test (must be exactly 0, Gaussian <z·z^2>=0): 9.2e-4 ± 6.5e-4 — consistent
  with 0 but with a variance comparable to the F-ON "signal".
* n_lambda convergence (same seed, F on): ratio = 39, 78, 49, 107 at n_lambda =
  300, 600, 1200, 2400. NON-CONVERGENT — the "signal" is noise.
* sigma_lambda scan (gamma=10'): FK = 4.6e-3, 3.8e-3, 2.0e-3, 7.6e-4 at
  sigma = 15, 30, 60, 120 — a small-sigma plateau at ~90x (so it is a normalisation
  / variance pathology, not a clean 1/(2 sigma) factor).

### Root cause (NOTE: "Q~3e10 / df~1000x z" below was an ARITHMETIC ERROR;
### corrected in the VR section -- real Q~-1.4e7 (sigma-indep), df/z~0.7-1.3,
### skewness~4-6. The naive estimator's failure is the SAME catastrophic-
### cancellation bug fixed in simulate_fk_vr, plus the genuine skewness-driven
### Riccati blow-up of the FULL nonlinear skewed arm.)
The equal-time zeta is Levy-peaked: V ~ 1e-14, zeta ~ 1e-17, so
`Q = solve_Q ~ zeta/(3 V^2) ~ 3e10` and `df = 0.5 Q(zz - V) ~ 1e-4` is ~1000x the
field `z ~ 1e-7`. The field 3-cumulant injection is CORRECT (int int <fff> = Zeta6,
verified analytically: AR(1) gives (2 sigma)^2 * node-mu3, and node-mu3 =
Zeta6/(2 sigma)^2). Path A is exactly LINEAR in df, so its MEAN is the correct
leading-Q FK and it does NOT blow up. But `Var(kappa_d) ~ Q^2 V^2` is enormous, so
`<kappa_g kappa_d>` is a tiny mean on a gigantic variance — it neither converges in
n_lambda nor beats its own error bar. The raw deformed field f_s = z + 0.5Q(zz-V)
also has a node <f_s^3> measured 3x its target (Q^2/Q^3 contamination), confirming
the deformation is non-perturbatively large (but that contamination is in the FULL
field, NOT in Path A's strictly-linear mean).

### Variance-reduced Path A — IMPLEMENTED + tested. Verdict: deformation-MC FK is
### INFEASIBLE for this Levy zeta (a real negative result).
`simulate_fk_vr` in `sachs_mc_core.py` pulls Q OUTSIDE the MC average. The FK
diagram connects obs_b to one zeta leg F-FREE, so the no-F windowed response
`kappa_d^(0)_B = int W delta_f_B` suffices (no s_d ODE). With
`delta_f_B = 0.5 sum_ij Q[k]_Bij w_ij(k)`, `w_ij = z_i z_j - V_ij`,

    <kappa_g_A kappa_d^(0)_B> = 0.5 dlam sum_k W(k) sum_ij Q[k]_Bij <kappa_g_A w_ij(k)>,
    FK6 = T + T^T  (x<->y symmetrisation; cross-ray block = FK).

Q multiplies the AVERAGED Cc, never the per-sample w. Magnitudes (CORRECTED):
Q = solve_Q = Zeta6/(3 Sigma2^2) ~ -1.36e7, SIGMA-INDEPENDENT (the earlier "Q~3e10"
was an arithmetic error). df = 0.5 Q(zz-V) has rms COMPARABLE to z (df/z ~ 0.7-1.3),
NOT 1000x; the field SKEWNESS ~ -4 to -6 (grows as sigma^-0.5, the Levy signature).

### The 35x systematic was a BUG (catastrophic cancellation) — FIXED.
The connected signal Cc = <kappa_g_A w_ij(k)>_c = Cov(kappa_g_A, z_i z_j) is ~1e7x
SMALLER than the disconnected <kappa_g_A><z_i z_j> = meankap * <z_i z_j>. The first
version subtracted the THEORETICAL V[k]; the AR(1) SAMPLE <z_i z_j> differs from
V[k] by transient/discretization terms, so meankap*(<zz>_sample - V[k]) ~ 1e7x the
signal leaked through -> a ~35x bias + huge scatter (the buggy grid sweep gave
ratios 126,33,88,12 at sigma=15,30,60,120 with std~mean -> pure noise). FIX
(`simulate_fk_vr`): subtract the MATCHED SAMPLE covariance
Cc_k = <kg_A z_i z_j>_batch - meankap * <z_i z_j>_batch (both from the same batch).
After the fix the magnitude is RIGHT: |FK_vr| ~ fk_analytic (e.g. -4.5e-5 vs
+5.4e-5 at gamma=10', sigma=30, nl=600; -6.3e-5 at sigma=60) -- a 35x -> ~1x
improvement. (User's "small correlation length -> denser grid" prompt is what
triggered finding this.)

### Precision fix: F-INDUCED CRN (the win that made it work).
The covariance noise was dominated by kappa_g's LINEAR part L1 (Cov(L1, zz) = 0 in
mean but adds variance ~ Var(L1) ~ Order-0 ~ (0.03)^2). Using the F-INDUCED part
`kappa_g(F-on) - kappa_g(F-off)` (shared-noise CRN) keeps the SAME signal (only L2+
couples to zz) while dropping the L1 variance ~1000x. `simulate_fk_vr` now evolves
both arms and uses the difference. Result at gamma=10' (4-seed):
    sigma=15 nl=1000: FK = 5.18e-5  (ratio 0.95, SE 4.8e-6)
    sigma=30 nl=1000: FK = 6.09e-5  (ratio 1.12, SE 4.0e-6)
    sigma=30 nl= 600: FK = 8.52e-5  (ratio 1.56, SE 1.4e-6)
    sigma=60 nl= 600: FK = 6.72e-5  (ratio 1.23, SE 1.2e-6)
SIGN is now POSITIVE (matches fk_analytic) -- the earlier negative was a noise
artifact of the L1-contaminated covariance, NOT a convention flip. SE ~ 1-5% (was
~100%).

GAMMA SWEEP (sigma=18, nl=1000, 4-seed; `fk_mc_vs_analytic.py`,
`figures/fk_mc_vs_analytic.pdf`): MC/analytic ratio is FLAT ~1.04-1.23 across
gamma = 1'..114' (1.14, 1.11, 1.10, 1.04, 1.10, 1.23), SE ~2-5%. The MC tracks the
analytic gamma-shape and magnitude to ~10%.

GRID CONVERGENCE -- HONEST, NON-MONOTONIC (do NOT over-claim "denser -> 1"). At
FIXED sigma, nl 600 -> 1200 nudges the ratio DOWN slightly (sigma=20: 1.22->1.20;
sigma=30: 1.36->1.26), but nl=2400 JUMPS UP at BOTH sigma (1.82, 1.90) -- a
consistent fine-grid UPWARD bias (SE ~5% but the jump is ~50%, so real, not one bad
seed). Likely heavy-tailed statistics (skewness ~4-6 -> rare large-positive
excursions sampled more at fine dlam) and/or a fine-grid discretization effect; the
5-seed SE underestimates the true (heavy-tailed) spread there. So the clean regime
is MODERATE grid (nl ~ 600-1200, sigma ~ 15-30) -> ratio ~ 1.0-1.35.

### Bottom line  (REVISED -- the MC FK WORKS to ~10-25%)
The variance-reduced MC FK (`simulate_fk_vr`: Q-outside + sample-covariance +
F-induced CRN) reproduces `fk_analytic` to ~10-25% with SIGN correct and a roughly
FLAT gamma-trend (~1.2 across 1'..114', at the figure's n_real=24000). This is a GENUINE independent MC confirmation
of the FK (3-cumulant) channel from the stochastic-ODE route -- the earlier
"infeasible/35x" conclusion was WRONG (it rested on the Q~3e10 arithmetic error +
the catastrophic-cancellation bug). BUT do not over-claim a clean ~1.000: there is a
persistent ~10-25% high offset (finite-sigma + a fine-grid heavy-tailed upward
bias), and naive grid-refinement does NOT monotonically -> 1. A precise (<10%)
sigma->0 limit needs (a) understanding/controlling the fine-grid bias, (b) far more
samples (heavy tails), or (c) a better estimator. `fk_analytic` (0.992 vs analysis-3)
remains THE precise FK number; the MC corroborates the channel's reality, magnitude,
and gamma-shape to ~10-25%.

## 3. FK MC convergence -- the knob is sigma_lambda, NOT n_lambda (2026-06-10)

Resolving the "denser grid -> ratio -> 1" question from section 2. The two grid
scales pull in OPPOSITE directions, so they must NOT be swept together:

* `n_lambda` (sets dlambda): at FIXED sigma_lambda a FINER grid makes the FK ratio
  WORSE. Decisive sweep at gamma=10' (6 seeds, n_real=9000): nl=2000 (dlambda~0.95)
  gives ratio 21.2, 20.0, 16.6, 14.3, 9.8 at sigma_lambda=18,14,10,8,6 -- a huge
  fine-grid UPWARD inflation (the heavy-tailed Levy zeta, skewness ~4-6, samples
  rare large excursions more as dlambda shrinks). This CONFIRMS section 2's nl=2400
  jump; n_lambda is NOT the convergence knob.

* `sigma_lambda` (AR(1) correlation length): at a FIXED MODERATE grid nl=1000
  (dlambda=1.91, sigma_lambda/dlambda ~ 4-9, AR(1) resolved) decreasing
  sigma_lambda drives the FK ratio MONOTONICALLY through 1. gamma=10' (6 seeds):

      sigma_lambda  18    14    11     9     7     5     4
      sigma/dlambda 9.4   7.3   5.8   4.7   3.7   2.6   2.1
      FK MC/analytic 1.99  1.69  1.35  1.08  0.79  0.52  0.42

  i.e. a clean monotone descent crossing 1 near sigma_lambda~9 (sigma/dlambda~4.7).
  Below sigma_lambda~7 it undershoots (genuine sigma->0 over-shoot + AR(1) under-
  resolution at sigma/dlambda<4). sigma_lambda IS the convergence knob: the analytic
  FK is the sigma_lambda->0 (delta / white-noise) limit, and the MC approaches it as
  the colored-noise width shrinks.

ADOPTED FIGURE SETTING: sigma_lambda=8, nl=1000 (sigma/dlambda=4.2, still AR(1)-
resolved), 8 seeds, n_real=24000/seed -- the `fig_xi_channels.py` default that builds
the published cache `outputs/appendix_mc_curve.npz`. FK MC/analytic across the full
gamma grid:

    gamma'        1.0   2.6   6.7  17.3  44.4 114.3 293.9   median
    FK MC/analytic 1.26 1.18 1.23 1.08 1.12 1.31 1.36    1.23

i.e. FK MC/analytic ~ 1.2-1.25 and roughly flat in gamma (range [1.08, 1.36] over
gamma<300'), down from ~2.4 at the old sigma_lambda=18 -- a ~20-25% finite-
sigma_lambda overshoot (matches the paper caption's "~25%"). The ratio drifts UP
with sample count (a lower-stat 9000/seed `_fk_converge_validate.py` run gave
median ~1.12); the near-delta Levy zeta has skewness ~4-6, so the estimator is
heavy-tailed with ~+/-40% inter-seed spread -- an order-of-magnitude corroboration,
not a converged precision number (`fk_analytic` = 0.992 vs analysis-3 remains THE
precise FK value). The monotone sigma_lambda-trend above is the validation: the
analytic FK workflow IS the converged (sigma_lambda->0) limit of the direct
stochastic-Sachs MC. This is the figure `figures/xi_kappa_channels.pdf` (regenerated
at this setting; the old sigma_lambda=18 version archived as
`xi_kappa_channels_archive_2026-06-10_*.pdf`).
