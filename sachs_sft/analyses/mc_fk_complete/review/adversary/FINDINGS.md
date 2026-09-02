# Adversarial review: "MC and exact agree only because they share the author's bugs"

Scripts a1..a8 in this directory; each prints its own falsification controls.
All runs: PyCCL env, OMP_NUM_THREADS=2, <= 2 processes, n_real <= 24000, batch 2000.

## Verdict
The NARROW form of the attack is CORRECT and is demonstrated in a3: MC/exact is
structurally blind to every convention carried by the Grid.  The BROAD form fails:
every named convention was tested against an object outside the shared code and
survived, and the load-bearing comparison (against the paper's fold) has teeth.

## a3 (central): shared-bug injection, gamma = 0.63319' node, sigma = 8, N = 1000
world                        exact/fold   ref(s->0)/fold   MC/exact
response ^2 (as coded)         0.9719        0.9992        0.9988 +- 0.0297
response ^2 -> ^2.06           0.9669        0.9943        0.9988 +- 0.0294
response ^2 -> ^2.4            0.9476        0.9756        0.9988 +- 0.0275
response ^2 -> ^4              1.0049        1.0376        1.0002 +- 0.0214
response ^2 -> ^1              1.1646        1.1915        1.0012 +- 0.0371
NOTE: at sigma = 8 the fold comparison alone PREFERS the wrong ^4 world.
Discrimination lives in the sigma -> 0 extrapolation, not at finite sigma.

## a2/a2b/a2c: Order-0 anchor (no F, no Q, no zeta, no Wick)
Richardson-extrapolated discrete white-noise Order-0 / converged quadrature of the
same formula = 0.99999 at 1.015', 5.30', 17.28'.
driver_stats.order0_mc's 48-node Gauss-Legendre is 0.9894/0.9913/0.9926 of its own
converged value -- the ~1% "residual" is that quadrature, not the discrete kernel.
Teeth (white limit / order0_mc): exponent 2->1 1.111-1.248, 2->4 1.094-0.937,
window reversed 847-1626, no window 1e-4, drop trapezoid ends 1.046-1.064.
V = Sigma2 (no /2sigma) -> 9.06-16.2 ; V = Sigma2/sigma -> 2.03-2.26.

## a8: the (D_la/D_t)^2 exponent, against sft-wick analysis-3 over 0.5'-145'
ratio spread (max/min - 1) of O0_eqtime/analysis3:  p=1 1.29 | p=2 0.029 | p=3 1.82 | p=4 115.5

## a1: F vertex, its linearisation, AR(1), propagator
_f_vertex6 vs F6 3e-18 ; _f_vertex6_lin vs FS6 3e-17 ; vs central finite difference 2e-13.
Controls: no symmetrisation 0.84, spurious factor 2 -> 1.0.
AR(1) R = rho^|k-p| A[min] vs direct simulation: 0.001-0.004 (noise 0.0013);
controls A->V 0.10-1.11, rho^{2 lag} 0.16-0.24.
Impulse response = dlam (D_p/D_k)^2 to 1e-15; exponent-4 control 0.58-0.91.

## a4: brute-force injection (sampling only, no Wick)
zeta_inj vs sym(zeta_tab): 8.4e-12 -- the sigma->0 answer is fixed BY CONSTRUCTION.
Central-difference O(h) third moment of the actual field / sym(zeta_tab), 4e6 samples:
0.9945-1.0035 at three nodes.  Controls: 0.5 -> 1.0 reads 1.995-2.006;
Q against Sigma2 (nominal) reads 0.640-0.846.
(T1+T2)[0,3] == (T1+T2)[3,0] exactly, so for the kk entry "+transpose" IS a bare
factor 2: it is justified by the derivation and by the fold, never by a symmetry.

## a5/a5b: chain <kappa1^3> (no F vertex), 12 x 40000
smeared, C matched to the simulated chain: 1.048/1.019/0.976/0.997 (+-0.018..0.025).
Controls: nominal calibration 3.0-5.5 ; C = V instead of the chain variance -0.52..-4.87.
_cholesky_psd floors 2 nodes of 1000; the exact expectation run on the FLOORED V
differs by 1.000136, i.e. the MC/exact model mismatch is 0.014% of the FK
(NOTES' "~6% of the kernel weight" overstates it).

## a6/a7: gate 4 at EXACT fold nodes (no gamma interpolation), N = 1000
sigma->0 extrapolated exact / fold: median 1.00069, range [0.99999, 1.00208] over 0.5'-17.3'
local reference (injected zeta) / fold: median 0.99979, range [0.99912, 1.00133]
My own MC, 16 seeds x 24000, gamma = 1.01546': MC/exact 1.0014 +- 0.0057 (s=8),
0.9928 +- 0.0101 (s=4); exact/fold 0.9729 / 0.9866; exact extrapolated 1.0002.

## Residuals worth reporting (none fatal)
1. Kernel not converged at N=1000: reference/fold runs 1.00146 (N=500) -> 0.99766
   (N=8000), still drifting.  The converged kernel sits ~0.23% BELOW the fold while
   N=1000 reads 0.99936.
2. run_gate2.analytic_fk interpolates the fold linearly in gamma; at 1.0' that is
   +0.09% vs log-log interpolation.
3. Raw vs symmetrised zeta legs is worth 0.9% at 1.015' (1.600698e-05 vs 1.586411e-05)
   and more at large gamma; gate 4 is flat only against the symmetrised object.
4. Discretisation choices still worth ~1% at N=1000: F at old vs new state 0.13%,
   trapezoid end weights 1.28%; both scale as dlam and vanish by N=8000.
5. The regulator is not weak at production settings: the discrete Order-0 two-point at
   N=1000, sigma=8 is 0.969 of the white-noise value.  FK is invariant under
   Sigma2 -> c Sigma2 so this does not bias FK, but the extrapolation does real work.
