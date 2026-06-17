# Jacobi flat-FLRW solution: D(lambda) = a(lambda) chi(lambda)

## Statement

On a flat FLRW background with matter, radiation, and a cosmological
constant, the Jacobi geodesic-deviation equation

    D''(lambda) = bar Phi_00(lambda) D(lambda),
    D(0) = 0,  D'(0) = 1

with the past-directed null-geodesic affine parameter (E_0 = 1) admits the
closed-form solution

    D(lambda) = a(lambda) chi(lambda),

where a is the scale factor and chi is the comoving radial coordinate
along the line of sight. Equivalently,

    theta^(sa)(lambda) := D'(lambda) / D(lambda) = (d/d lambda) ln D(lambda),

so

    integral_{lambda'}^{lambda} theta^(sa)(tau) d tau = ln D(lambda) - ln D(lambda'),
    R(lambda, lambda') = exp(-2 integral_{lambda'}^{lambda} theta^(sa) d tau)
                       = (D(lambda') / D(lambda))**2  for  lambda >= lambda'.

The full proof is symbolic, written in `scripts/derive_jacobi_flat_flrw.wl`
and runs without xAct; it uses only Mathematica's plain symbolic calculus.

## First principles

1. **FLRW metric in conformal time eta.**
   `ds^2 = a^2(eta) [-d eta^2 + dx_i dx^i]`. Source:
   `sections/cosmology.tex` Eq. `poisson gauge` with the Poisson-gauge
   perturbations `Psi, Phi, B, h_ij` set to zero.

2. **Past-directed null geodesic with E_0 = 1.** The radial null tangent
   has `k^eta = -1/a^2, k^chi = +1/a^2, k^perp = 0`, so
   `d eta / d lambda = -1/a^2` and `d chi / d lambda = +1/a^2`. Source:
   `sections/cosmology.tex` Eq. `dchi dlambda bg`; symbolic derivation in
   `scripts/ray_coord_map.wl` Phase 6 (`d chi/d lambda = EE/a = 1/a^2`).

3. **Sachs scalar.** `Phi_00 := -(1/2) R_{mu nu} k^mu k^nu` with the paper
   sign convention. Source: `sections/sachs_dynamics.tex` Eq.
   `Phi 00 from Ricci`. On flat FLRW this reduces to
   `bar Phi_00(lambda) = -(calH^2(eta) - calH'(eta)) / a^4(eta)`, with
   `calH := a' / a` (prime = d/d eta). Source: `sections/cosmology.tex`
   Eq. `driving bg`. Numerical implementation in
   `scripts/ccl/driving_bg.calH2_minus_calHprime_lcdm`.

4. **Jacobi (geodesic deviation) equation.** Source:
   `sections/sachs_dynamics.tex` Eq. `Jacobi D` and surrounding derivation.

5. **Friedmann pair.** On flat FLRW with matter, radiation, and Lambda,
   `3 calH^2 = 8 pi G a^2 (rho_m + rho_r + rho_L)` and
   `calH^2 - calH' = 4 pi G a^2 (rho_m + (4/3) rho_r)`. Equivalently in
   fractional form `calH^2 - calH' = calH^2 [(3/2) Omega_m(a) + 2 Omega_r(a)]`,
   the identity used by `scripts/ccl/driving_bg.calH2_minus_calHprime_lcdm`.

## Sign convention

`Phi_00 < 0` in matter and radiation eras, by the paper convention. The
Jacobi equation `D'' = bar Phi_00 D` is correct as written; no extra
minus sign. With `bar Phi_00 < 0` and `D > 0` we have `D'' < 0`, recovering
the focusing of a thin null beam in a matter-dominated universe.

## Symbolic-proof outline

The Mathematica script verifies five identities labelled T1 through T5.

**T1.** `D = a chi` satisfies `D'' = bar Phi_00 D`. Use the chain rule:
`d eta/d lambda = -1/a^2` and `d chi/d lambda = 1/a^2`, then
`d a/d lambda = -calH/a` and `d^2 a/d lambda^2 = (calH^2 - calH')/a^3`
(after substituting `a''(eta) = a (calH^2 + calH')`). The product rule
gives `D'' = -2 (calH^2 - calH') chi / a^3 + 0 = -(calH^2 - calH') (a chi) / a^4
       = bar Phi_00 D`. Residual after `Simplify` is exactly zero.

**T2.** Initial conditions. At `lambda = 0`: `a = 1, chi = 0`, so `D(0) = 0`.
And `D'(lambda) = -calH chi / a^2 + a / a^2`; at `lambda = 0` this collapses
to `0 + 1 = 1`.

**T3.** `calH^2 - calH' = calH^2 [(3/2) Omega_m + 2 Omega_r]`. Mathematica
solves the Friedmann pair for `(calH^2, calH')` symbolically in terms of
`(rho_m, rho_r, rho_L, G, a)`, computes the LHS, then uses
`rho_crit = 3 calH^2 / (8 pi G a^2)` and
`Omega_x = rho_x / rho_crit` to assemble the RHS. The residual collapses
to zero.

**T4.** `theta^(sa)(tau) - (d/d tau) ln D(tau) = 0`. Direct simplification.
Then by the fundamental theorem of calculus,
`int_{lambda'}^{lambda} theta^(sa) d tau = ln D(lambda) - ln D(lambda')`.

**T5.** `exp(-2 (ln D(lambda) - ln D(lambda'))) - (D(lambda')/D(lambda))**2 = 0`.
Direct algebraic simplification.

## Verbatim wolframscript output

```
=======================================================
  Jacobi flat-FLRW: verify D(lambda) = a(lambda) chi(lambda)
=======================================================
Phase 1:  da/dlambda      = -(Derivative[1][a][eta]/a[eta]^2)
          d2a/dlambda^2   = (-2*Derivative[1][a][eta]^2 + a[eta]*Derivative[2][a][eta])/a[eta]^5
          d chi/dlambda   = a[eta]^(-2)
          d2chi/dlambda^2 = (2*Derivative[1][a][eta])/a[eta]^5

Phase 2:  D'(lambda)  = (a[eta] - chi*Derivative[1][a][eta])/a[eta]^2
          D''(lambda) = (chi*(-2*Derivative[1][a][eta]^2 + a[eta]*Derivative[2][a][eta]))/a[eta]^5

Phase 3:  bar Phi_00 D = (chi*(-2*Derivative[1][a][eta]^2 + a[eta]*Derivative[2][a][eta]))/a[eta]^5
          T1 residual D'' - bar Phi_00 D = 0
          (PASS iff residual is 0)

Phase 4:  D(0)  = a(0) chi(0) = 0   (expected 0)
          D'(0) = 1   (expected 1)

Phase 5:  Friedmann pair => fractional form  ...
          calH^2 - calH'                              = (4*G*Pi*(3*rhoM + 4*rhoR)*a[eta]^2)/3
          calH^2 [(3/2) Omega_m + 2 Omega_r]          = (4*G*Pi*(3*rhoM + 4*rhoR)*a[eta]^2)/3
          T3 residual (LHS - RHS)                     = 0
          (PASS iff residual is 0)

Phase 6:  theta^(sa) integral and R-propagator  ...
  T4: theta^(sa)(tau) - d/dtau ln D(tau) = 0    (PASS iff 0)
      => by FTC,  int_{lp}^{l} theta^(sa) d tau
                   = ln D(l) - ln D(lp).
  T5: exp(-2 (ln D(l) - ln D(lp))) - (D(lp)/D(l))^2 = 0    (PASS iff 0)

=======================================================
Summary
  T1: D'' = bar Phi_00 D            -> PASS
  T2: IC D(0)=0, D'(0)=1            -> PASS
  T3: calH^2 - calH' fractional ID  -> PASS
  T4: int theta^(sa) = ln D(l)-ln D(lp) -> PASS
  T5: R = (D(lp)/D(l))^2            -> PASS

[PASS]  All five identities verified.
=======================================================
```

## Numerical cross-check note

The closed-form solution is the FLRW reference that
`scripts/ccl/sft_wick_analysis/saddle_point.theta_saddle_analytical`
computes. The numerical ODE solve in
`scripts/ccl/sft_wick_analysis/saddle_point.solve_jacobi_D` integrates
`D'' = bar Phi_00 D` with `solve_ivp(method="DOP853", rtol=1e-10)` using
the analytical `bar Phi_00` from `driving_bg.calH2_minus_calHprime_lcdm`,
and reaches relative agreement of better than 1e-7 against the analytical
`a chi` away from the integration's IC point. See
`scripts/ccl/sft_wick_analysis/tests/test_lightcone_and_jacobi_precision.py`
for the regression assertions.

## Downstream consequences

Once D is available numerically (cached in
`scripts/ccl/sft_wick_analysis/inputs/theta_sa.npz`), every theta^(sa)
integral in the analysis pipeline can be evaluated in closed form via T4
and T5. The canonical API is in
`scripts/ccl/sft_wick_analysis/D_callable.py`:

  - `D_at(lam)` returns the Jacobi amplitude.
  - `theta_sa_integral(lam_upper, lam_lower)` returns the integral via
    `ln D(lam_upper) - ln D(lam_lower)`.
  - `response_R(t, lam_prime)` returns
    `Theta(t - lam_prime) * (D(lam_prime) / D(t))**2`.

This eliminates the 1/lambda integrand singularity that made the
`cumulative_trapezoid(2 * theta^(sa))` path overshoot near the vertex
on coarse grids.
