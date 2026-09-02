# Why the two-ray Sigma2 is not positive semi-definite

*2026-08-28. Diagnosis only; no file outside `psd_fix/` was modified.*

## Root cause, in one paragraph

`corr_op` serves its 21x21 lambda table through a **bilinear**
`RegularGridInterpolator` (`interp_method` defaults to `"linear"` and the
callable does not override it), so on the equal-time diagonal `C_ab(t,t)` is
*exactly piecewise quadratic* and only **`C^0`**: the first derivative jumps, and
at `lambda = 1300` it flips sign, `+5.987e-13` to `-2.209e-13`. The true density
is therefore **discontinuous** at every interior node.
`Sigma2Builder._density_splines` fits a globally `C^2` interpolating cubic spline
through 160 samples of that `C^0` function and takes its **analytic derivative**.
A `C^2` spline cannot represent a jump, so it rings in the cells around each
node; the ringing overshoots the smallest eigenvalue and drives it negative. At
`lambda = 1308.93`, `gamma = 1'`, all six eigenvalues go negative, so
`_cholesky_psd` returns `L = 0` and the simulation injects **no driving noise at
all** at that node. This is a smoothness mismatch, not noise, and not a property
of the ray geometry.

## The chain, with the numbers that pin each link

**1. The defect is in the within-ray block, not the assembly.** At
`lambda = 1308.9` the smallest eigenvalue of the within-ray `A` is `-1.879e-14`,
identically at `gamma = 1'` and `17.3'` (`A` does not depend on gamma). The 6x6
inherits it; `A-C` and `A+C` are both negative there. So `_assemble_6x6` is
exonerated.

**2. It is a spline artefact: it moves with the sampling.**

| `n_nodes` | spacing | non-PSD lambda (of 1000) | worst eigenvalue / median |
|---|---|---|---|
| 80 | 24.3 Mpc | 0 | +0.018 |
| 120 | 16.2 | 0 | +0.018 |
| 160 (production) | 12.1 | 2 | **-0.158** |
| 240 | 8.0 | 3 | -0.677 |
| 320 | 6.0 | 5 | -1.026 |
| 480 | 4.0 | 2 | -0.119 |
| 640 | 3.0 | 1 | -0.278 |

Non-monotone, and absent at coarse sampling. Refining the grid makes it worse.

**3. The true density is PSD; the spline derivative is simply wrong.** Central
finite differences of the same corr_op data, at `h` = 0.5, 4 and 12 Mpc, agree
with each other to three digits and are positive everywhere:

| lambda | spline derivative, min eig | FD 0.5 Mpc | FD 4 Mpc | FD 12 Mpc |
|---|---|---|---|---|
| 1307.0 | +1.216e-14 | +5.399e-14 | +5.401e-14 | +5.418e-14 |
| **1308.9** | **-1.879e-14** | +6.053e-14 | +6.055e-14 | +6.072e-14 |
| 1312.7 | +2.752e-14 | +7.372e-14 | +7.374e-14 | +7.391e-14 |

The increments of `D^4 C` are PSD at every step size tested (0/280 negative from
0.5 to 12 Mpc), confirming the accumulated object is a genuine covariance.

**4. The source of the kink: the corr_op table's own 21-point lambda grid.**

```
[406, 550, 700, 850, 1000, 1150, 1300, 1447, 1575, 1700, 1900, 2037,
 2100, 2160, 2250, 2280, 2300, 2309, 2314, 2318, 2328]
```

`C_00(t,t)` fits a quadratic to a residual of `7e-16` on each side of a node,
confirming the bilinear interpolant. The **first** derivative jumps sign across
`lambda = 1300` (`+5.987e-13` to `-2.209e-13`) and across `lambda = 850`
(`+1.892e-13` to `-9.031e-14`). So `C` is `C^0` and the density `Sigma2` is
**discontinuous** at each table node, falling sharply across it.

**5. Every pathological window sits just past a table node.** Eleven windows,
matched one-to-one with the interior nodes:

| window [Mpc] | nearest table node | offset |
|---|---|---|
| 554.9 - 562.5 | 550 | +4.9 |
| 856.5 - 866.0 | 850 | +6.5 |
| 1003.5 - 1011.1 | 1000 | +3.5 |
| 1160.0 - 1165.7 | 1150 | +10.0 |
| 1303.2 - 1314.6 | 1300 | +3.2 |
| 1450.2 - 1459.7 | 1447 | +3.2 |
| 1581.9 - 1589.5 | 1575 | +6.9 |
| 1704.1 - 1713.6 | 1700 | +4.1 |
| 1910.2 - 1912.1 | 1900 | +10.2 |
| 2041.9 - 2047.7 | 2037 | +4.9 |
| 2104.9 - 2108.7 | 2100 | +4.9 |

**6. The ringing itself.** Relative error of the production spline derivative
against a 4 Mpc finite difference, cell by cell across `lambda = 1300` (the
160-node spline cells there are at 1288.4, 1300.5, 1312.6, 1324.7):

| lambda | 1292 | 1298 | 1301 | 1305 | 1307 | **1309** | 1313 | 1317 | 1325 |
|---|---|---|---|---|---|---|---|---|---|
| error | +13% | -22% | +47% | +7% | -77% | **-132%** | -56% | -2% | +10% |

The alternating-sign overshoot of a cubic spline forced through a kink.

## Candidate repairs, measured

Implemented as drop-in `Sigma2Builder` subclasses in `sigma2_repaired.py`;
compared over 600 lambda at `gamma = 1'` and `17.3'` against `R0Reference`
(dense per-cell finite differences, no interpolation).

| | non-PSD | median error vs R0 | Order-0 anchor (panelised) |
|---|---|---|---|
| **current** | 2/600 | 6.1e-3 | 1.00000 |
| **R1** knots at the table nodes | **0/600** | **1.0e-5** | 0.99970 |
| R2 density + PCHIP | 0/600 | 1.1e-4 | -- |
| R3 smoothing spline | 2/600 | 6.7e-3 | -- |
| R0 reference | 0/600 | 0 | -- |

*Max-error figures are omitted deliberately: `R0Reference` is not a valid
yardstick within ~1 Mpc of a table node, where the true density is
discontinuous, so any max taken over a grid that includes those points measures
the yardstick, not the candidate. The medians are robust and are what decide it.*

**R3 failing is diagnostically important.** A penalised smoothing spline is the
remedy for noise, and it does not help at all. That independently confirms the
defect is a smoothness mismatch.

**R1 is the recommended repair**: the same interpolate-then-differentiate recipe
as production, but with the spline broken at the 21 corr_op table nodes so no
cell straddles a kink. It removes every PSD violation, is 200x more accurate than
production in the worst case, and tracks the exact reference to 0.006% on the
Order-0 anchor.

## Effect on the validated Order-0 anchor: -0.03%, not +1.61%

An early version of this note reported `+1.61%`. That was a **quadrature
artefact**, not a property of the repair. `driver_stats.order0_mc` applies one
global Gauss-Legendre rule to the integrand; the repair makes that integrand
discontinuous at 21 points, and a global rule cannot integrate a discontinuous
function. Panelising the same rule at the corr_op table nodes, so no panel
straddles a jump:

| gamma | 0.5' | 1' | 5' | 17.3' | 114.3' | 293.9' |
|---|---|---|---|---|---|---|
| R1/stock, global GL n=256 | 1.01558 | 1.01564 | 1.01663 | 1.01766 | 1.01991 | 1.02994 |
| **R1/stock, panelised GL** | **0.99970** | **0.99970** | **0.99973** | **0.99978** | **0.99995** | **0.99983** |

The global rule's answer is both wrong and gamma-dependent; the panelised one is
`-0.03%`, flat in gamma. Consequences:

* `ANCHOR_RATIO = 0.881` does **not** need re-deriving.
* `order0_mc` **does** need its quadrature panelised at the table nodes, in the
  same commit as the repair. Without that it will read 1.6-3.0% high, by a
  gamma-dependent amount.

## What the repair does and does not justify

**"The negative eigenvalue disappeared" is not a justification.** Five of seven
deliberately wrong break placements (nodes shifted +/-30 Mpc, cell midpoints, 19
uniform breaks, +/-2 Mpc jitter) also give 0 non-PSD: breaking a global spline
anywhere shortens its support and damps the ringing.

The real justification is that the density is **analytic** given the bilinear
interpolant: `C(t,t)` is piecewise quadratic, so the density is piecewise linear,
with zero free parameters. R1 reproduces it to a median `8.4e-6`, and the true
table nodes are a sharp optimum: the median error degrades 10x at `+/-2 Mpc`,
14x at `+/-8 Mpc`, and 280-730x at `+/-30 Mpc` or at cell midpoints, where it is
as bad as the unrepaired builder. That is a structural signature, not a fitted
one.

## The larger defect this uncovered, which the repair does NOT address

The same bilinear interpolation is applied to the full 21x21 table, and `C` has a
**ridge along the equal-time diagonal**. At a cell midpoint the bilinear value
averages the two on-diagonal corners with two off-diagonal ones and undershoots.
Measured against a 1-D interpolation through the 21 exactly-stored diagonal
values:

| cell | 406-550 | 700-850 | 1000-1150 | 1300-1447 | 1575-1700 | 1900-2037 | 2160-2250 |
|---|---|---|---|---|---|---|---|
| bilinear / nodal at midpoint | **0.730** | 0.789 | 0.824 | 0.846 | 0.872 | 0.884 | 0.977 |

The mechanism, at the midpoint of `[1300, 1447]`: on-diagonal corners
`9.281e-11` and `1.449e-10`, off-diagonal corners `7.657e-11` twice; the bilinear
average is `9.772e-11` against `1.156e-10` for the diagonal's own interpolation.

So the equal-time `Sigma2` is systematically **low by 4-27% between table nodes**,
far larger than anything the ringing repair fixes. Refining the table is the
wrong lever (the error falls only 1.4-1.7x per halving of the cell). The right
one is a ridge-aware, one-dimensional reading of the diagonal, which the table
already stores exactly at all 21 nodes.

**Scope of that larger defect.** The FK diagram carries four `R` propagators and
**zero `C` propagators**, so it never touches this table and is immune. Order-0
and FF do use `C`, at general `(lambda_1, lambda_2)` rather than only on the
diagonal, so they could carry a percent-level bias. **That has not been
measured** and is the open question worth pursuing.

## Reproducing

```
export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
cd driver_field_emulators/code/mc_fk_complete
$PY psd_fix/d1_localise.py    # within-block, n_nodes dependence
$PY psd_fix/d2_noise.py       # true density is PSD; FD agrees across h
$PY psd_fix/d3_grid.py        # the corr_op 21-point lambda grid
$PY psd_fix/d4_kink.py        # curvature jump, window-to-node matching, ringing
$PY psd_fix/d5_compare.py     # candidate repairs: PSD and accuracy
$PY psd_fix/d6_anchor.py      # Order-0 anchor before and after
```
