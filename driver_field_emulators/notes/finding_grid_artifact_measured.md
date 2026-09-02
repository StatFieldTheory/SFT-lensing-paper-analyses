# Finding: the deployed cosine grid freezes the FK vertex at the (1,1,1) corner

Measured 2026-08-25 by building a kappa3 table densified on the collapsed
family and comparing it, at the vertex level, against the deployed
covgrid16 table as the callable actually reads it.

## Method

Both tables carry the SAME physics: tree-level bispectrum, ell_cut = 60,
ell_high_max = 1000, identical 16 lambda shells, identical cosmology. Only
the cosine grid differs:

* deployed: 15 uniform cosines (spacing 1/7), read through the callable's
  canonicalisation and 4-neighbour inverse-distance interpolation
* dense: the collapsed family (1, c, c) containing the 40 production cosine
  values exactly (nearest-row distance 0.0), read directly

So any difference is interpolation, not physics.

## Result: zeta_TTT on the collapsed family, mid shell

| gamma [arcmin] | deployed (IDW) | dense (direct) | ratio |
|---|---|---|---|
| 0.50 | -1.3148e-15 | -1.3122e-15 | 0.998 |
| 1.63 | -1.3148e-15 | -1.2933e-15 | 0.984 |
| 5.30 | -1.3148e-15 | -1.1021e-15 | 0.838 |
| 17.28 | -1.3145e-15 | -2.0736e-16 | 0.158 |
| 56.27 | -1.3119e-15 | -1.3920e-17 | 0.011 |
| 183.26 | -1.2840e-15 | +1.7772e-18 | -0.001 |
| 596.89 | -1.0319e-15 | +4.9479e-19 | -0.000 |
| 1944.08 | +1.1726e-18 | -3.5185e-18 | -3.001 |
| 5000.00 | +9.9370e-19 | -3.2350e-18 | -3.255 |

The deployed column is essentially CONSTANT from 0.5' to about 600': it is
the (1, 1, 1) corner value being returned again and again, exactly as the
weight analysis predicted. The true collapsed vertex meanwhile decays by
about three orders of magnitude and changes sign near 183'.

The two agree to 0.2% at 0.5', where the corner value IS the right answer.
That agreement is the control: it shows the two tables are mutually
consistent where the grid is adequate, so the divergence at larger angles
is the interpolation and not a build difference.

## The observable inherits it: measured

Both tables were folded through the production FK pipeline (the same
`run_FF_single.py`, same geometry, same n_gauss = 24). xi_FK for the kk
pair at order 2:

| gamma ['] | deployed | dense / deployed |
|---|---|---|
| 0.50 | +4.287723e-06 | 0.9985 |
| 5.30 | +4.287638e-06 | 0.8692 |
| 10.77 | +4.287371e-06 | 0.5706 |
| 17.28 | +4.286816e-06 | 0.2855 |
| 44.43 | +4.281728e-06 | 0.0735 |
| 90.24 | +4.263069e-06 | 0.0178 |
| 144.71 | +4.224707e-06 | 0.0027 |
| 596.89 | +3.365458e-06 | -0.0006 |

The published baseline at 0.5 arcmin is safe AGAINST THE COSINE-GRID
EFFECT: 4.2815e-06 against 4.2877e-06, a 0.15% difference. (A distinct
systematic, the linear-in-lambda interpolation of the sparse shell grid,
inflates both numbers by +11% at this gamma; the grid-converged tree
baseline is +3.84e-6. See `theory_redshift_dependence.md` section 2.2.
The two effects live in different table dimensions: cosine rows versus
lambda columns.) Everything beyond a few arcminutes is not safe even
against the cosine grid alone: the deployed FK curve barely decays out
to 20 degrees, while the corrected one falls by three orders of
magnitude by 145 arcmin.

## Consequence for the FK-vs-Order-0 crossover

Comparing against the Order-0 run on the same gamma grid:

| gamma ['] | O0 | \|FK/O0\| deployed | \|FK/O0\| dense |
|---|---|---|---|
| 0.50 | +8.4432e-04 | 0.0051 | 0.0051 |
| 8.51 | +3.1318e-04 | 0.0137 | 0.0096 |
| 56.27 | +2.8539e-05 | 0.1499 | 0.0078 |
| 144.71 | +2.1948e-06 | 1.9249 | 0.0053 |
| 372.19 | -3.7454e-07 | 10.3949 | 0.0088 |
| 957.24 | -8.6422e-08 | 27.3560 | 0.0038 |

The deployed FK overtakes Order-0 at 145 arcmin, about 2.4 degrees, which
is the crossover the paper reports. With the grid corrected the ratio is
essentially FLAT at 0.4 to 1% of Order-0 across the entire range: the two
terms decay together, and no crossover occurs anywhere in the physical
range (the formal crossing at 1535 arcmin is between two curves already at
the numerical floor and changing sign, so it is not meaningful).

The apparent wide-angle FK enhancement was therefore the frozen corner
value being compared against an Order-0 term that genuinely decays. This
is a like-for-like comparison: both tables are tree level at
ell_high_max = 1000 on identical shells, so only the cosine grid differs.

## The evidence internal to the deployed table

The comparison above uses a second table, so it invites the question of
whether the dense build, rather than the interpolation, is what differs.
Two checks answer it without leaving the deployed file.

**The deployed grid has no sampling inside the plotted range.** Its cosine
spacing is 1/7, so the first two grid points on the collapsed family are
c = 1 (gamma = 0) and c = 0.857 (gamma = 31 degrees). The paper's FK
figure spans 0.5 to 2000 arcmin, that is 0.008 to 33 degrees: the whole
curve lies inside the FIRST grid gap. Interpolating the deployed table's
OWN (1,c,c) rows in one dimension, bypassing the k-nearest-neighbour
lookup entirely, still gives -1.3148e-15 at 0.5', -1.3067e-15 at 145' and
-7.48e-16 at 1212'. The flatness is therefore not an artefact of the
inverse-distance scheme; the grid simply carries no angular information
where the figure lives.

**The deployed table's own next sample is three orders down.** On the same
shell, its rows read -1.3148e-15 at c = 1 and +8.64e-19 at c = 0.857. The
deployed file itself therefore says the vertex falls by more than three
orders of magnitude, and changes sign, somewhere inside the gap. What it
cannot say is where, which is what the dense grid measures.

Note also that neither table's measure is in question: both are produced by
the same canoes Limber branch with the same 1/chi^4, the same (1+z)^-4
collapse Jacobian and the same h^4. The paper's wide-angle FK curve is not
an independent analytic evaluation either; it is this table folded through
the same pipeline. The only difference is angular sampling.

## Why the endpoints being right does not save the middle

A fair objection: the grid is coarse, but its two endpoints, c = 1 and
c = 0.857, are exact table entries. Interpolation between two correct
values should not be badly wrong.

It is not wrong in the interpolation variable. It is wrong because the
interpolation variable is cos(gamma), and the queries are crowded against
one endpoint. Position of each production separation inside the first grid
interval [0.857, 1], measured as (1 - cos gamma) / (1/7):

| gamma ['] | 1 - cos gamma | position in the interval |
|---|---|---|
| 0.50 | 1.06e-08 | 0.0000% |
| 21.88 | 2.03e-05 | 0.0142% |
| 56.27 | 1.34e-04 | 0.0938% |
| 144.71 | 8.86e-04 | 0.62% |
| 372.19 | 5.86e-03 | 4.10% |
| 1212.23 | 6.15e-02 | 43.07% |

Twenty-six of the forty production separations sit within 1% of the c = 1
endpoint, so linear interpolation returns at least 99% of the endpoint
value there.

Two facts then combine. First, the c = 1 endpoint is the vertex at ZERO
separation, which is the maximum of a steeply falling function. Second,
small-angle geometry compresses: 1 - cos gamma is approximately
gamma^2 / 2, so three decades in angle, 0.5 to 145 arcmin, cover only
0.62% of the first cosine interval.

The result is that the deployed curve equals the zero-separation vertex
value almost everywhere in the plotted range, and the direction of the
error is fixed: it is an overestimate, because it returns the maximum. At
0.5 arcmin the query really does sit on the endpoint, which is why the two
tables agree there to 0.15%.

So both endpoints are indeed correct. The failure is that the plotted
range does not span the interval between them; it hides inside one end of
it.

## Caveats

* The dense table carries the collapsed family only. That is valid for the
  FK two-point fold, which queries nothing else (verified by an
  instrumented probe), but such a table must not be used for the order-1
  three-point vertex, which needs open triangles.
* Beyond about 1500 arcmin both curves sit at the numerical floor and
  change sign, so ratios there carry no information.
* The dense values at wide angles (beyond roughly 10 degrees) are not
  claimed to be accurate in detail. Cross-checking the two builds at
  comparable geometries there (deployed 31.0 degrees against dense 32.4
  degrees) gives ratios between -70 and 0.3, because the vertex is
  oscillatory and cancellation-dominated in that regime; the deployed
  table's own permutation duplicates of one canonical triple differ by up
  to a factor 15 there for the same reason. The conclusions rest on the
  range 0.5 to about 600 arcmin, where the dense grid is log-refined, the
  vertex decays smoothly, and the two builds agree to 0.2% at the corner.
* This says nothing yet about the nonlinear bispectrum. A nonlinear B
  raises the FK amplitude, most strongly at small angles and near shells,
  but the angular SHAPE measured here comes from the vertex's angular
  dependence, which the grid was hiding. The nonlinear run is the next
  measurement.

## Next

1. Rebuild the production 1671-triple table with the collapsed family
   densified, keeping the open triangles for the order-1 vertex, and
   regenerate the affected figures (Fig 11 and Fig 12 families, the zeta
   slices, and any wide-angle FK statement in the text).
2. Re-examine every published wide-angle FK claim against the corrected
   curve before quoting it again.

## Why there is no model difference: the same question in real space

A natural objection is that the draft's wide-angle FK must come from a
different or better treatment, since its Limber projection is analytic and
its measure collapses correctly. It does collapse correctly, and there is
no model difference: the two curves come from the same canoes Limber
branch, the same tree-level SPT bispectrum, the same ell_high_max = 1000,
the same measure factors. The draft's curve is this table folded through
the pipeline, not a separate analytic evaluation.

The mechanism becomes obvious once the angle is translated into the
separation the vertex actually depends on. The collapsed family
(1, cos gamma, cos gamma) is the driving-field three-point cumulant with
two points coincident and the third displaced TRANSVERSELY by
r_perp = chi * gamma on the same shell. In real space it is
<delta^2(x_a) delta(x_b)> at separation r_perp. At the mid shell
(chi = 6623 Mpc):

| gamma | r_perp at chi = 436 Mpc | 6623 Mpc | 8318 Mpc |
|---|---|---|---|
| 0.50' | 0.06 | 1.0 | 1.2 |
| 21.88' | 2.8 | 42 | 53 |
| 56.27' | 7.1 | 108 | 136 |
| 144.71' | 18 | 279 | 350 |
| 372.19' | 47 | 717 | 901 |

Matter correlations die beyond a few tens of Mpc, so a three-point cumulant
at 279 Mpc transverse separation is orders of magnitude below its value at
1 Mpc. The vertex MUST decay across the plotted range; that is not a
numerical artefact but the clustering of matter.

Now the grid. The deployed cosine spacing of 1/7 puts its first two
collapsed-family samples at c = 1 and c = 0.857, that is r_perp = 0 and
r_perp = 3584 Mpc at the mid shell. The grid's resolution in the variable
that matters is therefore about 3600 Mpc, while the physics happens between
1 and 100 Mpc. Every separation the paper plots falls inside that single
first interval, so the interpolation can only return something close to the
r_perp = 0 value.

Read this way the two facts fit together without any model difference:

* at gamma = 0.5' the queried configuration is r_perp ~ 1 Mpc, close enough
  to the coincident corner that the corner value is right to 0.15%, which
  is why the published small-angle amplitude survives
* at gamma = 145' the queried configuration is r_perp ~ 279 Mpc, but the
  table still returns the r_perp ~ 0 value, which is a physically different
  configuration by orders of magnitude

The query pattern itself was re-verified directly for this note: an
instrumented run of the production FK fold logs exactly one sorted cosine
triple per call, and the triples are (1, cos gamma, cos gamma) at the
external separation, with the three leg times bit-identical (equal-time
collapse). So the vertex is genuinely sampled at the observed separation,
and the r_perp reading above is the right one.
