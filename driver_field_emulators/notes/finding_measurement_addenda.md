# Measurement addenda: three findings not covered by the other notes

Collected 2026-08-25 while auditing whether everything measured had been
written down. Each of these came out of a check run for another purpose,
and each is worth keeping.

## 1. The deployed build's own quadrature under-resolves at wide angles

Octave banding was introduced to make a high cutoff affordable. Summing
bands over (60, 1000] should reproduce a single direct run over the same
window, since the windows are disjoint in max-leg. It does not, exactly,
and the difference is informative:

| triple | shell | single run | band sum | difference |
|---|---|---|---|---|
| (1, 0.99987, 0.99987), gamma = 56' | 396.6 Mpc | -1.6297e-19 | -1.6332e-19 | 0.21% |
| same | 2326.6 Mpc | -2.1291e-16 | -2.1424e-16 | 0.63% |
| (1, 0.90194, 0.90194), gamma = 1535' | 2326.6 Mpc | +2.6887e-16 | +2.4152e-16 | **10.2%** |
| same | 396.6 Mpc | +5.7160e-21 | +2.4675e-21 | 57% (negligible magnitude) |

The banded result is the more accurate one: each band carries its own
96-node rule, so the node density in ln(ell) is uniform, whereas the single
run spreads 96 LINEARLY mapped nodes over the whole interval and therefore
under-samples the low-ell decade. At small separations the two agree to
under a percent, but at gamma = 1535 arcmin the deployed quadrature is off
by about 10%.

This is a property of the deployed build, independent of the cosine grid
and of the cutoff. It does not change any conclusion here, since the
conclusions rest on separations below 600 arcmin, but it means wide-angle
vertex values from a single-run build carry a percent-to-ten-percent
quadrature error on top of everything else.

## 2. What the nonlinear bispectrum does at the vertex, before folding

Measured on the canoes HIGH branch with the same quadrature and shells,
tree against BiHalofit, for zeta_TTT:

* at the nearest shell the vertex is amplified by 20 to 46 times
* some cosine triples CHANGE SIGN between the two models
* at the farthest shell the ratio is 1.00 to within 2%

The sign changes are not a defect. zeta is a band-cancelling sum: the
archived ell-band table shows the low-ell and high-ell bands carrying
opposite signs at intermediate angles, so a model that boosts the high-ell
band strongly can overturn the total. This is why a nonlinear-versus-tree
"ratio" is not bounded below by one, and why the integration test in
`code/b_model/tests/test_vertex_integration.py` asserts on amplification of
the magnitude rather than on the signed ratio.

It also sets the expectation for the fold: a 20 to 46 times amplification
at the near shell turned into only 1.66 in xi_FK, because the near shells
carry little lensing weight.

## 3. An upstream bug in fastnc 1.2.3

`Halofit.__init__` calls `set_pklin` before `set_cosmology`, and
`set_pklin` runs the sigma8 renormalisation, which reads `self.cosmo`.
Constructing with all arguments at once therefore raises
`AttributeError: 'Halofit' object has no attribute 'cosmo'`.

The workaround, in `code/b_model/halofit_backend.py`, is to construct empty
(the renormalisation short-circuits on `self.k is None`) and then set
cosmology, growth and P(k) in dependency order. Worth reporting upstream if
anyone is in touch with the fastnc author; the fix is a two-line reorder.

A second fastnc detail that matters for correctness here: its
`get_interpolated_pklin` EXTRAPOLATES above the table's k_max via the
spline rather than truncating to zero. Measured slope beyond k = 206 h/Mpc
is about -2.9, which is physically sensible, so BiHalofit remains usable
above the table limit where the canoes tree branch fails loud. The two
branches therefore have different high-k behaviour, which is why a
like-for-like tree-versus-nonlinear comparison must cap ell_max at the
tree branch's limit.
