# Route 1 references: convergence bispectrum measured on ray-traced lightcones

Collected 2026-08-25. Route 1 is the proposal in `emulator_candidates.md`
section 1d: predict the convergence bispectrum B_kappa with the same
machinery that builds the kappa3 vertex, and compare against a direct
measurement on simulations, which needs no fitting formula.

The literature turns out to be further along than the proposal assumed.
The comparison has already been published, on the simulation sets that are
public, and it INCLUDES squeezed configurations, which is the shape family
the FK vertex needs.

## The public simulation data

**Takahashi, Hamana, Shirasaki, Namikawa, Nishimichi, Osato, Shiroyama
2017**, "Full-sky Gravitational Lensing Simulation for Large-area Galaxy
Surveys and Cosmic Microwave Background Experiments", arXiv:1706.01472,
ApJ 850, 24.
108 full-sky multiple-lens-plane ray-tracing realisations, convergence and
shear maps from z = 0.05 to 5.3 in steps of 150 Mpc/h comoving (about
dz = 0.05 nearby), HEALPix Nside = 8192 and 4096 (0.48 and 0.96 arcmin
pixels), post-Born effects included. Public at
http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/

This is the better match for our purposes. Its shell spacing is a fixed
comoving interval, like our lambda-shell structure, and it reaches z = 5.3,
which brackets the z_s = 5 source plane the paper uses.

**Liu, Bird, Zorrilla Matilla, Hill, Haiman, Madhavacheril, Petri,
Spergel 2018**, "MassiveNuS: Cosmological Massive Neutrino Simulations",
arXiv:1711.10524, JCAP 03, 049.
101 cosmologies varying sum m_nu, Omega_m, A_s; 512 Mpc/h box, 1024^3
particles; 10,000 convergence map realisations per model, 12.25 deg^2 at
0.1 arcmin, source redshifts z_s = 1, 1.5, 2, 2.5 plus CMB lensing. Public
at http://columbialensing.org

## The measurements

**Namikawa, Bose, Bouchet, Takahashi, Taruya 2019**, "CMB lensing
bispectrum: assessing analytical predictions against full-sky lensing
simulations", arXiv:1812.10635, Phys. Rev. D 99, 063511.
Measures the convergence bispectrum in the Takahashi 2017 full-sky maps
and tests analytic predictions against it. Reports that tree-level
perturbation theory follows the simulation only to about ell = 200, and
one-loop to about ell = 600. That number is directly relevant to us, since
our vertex is built on a tree-level bispectrum and integrated well beyond
those multipoles.

**Sato, Hamana, Takahashi, Takada, Yoshida, Matsubara, Sugiyama 2009** and
**Kayo, Takada, Jain 2013**, arXiv:1207.6322.
1000 ray-traced convergence maps, 5 x 5 deg^2, source redshift z_s = 1.
The measured convergence bispectrum is valid to within 5% up to
ell = 4000.

## The comparison already exists, and it covers squeezed shapes

**Takahashi, Nishimichi, Namikawa, Taruya, Kayo, Osato, Kobayashi,
Shirasaki 2020** (the BiHalofit paper), arXiv:1911.07886, ApJ 895, 113,
section 6.
Section 6.1 compares BiHalofit, GM12 and SC01 against the Namikawa 2019
CMB-lensing convergence bispectrum in four configurations: equilateral
B_kappa(l, l, l), flattened B_kappa(l, l/2, l/2), **squeezed
B_kappa(l, l, 50)** and isosceles B_kappa(l, 1000, 1000), binned with
dl = 100 out to l = 2000.
Section 6.2 does the same for cosmic shear at z_s = 1 against the Kayo
2013 maps, in equilateral, flattened and **two squeezed configurations,
B_kappa(l, l, 84) and B_kappa(l, l, 510)**, binned dlog10 l = 0.13.

The squeezed panels are the useful ones here. They hold the soft leg fixed
at l_soft = 50, 84 or 510 and let the hard pair run, which is exactly the
integrand structure of our collapsed vertex: one soft leg selected by the
separation, two hard legs summed. So a published, simulation-measured
target already exists for the shape family the FK term integrates over.

## What this buys, and what it still does not

It tests, independently of any fitting formula, the chain from B_delta
through the Poisson factor, the screen-Hessian projection, the Limber
integral and the measure, on a quantity our machinery can produce.

It does not reach the collapsed corner itself. A measured B_kappa lives on
resolvable triangles with a finite soft leg, so it cannot probe the
coincident limit the FK vertex actually samples. What it can do is bound
the error at finite l_soft, and check whether our predicted squeezed
B_kappa matches the simulations over the l range where they overlap.

The Namikawa result is the sharpest warning already in print: tree-level
perturbation theory tracks the measured CMB-lensing bispectrum only to
l ~ 200. Our vertex integrates a tree-level bispectrum to l = 1000 as
deployed, and to 8192 in the converged variant.

## Suggested order of work

1. Reproduce the published squeezed comparison with our own projection
   code as a validation of the machinery: predict B_kappa(l, l, 84) and
   B_kappa(l, l, 510) at z_s = 1 with our Poisson factor, screen Hessian,
   Limber integral and measure, and check against Figure 13 of the
   BiHalofit paper. Agreement validates the transfer chain; disagreement
   localises a normalisation error.
2. Only then extend toward smaller l_soft, tracking where the measurement
   loses validity, to see how close to the collapsed corner a
   simulation-anchored statement can be pushed.
3. The remaining gap, the coincident limit itself, still requires a direct
   measurement of zeta on lightcones.
