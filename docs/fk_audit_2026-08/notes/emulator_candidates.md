# Verified 3-point emulator landscape (survey 2026-08-25)

Method: three parallel literature sweeps (Euclid-specific; general matter
bispectrum; potential/convergence level), 39 raw candidates, per-candidate
adversarial verification against arXiv full text and code repositories, plus a
completeness critic pass. Every entry below was source-verified.

## Genuine 3-point emulators

### 1. Euclid/BACCO baryonic matter-bispectrum boost emulator  <- THE relevant one

- Paper: "Euclid: An emulator for baryonic effects on the matter bispectrum",
  P. A. Burger et al. (Euclid Consortium), arXiv:2506.18974, A&A 705, A170.
- Emulates: R(k1,k2,k3; a) = B_hydro / B_gravity-only of the 3D COLD-matter
  density field, equal-time, arbitrary triangles (k values are NN inputs);
  companion P(k) boost co-emulated on the same footing.
- Domain (implemented module): k in [0.017, 17.655] h/Mpc; a in [0.24, 1.01]
  (z from 0 to ~3.2). Params: omega_cold, omega_baryon, sigma8_cold + 5
  baryonification params (log10 Mc, M1, beta, eta, theta_inn).
- Training: CosmoPower deep NN (4x32), 1200 Latin-hypercube nodes of BACCO
  high-res N-body + baryonification + cosmology rescaling; validated vs
  FLAMINGO hydro; <2% at 68th percentile.
- Code: `pip install baccoemu`, class `Matter_bispectrum()`
  (Bitbucket rangulo/baccoemu; docs baccoemu.readthedocs.io).
- Caveats: boost ONLY (absolute B needs BiHalofit underneath, whose ~10%
  systematic propagates); cold-matter parameter conventions (sigma8_cold vs
  total-matter, relevant with massive neutrinos); baryonified N-body, not
  full hydro, as training truth.
- Status for us: the baryonic-correction factor of the recommended hybrid
  input model. The ONLY Euclid-affiliated 3-point emulator; EuclidEmulator3
  does not exist (confirmed).

### 1b. What BiHalofit is, and how independent of SPT it is

Since it is the gravity-only backbone of every option here, its status
matters. BiHalofit (Takahashi et al. 2020) is a FITTING FORMULA for the
nonlinear matter bispectrum, calibrated on 41 wCDM N-body models, quoted
at about 10% for k < 3 h/Mpc and 15% for k < 10, over z = 0 to 3.

It is not independent of perturbation theory. Its structure is halo-model
inspired with an SPT core:

* the three-halo term is 2 (F2_tree + d_n q) PE_1 PE_2 + cyclic, times a
  damping factor. It uses the tree-level F2 kernel directly, with the
  linear spectrum replaced by a fitted effective spectrum PE and a small
  additive correction
* the one-halo term is a pure fit, a product over legs of
  1 / (a_n q^alpha + b_n q^beta) / (1 + 1/(c_n q)), which vanishes on
  large scales

So on large scales it reduces to the tree-level structure, and the genuine
new information is at small scales, where the one-halo term is calibrated
on simulations. Measured here, the reduction is not exact: BiHalofit sits
2 to 4% above the exact tree assembly for every k below 1e-3 h/Mpc, which
is inside its quoted accuracy but is the reason the blend seam in
`code/b_model` is placed at the crossover rather than deeper.

The practical consequence for this project: comparing tree against
BiHalofit measures how much the FK term grows once small-scale nonlinear
power is included. It is NOT an independent theoretical check of the
tree-level calculation. An independent check would be either a direct
measurement of zeta on N-body lightcones, which nobody has done, or the
Monte-Carlo cross-check, which `finding_mc_crosscheck.md` shows validates
the fold rather than the input.

### 2. Lensing-observable-level surrogates (methodology references only)

All are NN surrogates of BiHalofit-based theory codes for projected,
survey-specific statistics; none outputs a 3D B(k1,k2,k3,z); none is a
candidate input for the SFT driving-field vertex.

- i3PCF (Gong, Halder et al., arXiv:2304.01187, JCAP): integrated shear 3PCF
  zeta_+/- kernels on 100 ell bins in [2, 15000], z in [0, 2.1]; wCDM +
  HMcode c_min; theory pipeline public (github.com/anikhalder/i3PCF), trained
  weights not.
- KiDS-1000 <Map^3> + COSEBIs (Burger et al., arXiv:2309.08602): CosmoPower
  NN over {Om, S8, A_IA, dz}; apertures 4'-32'; BiHalofit + HMcode2020
  training data (github.com/sheydenreich/threepoint).
- DES-Y3 <Map^3> via fastnc (Gomes et al., arXiv:2503.03964; data paper
  2508.14018): CosmoPower NN, 5-param Sobol training, <0.3% NN error;
  fastnc itself is public and pip-installable (github.com/git-sunao/fastnc),
  the trained weights are not.
- HSC-Y3 <Map^3> (Sugiyama, Gomes & Jain, arXiv:2508.14019): emulates
  kernel-stripped Map^3(R, z) on 8 radii x 130 z nodes (z up to 3), 6-param
  wCDM box; fastnc public, analysis repo 404s.
- KiDS-Legacy <Map^3>_z (Linke et al., arXiv:2606.12389): redshift-resolved
  integrand emulator (arbitrary source n(z) reusable), JAX framework
  cosmoemuJAX (github.com/pburger112/cosmoemuJAX) + threepoint_py
  (github.com/llinke1/threepoint_py, v1.0.0, Zenodo). The most reusable
  design pattern of the family; 7-param box includes BCM Mc, eta.
- Design pattern worth copying: KiDS-Legacy/HSC decouple the emulator from
  the source n(z) by emulating the per-lens-redshift integrand; the SFT
  analogue would be emulating zeta_XYZ(gamma, chi) per shell.

### 3. Simulation-trained but unreleased / other fields

- MassiveNuS tomographic convergence bispectrum GP (Coulton et al. 2019,
  arXiv:1810.02374): the one genuinely simulation-trained lensing-bispectrum
  surrogate; ell in [150, 3150], 3 params (Mnu, Om, As), 101 cosmologies;
  ~0.5% interpolation claim; NEVER released (maps public at
  columbialensing.org, emulator internal).
- 21-cm EoR bispectrum ANN "EmuPBk" (Tiwari et al., arXiv:2108.07279):
  wrong field (21-cm), single z = 8, 3 astro params; repo dormant.

## Closest relatives (NOT emulators)

- BiHalofit (arXiv:1911.07886): nonlinear matter bispectrum FITTING FORMULA,
  41 wCDM N-body models, ~10% (15%) at k < 3 (10) h/Mpc, z = 0-3; public
  C/Fortran/Python. The de facto standard; underlies every <Map^3> emulator
  above. -> the gravity-only base of our recommended hybrid.
- SC01 (astro-ph/0009427), GM12 (arXiv:1111.4477), three-shape halo model
  (arXiv:1511.02022): older fitting formulas, superseded by BiHalofit.
- GEO-FPT (arXiv:2303.15510): galaxy (redshift-space) bispectrum fitted
  model; wrong tracer.
- COMET (arXiv:2208.01070): emulates one-loop P(k) multipoles; its bispectrum
  output is direct tree-level PT, not emulated.
- Field-level emulators (Jamieson et al. 2206.04594, MG-NECOLA 2510.20086,
  Learning-the-Universe 2502.13242, flow-matching kappa maps 2605.23114):
  emulate the field; bispectrum only a validation metric. Measuring B on
  emulated maps is the fallback route to a "nonstandard" bispectrum, at high
  cost and noise.
- SBI pipelines (SimBIG 2310.15243, SBi3PCF 2510.13805): normalizing-flow
  posteriors FROM 3-point data, not surrogates OF the statistic.

## Confirmed negative results

- No absolute gravity-only matter bispectrum emulator (any collaboration).
- No matter/galaxy 3PCF emulator (Aemulus VI "beyond-standard" stops at
  2-point-ish statistics).
- No Weyl/gravitational-potential bispectrum emulator; no CMB-lensing
  bispectrum emulator (analytic + fitting formulas only, ~10% accuracy noted
  as the limitation by 2507.20262).
- No Euclid HOS emulator: HOWLS (2301.12890) and DR1 HOS (2510.04953) are
  Fisher forecasts; the convergence bispectrum is not among their statistics.
- Euclid spectroscopic P+B pipeline is analytic PT (2603.27966, 2605.21436);
  OU-LE3 3-point products are estimators, not surrogates (2605.03012).
- No symbolic-regression bispectrum (syren family is P(k) only; the 2025
  review 2510.18749 lists the bispectrum as future work).
- Suites with no 3-point product: DarkQuest I/II, Aemulus, EuclidEmulator1/2,
  CSST/Kun, FORGE, e-MANTIS, Mira-Titan/CosmicEmu, Matryoshka, PyBird-JAX,
  COBRA (PyBird-JAX and CSST explicitly list bispectrum as future work).
- Completeness critic re-checked (a)-(f) angles (Euclid 2024-2026, baccoemu
  docs, EE3, ML4PS, CosmoPower roster, halo-model hybrids): nothing
  substantive missed.

## Vault status (Obsidian)

The research vault holds no note on any matter-bispectrum emulator
(BiHalofit/EuclidEmulator/bacco: zero hits). Nearest holdings: COBRA one-loop
galaxy bispectrum factorization [[bakxOneLoopGalaxyBispectrum2025d]], WL
bispectrum information content (Kayo/Takada), [[zhangGeneralPolynomialEmulator2026]]
(MomentEmu), and the 57-paper collection
`ResearchLib/collections/Search-CosmologyBispectrumPDF.md`. Candidates for a
future paper-ingest: 2506.18974, 1911.07886, 1810.02374, 2304.01187.

### 1c. What a tree-versus-BiHalofit comparison does and does not establish

Asked directly: which part of the SPT model does the cross-comparison
validate, and can it support "our model is not far off"?

**It validates nothing about correctness.** Both models are evaluated
through the same `b_delta_fn` interface, the same quadrature, the same
shells, the same Poisson factor, the same screen-Hessian projection, the
same measure and the same spin-channel algebra. Any error in that chain
appears identically on both sides and cancels in the ratio. Worse, the
shared structure extends into the bispectrum itself: BiHalofit's three-halo
term uses the tree-level F2 kernel, so even a mistake in the mode-coupling
kernel would be common to both. The comparison is a sensitivity
measurement, not a check.

**It does support one bounded statement, under a condition that the
deployed setup happens to satisfy.** The hard legs reach k = ell / chi_h,
so at ell_max = 1000:

| shell z | chi_h [Mpc/h] | k at ell_max = 1000 | k at ell_max = 8192 |
|---|---|---|---|
| 0.10 | 292.6 | 3.42 | 28.00 |
| 1.00 | 2291.0 | 0.44 | 3.58 |
| 3.10 | 4445.0 | 0.22 | 1.84 |
| 5.70 | 5582.2 | 0.18 | 1.47 |

At the deployed cutoff the vertex never leaves BiHalofit's calibrated range
(about 10% for k < 3 h/Mpc). Since BiHalofit is N-body calibrated there,
the true gravity-only bispectrum sits within roughly 10% of it, and tree
level sits a factor 1.66 away in the folded observable. So the honest
version of "not far off" is:

> switching from tree level to an N-body-calibrated bispectrum moves the
> FK amplitude by a factor 1.66, not by orders of magnitude, and the truth
> should lie near the BiHalofit end with tree level as a lower bound.

**The condition fails once the cutoff is raised.** At ell_max = 8192 the
nearest shell reaches k = 28 h/Mpc, nearly three times beyond BiHalofit's
calibration, where baryons also dominate. So the bounding argument weakens
exactly as the integral approaches convergence. That is uncomfortable but
it is what the numbers say: the closer to converged, the less controlled
the model error.

**What actually supports "not far off" is elsewhere**, and it covers
precisely what the tree-versus-BiHalofit ratio cancels:

* Order-0 against PyCCL, an independent code, validates the whole
  projection chain (Poisson factor, screen Hessian, Limber, measure) at
  0.88 to 1.10 across four observables
* the fastnc closure test validates the three-point machinery against an
  independent implementation
* the Mathematica/xAct derivations validate the transfer formulas
* the Monte-Carlo cross-check validates the diagram assembly and the fold

**The gap.** None of those covers the collapsed slice the FK term actually
uses. The fastnc check runs on OPEN triangles, since fastnc cannot reach
the collapsed corner; `fk_analytic.py` states in its own docstring that it
does not independently re-derive the zeta table's internal units, because
both sides read the same table; and the Order-0 comparison is a two-point
check that never touches the three-point vertex. So the projection chain,
the diagram assembly and open-triangle zeta all have independent support,
while the absolute normalisation of the collapsed slice does not. The only
route to closing that gap is a direct measurement of zeta on N-body
lightcones, which nobody has done.

### 1d. What an emulator could be compared against, and the one route that closes the gap

Follow-up question: if BiHalofit cannot validate the model, can the
emulator?

**baccoemu cannot either, for the same reason.** It emulates the RATIO
B_hydro / B_gravity-only and must be multiplied onto a gravity-only model,
which in practice is BiHalofit. The gravity backbone, and with it the
shared F2 kernel and the whole projection chain, is unchanged. Comparing
with it is again a sensitivity measurement.

What it does buy is a new systematic axis that neither tree nor BiHalofit
contains: the baryonic response, calibrated on baryonification and
validated against FLAMINGO hydrodynamics. That axis lands exactly where it
is most needed, since the converged vertex is dominated by k of ten to a
hundred h/Mpc, which is where baryons dominate. So baccoemu quantifies HOW
UNTRUSTWORTHY a converged FK amplitude is, rather than making it
trustworthy. Not yet measured here: baccoemu is installed nowhere.

**There is, however, a genuinely independent comparison available, and it
uses simulations rather than an emulator.** The pipeline maps B_delta
through the Poisson factor, the screen-Hessian projection, the Limber
integral and the measure bookkeeping. That same chain can produce the
CONVERGENCE bispectrum B_kappa(l1, l2, l3), which is measurable directly on
ray-traced lightcones without any fitting formula. Public data exist: the
MassiveNuS convergence maps (Columbia Lensing; the Coulton et al. 2019
emulator was never released but the maps are public) and the Takahashi
et al. 2017 full-sky lensing simulations.

Predicting B_kappa with our machinery and comparing against a direct
measurement would independently test:

* the Poisson factor and screen-Hessian normalisation
* the Limber projection and the measure bookkeeping
* the accuracy of the input B_delta over the relevant k range

What it would NOT test is the collapsed slice itself, since a measured
B_kappa lives on resolvable triangle configurations and cannot reach the
phi to zero corner. That is the same structural limitation the fastnc
check ran into. Still, it is closer to the gap than anything currently
available: Order-0 is a two-point check, fastnc tests open-triangle
machinery, and the Monte-Carlo tests the fold.

**The only route that closes the gap fully** is a direct measurement of
zeta on N-body lightcones: construct Phi_00 and Psi_0 from a simulated
density field and measure their three-point cumulant in the collapsed
configuration. Nobody has done this. The project already has the
infrastructure to start from: `/Users/zzhang/Workspace/stf-transfer`
computes driving-field two-point spectra C_l^AB from gevolution metric
perturbations using the same symbolic transfer catalogue, so extending it
to the three-point cumulant is the natural continuation.
