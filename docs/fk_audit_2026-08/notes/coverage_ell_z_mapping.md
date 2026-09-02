# What (ell, z) coverage do we actually need from a bispectrum emulator?

The pipeline never asks an emulator for an "angular scale" directly. The
Limber-reduced vertex evaluates B_delta at k_i = ell_i / chi(z) per shell, so
the emulator's (k, z) training domain induces the (ell, z) coverage. Fiducial
cosmology (Om = 0.3161, h = 0.6711):

| shell z | chi [Mpc] | ell=10          | ell=60          | ell=1500        | ell=5000        |
|---------|-----------|-----------------|-----------------|-----------------|-----------------|
| 1       | 3414      | 4.4e-3 h/Mpc    | 2.6e-2 h/Mpc    | 0.65 h/Mpc      | 2.2 h/Mpc       |
| 2       | 5330      | 2.8e-3 h/Mpc    | 1.7e-2 h/Mpc    | 0.42 h/Mpc      | 1.4 h/Mpc       |
| 3       | 6527      | 2.3e-3 h/Mpc    | 1.4e-2 h/Mpc    | 0.34 h/Mpc      | 1.1 h/Mpc       |
| 5       | 7971      | 1.9e-3 h/Mpc    | 1.1e-2 h/Mpc    | 0.28 h/Mpc      | 0.93 h/Mpc      |

## The (z=5, ell=10) corner: soft leg tree-exact, hard pair NOT

The wide-angle vertex is a squeezed object with one soft and two hard legs,
and the two must be kept separate.

**Soft leg.** The FK collapsed configuration is (n_a, n_a, n_b) with
n_a . n_b = cos gamma: the leg at angle gamma carries the Legendre factor
P_{l3}(cos gamma), so at gamma >~ 2 deg only l3 <~ 60 survives. At the z = 5
shell that is k_soft <~ 1.1e-2 h/Mpc (k = 1.9e-3 h/Mpc at l3 = 10), deep in
the linear regime: tree-level/linear theory is exact for THIS leg, and every
squeezed-limit model uses P_lin(k_soft) for it anyway, so emulator k_min
coverage of the soft leg is a non-issue.

**Hard pair.** The two COINCIDENT legs have relative angle zero, so
P_{l1}(1) = 1: the pair sum over l1 ~= l2 (forced near-degenerate by the
triangle inequality |l1 - l2| <= l3) runs un-damped up to the build cutoff
(currently ell = 1000 in the HIGH/Limber branch). These are hard modes,
k_hard = l/chi(z_shell), integrated at the vertex like
INT dl l B(k_soft, l/chi, l/chi): the nonlinear response of small-scale
power to the long mode. Tree-level SPT does NOT describe them at the shells
where the lensing kernel peaks. Leading squeezed boost factor
(B_sq ~ R_1(k_h) P_lin(k_s) P(k_h), so boost ~ P_NL/P_L at the hard leg;
pyccl halofit, fiducial cosmology with sigma8 = 0.810 pinned by direct
top-hat integration of PCAMB_pyccl_stf_fid_z0.txt, whose header also states
0.81):

| shell z | chi [Mpc] | l=60 | l=250 | l=500 | l=1000 | l=2000 | l=5000 |
|---|---|---|---|---|---|---|---|
| 0.1 | 436 | 1.18 | 4.42 | 10.5 | 20.4 | 31.9 | 49.6 |
| 0.3 | 1242 | 0.99 | 1.40 | 2.56 | 5.98 | 13.3 | 26.9 |
| 0.5 | 1959 | 0.99 | 1.10 | 1.55 | 3.02 | 7.04 | 18.1 |
| 1.0 | 3414 | 0.99 | 1.01 | 1.11 | 1.49 | 2.74 | 8.01 |
| 2.0 | 5329 | 1.00 | 1.00 | 1.02 | 1.10 | 1.39 | 2.88 |
| 3.0 | 6526 | 1.00 | 1.00 | 1.01 | 1.05 | 1.16 | 1.72 |
| 5.0 | 7969 | 1.00 | 1.00 | 1.01 | 1.02 | 1.05 | 1.18 |
| 5.7 | 8315 | 1.00 | 1.00 | 1.01 | 1.02 | 1.04 | 1.13 |

**Consequences for the wide-angle FK amplitude.**
1. The wide-angle vertex value at a given shell splits into fully-soft
   triples (all three ell <= 60) and squeezed high-ell bands. The archived
   ell-band table (appendix draft Table II, z = 5.7 shell) gives the
   high-ell share of zeta_TTT as 40% at gamma = 42', 26% at 117' (~2 deg),
   6% at 326', AT THAT SHELL; mid-z shells reach larger boosts at the same
   ell (table above), so the kernel-weighted effect on the FK fold is
   plausibly tens of percent and is NOT tree-exact. This is the same channel
   behind the intro's statement that the tree-level FK amplitude is
   conservative; note however that the bands carry mixed signs at
   intermediate angles (Table II), so the sign and size of the wide-angle
   shift await the numerical probe of point 3.
2. The hard-pair sum's ell_max convergence must be re-examined, for tree AND
   boosted. A toy diagnostic with the faithful squeezed scaling (hard-pair
   sum ~ INT dl l P(l/chi)) gives S(4000)/S(1000) = 1.44 / 1.63 / 1.87
   (tree) and 3.25 / 2.70 / 2.2 (halofit-boosted) at z = 0.5 / 1 / 2: under
   this monotone proxy neither is converged at the current
   ell_high_max = 1000. The proxy lacks the Gaunt structure and the band
   sign cancellations of the real sum, so it mandates rather than replaces
   an explicit ell_high_max sweep ({1000, 1500, 2000, 3000, 4000}) on the
   (1, cos gamma, cos gamma) family.
3. Decisive numerical probe (not yet run): rebuild zeta(gamma, chi) at a few
   wide gamma with the hard legs boosted (halofit-ratio or response model)
   and compare the FK fold against tree; band-decompose per shell as in the
   archived Table II.

What remains true from the naive reading: the OBSERVABLE'S soft angular
scale never needs emulator coverage at its own k = ell_obs/chi, and the
fully-soft triple population is tree-safe except at the lowest shell
(z = 0.1, chi = 436 Mpc: ell = 60 reaches k = 0.138 /Mpc = 0.21 h/Mpc with
an 18% one-leg boost, carried at ~0.5x the peak of the a^2 chi
(chi_s - chi)/chi_s kernel weight; the ell_cut = 60 seam probe in the
roadmap decides whether the LOW branch keeps tree-level there).

## Where emulator coverage actually bites

1. High ell at each shell: ell = 1500-5000 maps to k ~ 0.3-2.2 h/Mpc for
   shells z >= 1 (within typical emulator domains, k_max ~ 2-10 h/Mpc), but
   for NEARBY shells (chi ~ 50-500 Mpc) k = ell/chi reaches 10-100 h/Mpc,
   beyond any training domain. Those shells are strongly lensing-suppressed,
   so blending to BiHalofit/tree-level outside the domain is benign, but the
   blend must be smooth (no hard clip; see the power-conservation note).

2. Redshift ceiling: the vertex needs every shell z in [0, z_s] with z_s up
   to 5. Sim-trained emulators typically stop at z ~ 1.5-3. Beyond the
   emulator's z_max the fallback is tree-level, and the fallback error shrinks
   exactly where it is used: at fixed k the nonlinear-to-tree ratio tends to 1
   as z grows (growth suppression). At z = 3-5, even the highest-k calls
   (~0.3-1.1 h/Mpc at ell = 1500-5000) are only mildly nonlinear.

3. Squeezed configurations: the FK collapsed geometry weights ell2 ~= ell3
   with ell1 free, so triangle-shape coverage (squeezed corner) matters as
   much as the (k, z) box. Check each candidate's training-triangle sampling.

## Recommended hybrid evaluation scheme

    B(k1, k2, k3; z) =
      B_emu                      inside the training box (k and z)
      B_tree                     at low k (exact there)
      BiHalofit or B_tree        at high k / high z outside the box

with smooth (e.g. log-k window) blending at the seams, and boundary tests at
both seams per the boundary-validation methodology (probe at seam +/- delta,
including extreme shells: the z ~ 0.1 shell and the z = 5 shell).
