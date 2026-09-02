TASK: Independently cross-validate the STF_lensing driving-field 3-cumulant by computing the
COLLAPSED (single-source-plane) cosmic 3PCF two ways and comparing. STAGE 0 = scalar convergence
3PCF self-check; STAGE 1 = shear 3PCF vs the external code `fastnc`. Fix z_s = 5. Start
LIGHTWEIGHT (a single angular separation gamma) before any sweep.

================================================================================
WHAT IS BEING VALIDATED
================================================================================
- The object under test is the equal-shell driving-field 3-cumulant zeta_abc(gamma, lambda)
  (appendix `eq: appendix equal shell approx`, "D7"). This task certifies zeta's ABSOLUTE
  normalization, radial measure, and spin projection by building a GENUINE 3-point function
  (convergence in STAGE 0, shear in STAGE 1) directly from zeta and comparing it, at a single
  source plane z_s=5, against an independent reference (a hand-rolled SPT-Limber kappa-3PCF in
  STAGE 0; the external code fastnc in STAGE 1).
- This is a self-contained 3PCF validation of zeta. It is a genuine THREE-point calculation
  (Order-1 K-vertex with the three legs fanning out to the external triangle) and has NOTHING to do
  with any 2-point observable; do not reference the FK <kappa kappa> diagram anywhere.

SCOPE NOTE: the equal-shell reduction is RADIAL (equal-time, single source shell); the ANGULAR
triangle is arbitrary — canoes compute_kappa3_zeta_table accepts general cosine_triples=(cos g12,
cos g23, cos g31). The COLLAPSED config (1, cos gamma, cos gamma) is the DEGENERATE-triangle limit
of the 3PCF: two of the three correlated sky points coincide and the third sits at angular
separation gamma, so the 3PCF reduces to a single angle gamma. We use this collapsed slice (z_s=5,
single gamma) per the project lead's directive, plus a non-degenerate isoceles cross-check
(STAGE-1 geometry note).

================================================================================
WHY THIS MATTERS
================================================================================
zeta's absolute normalization currently rests on a structural argument plus a STALE (2026-05-17,
pre-spin-2-helicity-fix) hand-rolled SPT reference. STAGE 0 refreshes the scalar SPT cross-check on
the present build (an INTERNAL different-quadrature/normalization check — same tree+Limber+flat-sky
family + same P(k), so NOT approximation-independent). STAGE 1 vs fastnc is the only GENUINELY
external check (fastnc shares nothing with our pipeline).

================================================================================
REPO + ENVIRONMENTS
================================================================================
- Project root: /Users/zzhang/Documents/MyDrafts/STF_lensing  (git, branch master)
- canoes: /Users/zzhang/projects/canoes (importable from the sft-wick env). PCAMBz0.txt present at
  /Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt (h-units: k [h/Mpc], P [(Mpc/h)^3]).
- sft-wick: /Users/zzhang/projects/SFT/sft-wick
- Our-side Python: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python (canoes+sft-wick+jax)
    Alt: /opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python (pyccl+canoes).
    Subagent Bash CANNOT `conda activate`; call interpreters by absolute path.
- fastnc: install into a FRESH dedicated venv (NOT the sft-wick/PyCCL envs — avoid perturbing the
  validated jax/canoes stack). Bridge results between envs via npz files.
- Project memories (ABSOLUTE path):
    /Users/zzhang/.claude/projects/-Users-zzhang-Documents-MyDrafts-STF-lensing/memory/
    Read: project_kappa3_spin2_helicity_bug_2026-06-04.md, project_kappa3_radial_measure_convention_2026-05-17.md.

================================================================================
READ FIRST (our formalism — refer to equations by \label; numbering may drift)
================================================================================
- sections/appendix.tex, the "Driving-Field Three-Point Cumulant" appendix:
    `eq: appendix bdelta tree`               B_delta = 2 F2 P P + 2 cyc  (tree matter bispectrum)
    `eq: appendix screen hessian`            screen-Hessian multipliers (spin-0 L^2/chi^2; spin-2 sqrt(L^2(L^2-2))/chi^2)
    `eq: appendix zeta squeezed`             zeta_XYZ = sum G P_{l3}(cos gamma) B^XYZ(lambda)
    `eq: appendix limber reduced bispectrum` B^XYZ = D(z)^4/chi^4 K^XYZ B_delta(l_i/chi)
    `eq: appendix modulus reconstruction`    zeta_B=<Phi00 Psi0 Psi0*>, zeta_D=<Psi0 Psi0 Psi0*>; <Psi+^2>,<Psix^2> from {zeta_PPP, zeta_B/D}
    `eq: appendix kappa3 cross shell`        folded lensing 3-point xi^(3)
    `eq: appendix equal shell approx`        kappa^(3) ~ zeta_abc(gamma,lambda) delta delta   (D7, the 3-cumulant)
  Real basis: a,b,c in {Phi00 (=T, spin-0, sources kappa), Psi+ = Re Psi0, Psix = Im Psi0 (=P,
  spin-2, sources shear gamma)}.
- sections/cosmology.tex: driving fields; shear vs kappa; the response propagator `eq: Resp Op cosmo`:
    R(n,lambda; n',lambda') = delta(n-n') * Theta(lambda'-lambda) * (Dbar(lambda)/Dbar(lambda'))^2,
    Dbar = a*chi, with the OBSERVED (smaller-lambda) point in the NUMERATOR. [Do NOT use the
    reciprocal.] kappa = -int_0^lambda Dsachs_1 dlambda' (cosmology.tex:923), kappa>0 <=> D<Dbar.
    c=1 geometric units: Phi,Psi dimensionless; A(a) = -3/2 Omega_m H0^2/a; H0 ~ 2.24e-4 Mpc^-1.
  The single-source-plane convergence efficiency to USE (cosmology.tex:1064-1068):
    K(lambda, lambda_s) = Dbar^2(lambda) * int_lambda^{lambda_s} dlambda1 / Dbar^2(lambda1)
                        = a^2(chi) * chi * (chi_s - chi) / chi_s.
- sections/path_int.tex: three-point diagrams; `eq: explicit resp op`.
- Invoke project skills `sft-wick` and `pyccl`.
- CRITICAL prior finding: a naive cross-shell integral int dlambda2 dlambda3 of the BARE kappa3
  BLOWS UP off-diagonal (FFTlog double-Bessel at unequal shells; this session: sigma3_GL = -2.6e3,
  unphysical). Use the equal-shell collapse (D7) + response folding; the responses damp the
  off-diagonal. Do NOT integrate the bare cross-shell.

================================================================================
STAGE 0 — scalar convergence 3PCF, collapsed, z_s=5  (clean scalar self-check; no fastnc)
================================================================================
OUR SIDE (canonical call recipe = scripts/sachs_sft/callables/kappa3_vertex/equal_time_limber/build_equal_time_limber_table.py):
- LOW : zeta_TTT via compute_kappa3_zeta_table(cosine_triples, lambda_shells, pk_matter_today, cosmo,
        ell_max=60, Nmax=128, spt_kind="tree_phi", radial_measure="lambda", units="physical",
        lambda_convention="project", channels=("TTT",)).
- HIGH: compute_kappa3_sigma3_high(..., ell_cut=60, ell_high_max=1000, n_ell=96, n_phi=64,
        radial_measure="lambda", units="physical", channels=("TTT",)).  <-- NOTE: HIGH does NOT
        accept spt_kind (its kernel is fixed tree-level). Combine: kappa3_combine_low_high(low, high).
- FOLD the collapsed convergence 3PCF in the AFFINE parameter lambda (NOT chi), single source plane z_s=5:
    Z_kappa(gamma) = int_0^{lambda_s} dlambda  K(lambda, lambda_s)^3  zeta_TTT(gamma, lambda),
    K as above. radial_measure="lambda" means zeta is a lambda-density 3-leg object: do NOT apply
    (dchi/dlambda)^3 (radial-measure memory; the per-leg (1+z)^4 already bakes in the Jacobian). If you
    integrate in chi instead, convert BOTH the kernel and zeta consistently (each leg picks up a^2(chi))
    and cross-check the two routes agree.
REFERENCE (independent IMPLEMENTATION, same tree+Limber+flat-sky family):
- Re-derive a hand-rolled SPT tree-level flat-sky Limber convergence 3PCF on the squeezed config at
  z_s=5: standard-SPT B_delta = 2 F2 P P + cyc with j0(L gamma)^3-type kernels (pure Limber, ~30-40
  lines, fast, no FFTlog blow-up). SAME cosmology + SAME PCAMBz0.txt P(k).
SIGN: the convergence 3PCF carries (-1)^3 from kappa = -int Dsachs_1. The +B_delta SPT j0^3 reference
  is sign-FLIPPED relative to the zeta fold (cf radial-measure memory: ref +1.15e-10 vs ours -9.66e-11).
  A single overall minus is EXPECTED — compare |magnitude| and note the sign explicitly.
COMPARE: single gamma first (e.g. gamma=1 arcmin) -> both numbers + |ratio| + sign, reported BEFORE
  anything else. Then a small gamma sweep (~6-10 pts, ~[0.5', few hundred ']) + plot.
GATE: PASS if |ratio-1| < ~0.5 (the standing SPT tree-level cross-check is ~30%, so 30% is expected
  agreement, NOT failure). STOP and debug before STAGE 1 only if disagreement exceeds a factor ~2 or
  the sign is wrong.

================================================================================
STAGE 1 — cosmic SHEAR 3PCF, collapsed, z_s=5  vs fastnc
================================================================================
OUR SIDE:
- The all-shear cumulant needs the spin-2 MODULUS channel (post-2026-06-04 helicity fix): dominated by
  zeta_D = <Psi0 |Psi0|^2> (the SUM Psi+^2+Psix^2, finite as gamma->0), with the un-conjugated zeta_PPP
  gamma^4-suppressed. Build the real-basis components <Psi+^3>, <Psi+ Psix^2> via the
  `eq: appendix modulus reconstruction` map from BOTH {zeta_PPP, zeta_D} (NOT zeta_D alone).
- HOW TO GET THE MODULUS CHANNELS (they are NOT valid `channels=` args — Channel literal is only
  TTT/TTP/TPP/PPP): EITHER (a) read zeta_Bmod, zeta_Dmod directly from the deployed npz
  .../equal_time_limber/equal_time_limber_kappa3_z5covgrid16_omega0316_h06711.npz (it already contains
  them; has_modulus_channels_flag=True), OR (b) call compute_kappa3_mod_zeta_equal_time (LOW,
  spt_kind="tree_phi") + compute_kappa3_mod_sigma3_high (HIGH, NO spt_kind) and sum, exactly as
  build_equal_time_limber_table.py lines ~157-198 do.
- Fold with the single-plane SHEAR response and PROJECT to the cosmic shear 3PCF natural components
  (Schneider-Lombardi Gamma^0..Gamma^3), flat-sky. Match fastnc's projection (see below).

STAGE-1 GEOMETRY (DO BOTH: collapsed primary + non-degenerate isoceles cross-check):
- fastnc parameterizes natural components as SAS with t1=t2 FORCED (FastNaturalComponents sets t2=t1)
  and an opening angle phi; third side t3 = 2 t1 sin(phi/2). Fix the two equal sides t1=t2=gamma and
  SWEEP the opening angle phi. This is ONE isoceles family our zeta matches exactly: evaluate zeta on
  cosine_triples = (cos t3, cos gamma, cos gamma) with t3 = 2 gamma sin(phi/2) (the base), i.e. two equal
  angular separations gamma and the third t3(phi). [canoes accepts general cosine_triples; the
  equal-shell collapse is RADIAL only, so this angular family is fully supported.]
- PRIMARY (collapsed): the phi->0 end (t3->0, cos_triple -> (1,cos gamma,cos gamma)) IS the collapsed
  config = two of the three 3PCF points coincident. It is a DEGENERATE triangle where shear
  natural-component spin phases are delicate; obtain it by EXTRAPOLATING the finite-phi curve to phi->0
  (do NOT evaluate fastnc exactly at
  phi=0). Test convergence in phi and in fastnc's Mmax (phi-resummation), Lmax_diag, epmu (default 1e-7).
- CROSS-CHECK (non-degenerate): also report the comparison at one or a few FINITE phi (well-posed, clean
  spin phases) on the SAME t1=t2=gamma family. This guarantees a delicate phi->0 limit cannot block the
  validation: if our zeta-derived Gamma^i and fastnc agree at finite phi and only diverge as phi->0, the
  issue is localized to the degenerate-limit numerics, NOT to zeta. The finite-phi agreement is the
  PRIMARY evidence that zeta is correct; the phi->0 extrapolation is the bonus that reaches the exact
  collapsed slice.

fastnc SIDE (install from scratch; commit d94fd2a / v1.2.3; method = Sugiyama+2024 arXiv:2407.01798):
- Install: `git clone https://github.com/git-sunao/fastnc`; `pip install .` into a FRESH venv. HARD
  dependency mpi4py is imported at import-time (coupling.py) even single-process: `brew install open-mpi`
  then `pip install mpi4py`. FFTLog/2DFFTLog are pure-Python (no compiled deps). First run builds a
  ~/.fastnc mode-coupling cache (minutes). Record the commit hash + the exact dep versions.
- INPUT MODEL — DO IT IN TWO STEPS (the second is the real validation):
  STEP 1 (quick run-check; NO subclass; built-in nonlinear): instantiate a SHIPPED bispectrum directly
    from cosmology + P(k) — BispectrumHalofit (BiHalofit) — set cosmology + the single source plane, run
    the FULL pipeline to Gamma^i. PURPOSE: de-risk the install (mpi4py, ~/.fastnc cache) and confirm the
    geometry/projection/units plumbing and ORDER OF MAGNITUDE before any subclass work. Compare to our
    side ONLY at LARGE gamma (>= ~100'), where tree-level ~ nonlinear to ~10-20% (appendix
    `zeta ell bands` Table: high-ell share falls to 6-26% there). Loose ~20% tolerance. This is a
    PLUMBING check, NOT the physics validation: at arcmin gamma tree vs nonlinear differ by ~2-10x
    (BiHalofit is nonlinear; our zeta is tree-level), so do NOT compare Step-1 BiHalofit at small gamma.
  STEP 2 (the clean validation; tree-level subclass; ALL gamma): there is NO entry point for a precomputed
    B_kappa(l1,l2,l3). You supply a MATTER bispectrum by subclassing fastnc.bispectrum.BispectrumBase and
    implementing matter_bispectrum_no_baryon(k1,k2,k3,z) [k in h/Mpc]; fastnc then applies its OWN built-in
    single-plane convergence Limber kernel (g = (3/2)(H0/c)^2 Om0 (1-chi_l/chi_s)) and 2DFFTLog. fastnc
    ships NO tree-level class (only BiHalofit, Gil-Marin, NFW); its internal F2_eff is the Gil-Marin
    EFFECTIVE F2 — do NOT reuse it. WRITE a ~30-line subclass implementing standard SPT:
    B_delta = 2 F2(k1,k2) P(k1) P(k2) + 2 cyc with the true
    F2 = 5/7 + (k1.k2/2)(1/k1^2 + 1/k2^2) + (2/7)(k1.k2)^2/(k1^2 k2^2), feeding the SAME PCAMBz0.txt P(k)
    via set_pklin/set_lgr. This MUST match canoes' tree-level SPT F2 exactly (NOT F2_eff), and gives a
    clean ratio at ALL gamma incl. arcmin. This is the run that validates zeta.
- Recipe (mirror the tutorial): subclass.set_pklin/set_lgr; set_cosmology(astropy wCDM with
  Omega_m=0.3160919980475834, h=0.6711, n_s=0.97); set_source_distribution([z_s],[1.0]) (single delta
  plane); compute_kernel(); interpolate(); decompose(); FastNaturalComponents.set_bispectrum(bs);
  then evaluate Gamma^i. Cosmology/P(k) MUST match our side exactly.
- z_s=5 CAVEAT: set_cosmology hardcodes the chi<->z spline on z in [0,5] (bispectrum.py:292), so z_s=5 is
  at the grid edge. Verify no extrapolation artifacts (check the kernel/los grids near the top); if
  flaky, patch the local clone to linspace(0, 5.5, ...) or run z_s slightly below 5 and document.
- PROJECTION: config key `projection` in {x (FFT-native default), cent (centroid, arXiv:2309.08601),
  ortho (orthocenter)}. Pick cent OR ortho directly from the native x-projection to MATCH our side; AVOID
  the cent<->ortho round-trip (ortho2cent is flagged unvalidated in-source). All Gamma^i are flat-sky.
COMPARE: at fixed gamma, do ONE clean FINITE-phi point FIRST (well-posed: our zeta at (cos t3,cos g,cos g)
  vs fastnc t1=t2=g at that phi) -> numbers + |ratio| + sign, reported before any sweep. Then the
  phi-sweep down toward phi->0 (the collapsed extrapolation), then a small gamma sweep. Plot Gamma^i vs
  phi (showing the phi->0 limit) and vs gamma.

================================================================================
CONVENTION CHECKLIST (verify EVERY item; an h^3/h^6 slip is the likeliest OOM cause)
================================================================================
[ ] Same cosmology + same PCAMBz0.txt P(k) + same z_s=5 on both sides and both stages.
[ ] Same bispectrum in the VALIDATION (Step 2): standard-SPT tree-level everywhere (our zeta tree_phi;
    fastnc your SPT subclass, NOT F2_eff). Step-1 BiHalofit is nonlinear and is ONLY the large-gamma
    (>=~100') plumbing check (loose ~20%); never compare Step-1 BiHalofit at arcmin.
[ ] Units: assert the loaded zeta npz cosmo_meta units == "physical" (canoes builds in h-units;
    units="physical" applies h^4 to the equal-shell zeta [2026-06-09 dimensional correction:
    native (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, was h^6] and h^-1 to chi/lambda). PCAMBz0.txt is h-units (k[h/Mpc],
    P[(Mpc/h)^3]); carry h=0.6711 through every chi, k, and arcmin->rad (gamma[rad]=gamma[arcmin]*pi/10800)
    conversion; print the unit of every intermediate.
[ ] Collapsed config matched: (1, cos gamma, cos gamma) on our side <-> fastnc t1=t2=gamma, phi->0.
[ ] Spin-2: helicity-corrected modulus (zeta_D dominant) + full reconstruction from {zeta_PPP, zeta_D};
    natural-component basis (Gamma^0..3) and projection (cent or ortho) matched to fastnc.
[ ] Sign: convergence 3PCF carries (-1)^3; expect a sign flip vs the +B_delta SPT reference. Compare
    |magnitude|, note the sign.
[ ] Flat-sky both sides.

================================================================================
DELIVERABLES (write everything under scripts/sachs_sft/analyses/<new-folder>/)
================================================================================
1. STAGE 0 lightweight single-gamma result (our vs SPT-Limber: numbers, |ratio|, sign) reported FIRST.
   Then the STAGE 0 gamma sweep + plot.
2. STAGE 1, in TWO steps. Step 1 (plumbing): the BiHalofit large-gamma (>=~100') run-check
   (order-of-magnitude agreement, loose ~20%), reported as confirmation the fastnc pipeline runs and the
   geometry/units are right -- BEFORE any subclass work. Step 2 (validation, tree-level subclass):
   lightweight one (gamma, finite phi) point FIRST (our vs fastnc: numbers, |ratio|, sign), reported
   before any sweep; then BOTH (a) the phi-sweep -> phi->0 collapsed extrapolation (primary slice) and
   (b) the finite-phi non-degenerate isoceles cross-check, with plots of Gamma^i vs phi (collapsed limit)
   and vs gamma.
3. Plots (our vs reference + ratio panel) -> <new-folder>/figures/ (PDF) and <new-folder>/outputs/ (PNG +
   npz of the numbers). Do NOT write to the project-root figures/ (reserved for paper \includegraphics).
4. <new-folder>/README.md: matched conventions, both single-gamma results, both sweep ratios, residual
   sizes, the phi->0 handling, the fastnc commit + deps, and a verdict (does each reference confirm
   the 3-cumulant zeta, to what %).
5. Reproducible scripts per path (STAGE-0 ours, STAGE-0 SPT reference, STAGE-1 ours, STAGE-1 fastnc),
   each with the exact env interpreter path in the header.

================================================================================
CONSTRAINTS
================================================================================
- English in code/comments/docs. Do NOT edit sections/*.tex or biblio.bib. Do NOT modify canoes or
  sft-wick source (consume them). Archive, do not delete, prior artifacts.
- Prefer analytical/transparent reductions; document every convention you assume.
- Report failures honestly: if the phi->0 limit (or any convention) cannot be matched, STOP and report
  rather than compare mismatched configs.
- Mandatory gates: the STAGE-0 |ratio-1|<~0.5 gate, and single-gamma-FIRST before any sweep.
