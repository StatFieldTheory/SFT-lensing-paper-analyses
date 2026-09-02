# Converting an emulator's matter bispectrum into SFT driving-field inputs

> Overview note. The implementation-level, source-verified specification
> (exact formulas, file:line, quadrature, query set, grid requirements) is
> `input_ready_b_to_zeta_spec.md` - code from that one.

Sources of truth: `sections/appendix.tex` (append: driving-field spectra),
`sections/cosmology.tex` (subsec: driving fields), and the deployed callable
`SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber/`.

## Why no "Ricci/Weyl field emulator" is needed

The driving fields (Phi_00, Psi_0) are per-multipole LINEAR transfers of the
single scalar potential: Poisson factor A(a)/k^2 per leg
(A = -(3/2) Omega_m H0^2 / a, c = 1), then the screen-Hessian multipliers
L^2/chi^2 (spin-0 trace, Phi_00) and sqrt(L^2(L^2-2))/chi^2 (spin-2 eth^2,
Psi_0), plus the per-leg Sachs weight (1+z)^4. Three-point physics enters the
pipeline through EXACTLY ONE function: the equal-time matter bispectrum
B_delta(k1, k2, k3; z), evaluated at k_i = ell_i/chi(z) per shell (Limber
branch, Eq. `appendix limber reduced bispectrum`) and on the FFTlog k-grid
(exact Wigner-3j branch, ell <= 60). Any matter-bispectrum source therefore
converts exactly; the transfer chain is untouched.

The swap replaces

    D(z)^4 * B_tree(k1, k2, k3)      [tree SPT: 2 F2 P P + 2 cyclic]

by

    B_model(k1, k2, k3; z)           [evaluated per shell at z(lambda)]

## Recommended B_model (given the 2026-08 landscape)

    B_model = B_tree(k; z)                                   k below emulator k_min (exact there)
            = BiHalofit(k; z) * baccoemu_boost(k; z)         inside [0.017, 17.7] h/Mpc, z <= 3
            = BiHalofit(k; z)  (or B_tree at high z)         outside the baccoemu box

with smooth log-k blending at the seams (no hard clips; hard clips alter the
k_par = 0 slice power the collapse conserves, see
cross_shell_vs_equal_time.md). BiHalofit is calibrated at z = 0-3; beyond,
the nonlinear-to-tree ratio shrinks toward 1, so the fallback error dies off
exactly where the fallback is used (coverage_ell_z_mapping.md).

**Squeezed corner (the FK collapsed vertex specifically).** The collapsed
configuration is one soft leg (selected by gamma) plus a coincident hard
pair summed to ell_max: a squeezed bispectrum integrated over the hard
scale. A 2026-08-25 topology probe of the production FK system confirms the
Order-2 <kappa kappa> fold queries the K table exactly on the
(1, cos gamma, cos gamma) family, cos gamma in [0.116, 1] (the table's
cosine_triples span [-1, 1] including the (1,1,1) corner), so the rebuild
must get this whole near-diagonal slice right. Two consequences:
- Among fitting formulas, BiHalofit is actually the BEST-performing in
  squeezed configurations (<= 10% vs N-body, where the older SC01/GM12
  formulas err by > 200%; Heydenreich et al. 2023, arXiv:2208.11686). The
  separate-universe RESPONSE-FUNCTION squeezed bispectrum,
  B_sq(k_s, k_h; z) = R_1(k_h, z) P_lin(k_s, z) P_NL(k_h, z), with
  R_1 = 1 + G_1 - (1/3) dln P_NL/dln k (plus a tidal (mu^2 - 1/3) R_K term
  that vanishes on angle average), G_1 tree value 26/21, is still preferred
  there, for exactness-by-construction in the squeezed limit and for
  propagating baryon response, not because BiHalofit fails. Primary refs:
  Barreira & Schmidt 2017 (JCAP 06, 053, arXiv:1703.09212); G_1 measured
  from separate-universe sims to k ~ 2 h/Mpc at z = 0-2 (Wagner et al.,
  arXiv:1409.6294, 1503.03487). This is the same model the i3PCF pipeline
  (arXiv:2304.01187, Eqs. 2.16-2.19) uses for its squeezed branch, stitched
  to Gil-Marin elsewhere. The soft leg always enters through P_lin, so
  emulator coverage at k_soft is never required.
- The hard legs are boosted by P_NL/P_L(l/chi, z) (up to x3-6 at the
  z = 0.3-0.5 shells for l = 1000; table in coverage_ell_z_mapping.md), so
  the nonlinear swap moves the FK amplitude at ALL angular scales, wide
  angles included, and the ell_max convergence of the hard-pair sum must be
  re-verified with the boosted integrand.

## Two concrete integration points

1. **Table-level (zero downstream churn).** Produce the same NPZ schema the
   `equal_time_limber` callable reads: six arrays
   `zeta_TTT/TTP/TPP/PPP/Bmod/Dmod`, each (1671, 16), on
   `cosine_triples (1671, 3)` x `lambda_shells_Mpc (16,)` (z = 0.1..5.7, the
   shared L2 grid), physical Mpc, lambda-density, with meta flags
   `already_R_contracted = False`, `equal_time = True`,
   `radial_density_measure = "lambda"`, `has_modulus_channels = True`.
   Deployment mechanism: TABLE_PATH is a hardcoded filename inside the
   callable (no env-var/config knob, by sachs_sft design), so either clone
   the callable folder with the new npz (and point the sft-wick YAML's
   coupling_module at the clone) or replace the npz in place under the
   canonical name; note the callable caches the table at module level.
   Either way sft-wick and the FK 2PCF drivers are otherwise untouched.

2. **Build-level (canoes API).** Replace the `B_delta = 2 F2 P P + cyc`
   assembly inside `compute_kappa3_zeta_table` / `compute_kappa3_sigma3_high`
   (and the `_mod_` variants for Bmod/Dmod) with a callable
   `B_delta(k1, k2, k3, z)`. CAVEAT: the tree-level growth scaling
   `B(z) = D(z)^4 B(0)` is hardwired (per-leg D promotion); a z-dependent
   nonlinear B must bypass that scaling and be evaluated per shell.
   Branch windows (kappa3.py:3552-3557): the HIGH (Limber) branch takes
   triangles whose LARGEST leg satisfies 60 < ell_max_leg <= 1000, but the
   individual legs run down to ell ~ 0 (Gauss nodes reach ell ~ 1e-3), so
   B_model must also be supplied at soft k INSIDE the HIGH branch (the
   P_lin soft leg of the squeezed model covers this). The LOW branch (all
   three ells <= 60, exact Wigner-3j) can keep tree-level except possibly
   at the one or two lowest shells (z <= 0.25, where ell = 60 reaches
   k ~ 0.06-0.21 h/Mpc, an 18% one-leg boost at z = 0.1); the ell_cut = 60
   seam boundary probe decides, and the b_delta_fn API extends into the
   LOW FFTlog path if needed.

## Bookkeeping that must NOT change (it lives downstream of B)

- Poisson factor A(a)^3 applied per leg internally
  (`spt_kind = "tree_phi"`); feed the emulator's B in MATTER form.
- h-convention: canoes computes in h-units from CAMB-format P(k) inputs
  (k in h/Mpc, P in (Mpc/h)^3) and converts with h^4 (native per-shell
  density A^3 P P / chi^4 ~ (h/Mpc)^4; the h^6 variant was a bug, fixed
  2026-06-09). baccoemu/BiHalofit work in h-units natively, so the natural
  place to plug them is BEFORE the h^4 conversion. Verify with the
  h-invariance boundary test (physical-vs-h fold must be 1.000).
- Radial measure: the two equal-shell collapse deltas carry
  (d lambda/d chi)^2 = (1+z)^-4 (fix of 2026-06-09/10); the per-leg Sachs
  weight (1+z)^4 and the 1/chi^4 Limber geometry are separate factors.
- Spin-2 channel algebra (TTT/TTP/TPP/PPP + modulus B, D channels and the
  parity zeros) is independent of the B model.

## Cosmology matching

Pipeline fiducial: Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97
(corr_op meta). baccoemu takes cold-matter conventions (omega_cold,
sigma8_cold); with massless neutrinos (as here) cold = total and the mapping
is direct. sigma_8 must be computed from the SAME z = 0 linear P(k) that
feeds the build being upgraded. Caution: the deployed 3-point table was
built from `PCAMBz0.txt`, NOT the 2-point side's
`PCAMB_pyccl_stf_fid_z0.txt` (a README-documented difference pinned by the
FK baseline); the latter integrates to sigma8 = 0.810 exactly (top-hat
integration, matching its header). Pick which file anchors the closure test
explicitly, and compute sigma8 for `PCAMBz0.txt` before evaluating the
boost (roadmap Phase 0). Verbatim baccoemu call (readthedocs-verified):
`k, boost, flags = baccoemu.Matter_bispectrum().get_baryonic_boost(
omega_cold=..., omega_baryon=..., sigma8_cold=..., expfactor=...,
M_c=..., eta=..., beta=..., M1_z0_cen=..., theta_inn=..., k1=..., k2=...,
k3=...)` with expfactor in [0.25, 1.0], i.e. z <= 3 (not 3.2).

## Validation recipe after the swap

1. Closure test: independent Limber projection of the SAME B_model compared
   against the pipeline zeta (the fastnc input-validation pattern).
2. Boundary tests at every blend seam and at the emulator domain edges,
   probing seam +/- delta including extreme shells (z ~ 0.1 and z = 5.7),
   per the boundary-validation methodology.
3. h-invariance fold and the existing FK kk baseline pin (re-derive the
   number end-to-end; the on-disk meta stamp is known to carry a stale
   template value).

## Environment status (checked 2026-08-25)

- canoes moved to `/Users/zzhang/projects/angular_statistics/canoes` (src
  layout, git HEAD d7b20ce, contains the 1ddc141/e19b3e4/e9a34aa fixes).
  The callable's `_CANOES_ROOT`, the build script, and the editable pip
  installs (`_editable_impl_canoes.pth`) in BOTH the sft-wick and PyCCL
  conda envs still point at the dead `/Users/zzhang/projects/canoes`, so
  `import canoes` is broken in both envs. Working import: canoes' own
  `.venv` with `PYTHONPATH=.../angular_statistics/canoes/src`. Fix before
  any rebuild (roadmap Phase 0).
- baccoemu is installed in NO environment yet.

## Expected physical impact

The swap raises the squeezed vertex at ALL angular scales, because even the
wide-angle vertex carries the coincident hard pair (see the squeezed-corner
paragraph above and coverage_ell_z_mapping.md): at arcminute scales an
O(1)-factor increase (BiHalofit vs tree at k ~ 0.3-2 h/Mpc), at wide angles
a kernel-weighted mix of the per-shell hard-pair boosts, plausibly tens of
percent (to be quantified by the rebuild probe), plus a few-percent baryonic
suppression from the boost at the highest k. The paper's intro statement
that the tree-level FK amplitude is conservative covers both regimes. Only
the fully-soft triple population (all three ell <= 60 at the distant
shells) is strictly unchanged.
