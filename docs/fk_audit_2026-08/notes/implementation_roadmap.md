# Roadmap: Nonlinear-Bispectrum FK Rebuild (tree SPT -> BiHalofit x bacco boost)

> **Correction, 2026-09-02.** The absolute FK amplitudes in this note are
> `ell_max = 1000` values obtained with the sorting callable and the production
> linear-in-lambda fold. The manuscript's FK is evaluated at `ell_max = 15360`
> with the permutation-aware callable, where the same quantity is `+1.9480e-5`
> at `gamma = 0.5'` (2.31% of Order-0), a factor 4.5 larger. Ratios between
> variants at a fixed cutoff are unaffected. See
> [`../FK_BASELINE_NUMBERS.md`](../FK_BASELINE_NUMBERS.md) for the full key and
> `sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/` for the
> current product.

Goal: recompute the FK term (three-point cumulant leakage into the lensing 2PCF) with a nonlinear matter bispectrum, replacing the tree-level SPT `2 F2 P P + cyc` inside the `equal_time_limber` kappa3 vertex. Options in increasing sophistication: (A) BiHalofit gravity-only; (B) BiHalofit x baccoemu baryonic boost; (C) response-function squeezed model `R_1(k_h) P_lin(k_s) P_NL(k_h)` as a targeted cross-check of the collapsed FK corner.

Design decision up front: build-level integration (patch canoes' B assembly) as the compute path, table-level as the deployment surface. The patched build emits the SAME NPZ schema the deployed callable already reads, so sft-wick, the callable, and the FK 2PCF drivers see only a new `TABLE_PATH`. Pure table-level post-processing is impossible anyway: the nonlinear B is not separable as `D(z)^4 B(0)`, and the table stores zeta (already ell-integrated), not B.

---

## Phase 0 - Prerequisites (half a day)

1. **Fix the stale canoes path.** canoes moved: it now lives at `/Users/zzhang/projects/angular_statistics/canoes` (src layout). The build script hardwires the dead path `/Users/zzhang/projects/canoes` at `build_equal_time_limber_table.py:47`, and `_local_cosmo_pk.py:30` points `_DEFAULT_CAMB_TABLE` at `/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt` (file exists at the new location: `/Users/zzhang/projects/angular_statistics/canoes/examples/data/PCAMBz0.txt`). Verified working import: `PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src` with `/Users/zzhang/projects/angular_statistics/canoes/.venv/bin/python` (no conda env has canoes; the .venv needs the PYTHONPATH or a `pip install -e`). Action: update both hardwired paths (constant at top of file, per no-hardcoded-values rule use one `_CANOES_ROOT`), smoke-test with `--smoke` before anything else.
2. **Environments.**
   - Vertex build: canoes `.venv` (jax lives there).
   - Standalone B_model module + tests: `PyCCL` env (`/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python`) — needs pyccl for `P_NL` (halofit) and linear P(k).
   - **baccoemu: installed NOWHERE currently** (checked PyCCL, sft-wick, sft-sachs, canoes .venv). Install into the **PyCCL env** first (that is where the standalone module and its tests run). For the vertex rebuild, ALSO install into the canoes `.venv` if TensorFlow coexists cleanly with canoes' jax; if it conflicts, fall back to pre-tabulating the (smooth, <2%-accurate) boost onto the exact per-shell (k1,k2,k3) query set in the PyCCL env and loading it as a plain interpolator inside the build (boost smoothness makes trilinear interpolation error negligible against BiHalofit's own 10-15%).
3. **Pin fiducial sigma8.** Compute sigma8 from `PCAMBz0.txt` (the P(k) the deployed table was built with; documented difference vs the 2-point side's `PCAMB_pyccl_stf_fid_z0.txt`, which integrates to sigma8 = 0.810 exactly). Record both numbers in a small `fiducial_sigma8.md` note; the baccoemu `sigma8_cold` input MUST be the PCAMBz0 value (the boost must be evaluated at the cosmology of the P(k) actually fed to the build). Fiducial: Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97. NOTE: canoes `FiducialCosmology` does not pin Omega_b — baccoemu needs `omega_baryon`; choose (e.g. Planck-like 0.049), record it, and quantify boost sensitivity to the choice.
4. **BiHalofit implementation choice.** Research-and-reuse order: (i) check the fastnc install at `/Users/zzhang/projects/fastnc` (already vetted against this pipeline in the 2026-06 fastnc closure work) for a usable BiHalofit module; (ii) else the i3PCF/threepoint public implementations; (iii) else vendor a small pure-numpy port of the Takahashi et al. 2020 reference code into `driver_field_emulators/code/`, pinned against reference-code output values in tests. Whatever the source, pin 5-10 published/reference (k1,k2,k3,z) values as regression numbers.

Exit criteria: `--smoke` build runs green from the fixed paths; baccoemu imports and reproduces its README example; sigma8 note written; BiHalofit source chosen with pin values in hand.

## Phase 1 - Standalone `B_model` module (1-2 days)

Location: `/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/b_model/` (small focused files, not one monolith).

```
b_model/
  __init__.py        # public API
  tree.py            # B_tree(k1,k2,k3, pk_lin) = 2 F2 P P + 2 cyc  (z=0 shape)
  bihalofit.py       # B_bihalofit(k1,k2,k3, z, cosmo)  [chosen impl or vendored port]
  bacco_boost.py     # boost(k1,k2,k3, z) -> R = B_hydro/B_grav  (lazy-load emulator)
  response.py        # B_sq(k_s, k_h, z) = R_1(k_h,z) P_lin(k_s,z) P_NL(k_h,z); R_1 = 1 + G_1 - (1/3) dlnP_NL/dlnk
  blend.py           # C^1 smoothstep windows in ln k / z; NO hard clips
  core.py            # B_model(k1,k2,k3, z, model={"tree","A","B","C"}) dispatcher
  tests/
```

**API contract (document in module docstring, enforce with guards):**
- All k in **h/Mpc**; P in (Mpc/h)^3; B returned in **(Mpc/h)^6**; **MATTER** form (no Poisson factor A(a)/k^2 — canoes' `tree_phi` applies it per leg internally); equal-time, one scalar z per call; vectorized over (N,3) triangle arrays.
- The h^4 native-density conversion, the (1+z)^-4 lambda-measure Jacobian, the per-leg (1+z)^4 Sachs weight, and the 1/chi^4 Limber geometry all live DOWNSTREAM in canoes and are not touched here (conversion_interface.md "Bookkeeping that must NOT change").
- Blending scheme (per coverage_ell_z_mapping.md): `B_tree` below the low-k seam (exact there); `BiHalofit [x bacco boost]` inside k in [0.017, 17.7] h/Mpc and z <= ~3.2; `BiHalofit` (or tree at high z, where NL/tree -> 1) outside the bacco box; smooth log-k smoothstep over ~half a decade at each seam. Rationale: hard clips alter the k_par = 0 slice power that the equal-shell collapse conserves.
- Domain guards: fail loud on non-finite output, non-positive P, triangle-inequality violation; never silently extrapolate the emulator.

**Tests (pytest, boundary-validation methodology at EVERY seam):**
- Seam agreement: evaluate both sides directly (bypassing the dispatcher) at each seam +/- delta; rel_err gates 1e-3 / 1e-1 per the methodology. Seams: tree/BiHalofit low-k, bacco k_min and k_max edges, z = 3.2 boost ceiling, high-z tree rolloff.
- Extreme corners always included: z corresponding to the extreme shells (z ~ 0.1 and z = 5.7), k endpoints, maximally squeezed triangles (k_s/k_h ~ 1e-3), equilateral, folded.
- Squeezed-limit cross-check: model A/B vs model C on `(k_s, k_h, k_h)` configurations at the shells where the boost table predicts x1.4-x13 (z = 0.3-2, l = 250-2000).
- Pins: BiHalofit reference values; boost -> 1 at bacco domain edges; `B_model("tree") == B_tree` exactly.
- h-bookkeeping self-consistency: evaluating the module through two internal h foldings must agree to machine precision (module-level h-invariance, separate from the Phase 3 pipeline fold).

Exit criteria: all tests green in the PyCCL env; a one-page plot script showing B_model/B_tree vs k for a few z (visual sanity, goes in `driver_field_emulators/code/`, not `figures/`).

## Phase 2 - Vertex rebuild (2-3 days incl. one overnight build)

**Entry point: the HIGH (Limber) branch of the canoes build, NOT the LOW 3j branch.**

- Patch sites (canoes `src/canoes/sachs/kappa3.py`): `compute_kappa3_sigma3_high` assembles `b_delta_today = 2*(f12 p1 p2 + f13 p1 p3 + f23 p2 p3)` at line 3618 and scales by `growth**4` at 3622 (per-shell `growth = M[0,i_lam]/(1+z)`); `compute_kappa3_mod_sigma3_high` repeats the pattern at 4032-4035. Add an optional `b_delta_fn(k1, k2, k3, z) -> B` parameter: when given, evaluate per shell at `z(i_lam)` with `k_i = ell_i / chi(z)` and **bypass the hardwired `growth**4` promotion entirely** (that scaling IS the tree-level `D^4 B(0)` assumption). When `b_delta_fn=None`, behavior is bit-identical to today (regression guarantee). Immutable-style: new keyword, no change to existing call sites; canoes work goes on a feature branch (`feat/nonlinear-bdelta-vertex`), committed in canoes' own repo with its test suite.
- **LOW branch (exact Wigner-3j, ell <= 60) stays tree-level.** Justification: at ell <= 60 all shells except the lowest map to k <= 2.6e-2 h/Mpc (coverage table; the z = 0.1 shell, chi = 436 Mpc, reaches k = 0.21 h/Mpc with an 18% one-leg boost at ~0.5x peak kernel weight), so `B_model = B_tree` by construction of the blend almost everywhere. Verify, do not assume: boundary probe at the ell_cut = 60 dispatch seam — evaluate the tree LOW and the boosted HIGH both AT ell ~ 60 (bypassing the LOW/HIGH combine) across all 16 shells including both extremes; the boost at ell = 60 must be within the seam tolerance (predicted ~1.00 at every shell; if the z ~ 0.1 shell fails, lower nothing — extend `b_delta_fn` into the LOW FFTlog path instead, which the API already permits since LOW consumes B on the FFTlog k-grid).
- **Build script:** extend `build_equal_time_limber_table.py` with `--b-model {tree,A,B}` (and `--b-model-tabulated <npz>` for the pre-tabulated-boost fallback). Reuse the existing LOW products where possible: LOW is unchanged, so cache LOW/`mod` LOW once and rebuild only HIGH + mod-HIGH per B-model variant (`kappa3_combine_low_high` takes them as-is). New meta keys: `b_delta_model`, `bihalofit_source`, `bacco_version`, `sigma8_pin`, `omega_baryon`, plus a fresh `fk_kk_baseline` note (the current on-disk stamp carries the tree value; the README warns the meta stamp history is error-prone — re-derive, never copy).
- **ell_max convergence sweep.** Treat convergence at `ell_high_max = 1000` as UNPROVEN for tree as well as boosted: the P^1 hard-leg toy diagnostic (coverage_ell_z_mapping.md, consequence 2) gives S(4000)/S(1000) = 1.44-1.87 (tree) and 2.2-3.25 (boosted) at z = 0.5-2; the real Gaunt sum has band sign cancellations the toy lacks, so only the sweep decides. Sweep `ell_high_max` in {1000, 1500, 2000, 3000, 4000} on the squeezed family `(1, cos g, cos g)` at ~6 gamma x all 16 shells, for BOTH tree and boosted (direct HIGH calls, no full-table rebuild); pick the value where zeta changes < 1% per shell and freeze it. Also re-check `n_ell = 96` log-spacing is adequate for the steeper integrand.
- **Wall-clock.** The full production kappa3 build is the ~hours-class item (historically ~14 h for full builds; the current script's LOW Nmax=128 over 1671 triples x 16 shells dominates). Since LOW is reused, expect the nonlinear rebuild to cost roughly the HIGH + mod-HIGH share plus B_model evaluation. GOOD NEWS from the source extraction (input_ready_b_to_zeta_spec.md sec 2): B is evaluated on the (96, 96, 64) = 589,824-node quadrature grid ONCE PER SHELL, shared across ALL cosine triples and all six channels — total ~9.4e6 B_model calls per full table, NOT per-triple. Vectorize `b_delta_fn` over the (96,96,64) batch; if baccoemu NN eval is the bottleneck, switch to the pre-tabulated boost. Sequence: `--smoke` first, then one overnight full build per B-model variant. Watch RAM: keep `n_workers` low per the loky OOM lesson.
- **Grid densification (REQUIRED, new finding).** The production FK fold queries only the (1, c, c) family, and on covgrid16 (cosine spacing 1/7) the cKDTree IDW lookup puts >= 96% weight on the single (1,1,1) corner cell for 27 of the 40 production gammas (0.5'..~230'), with 10-45% off-family contamination beyond — the CURRENT table cannot resolve gamma dependence at small/mid angles. The rebuilt table must add (1,c,c) rows at the 40 production c values (or log-spaced in (1-c) over [1e-8, 0.9]). Build BOTH a tree table and the nonlinear tables on the densified grid: the densified-tree-vs-covgrid16-tree comparison isolates the grid effect on the existing paper FK curves from the new physics. Also extend the Bmod/Dmod combine to use the guarded `kappa3_combine_low_high` checks, and decide explicitly about the lambda floor (queries below 396.6 Mpc are currently zeroed).
- **Output naming (archive, don't delete):** new NPZ e.g. `equal_time_limber_kappa3_z5covgrid16_nlB_bihalofit.npz` / `..._bihalofitXbacco.npz` alongside the tree table; the tree table stays in place untouched (it IS the current paper baseline). Nothing is overwritten; superseded intermediates go to a dated `_archive_<YYYY-MM-DD>/` folder.

Exit criteria: `b_delta_fn=None` path reproduces the existing table bit-for-bit (or to FP noise); ell_cut seam probe green; ell_max frozen with the sweep documented; full A and B tables built.

## Phase 3 - Validation (2 days)

Run in order; each is a gate.

1. **Tree-path regression through the NEW code.** Rebuild the tree table via the patched script (`--b-model tree`) on the ORIGINAL covgrid16 grid and re-run the FK 2PCF on the production geometry (40-pt grid, t_final = 2313.03, equal_time): must reproduce **FK kk(0.5') = +4.29e-6** to numerical noise (the deployed post-Jacobian baseline; the '+1.219e-4' in the BUILD docstring and the '+3.086e-5' in the on-disk meta are BOTH stale stamps — see input_ready_b_to_zeta_spec.md sec 5; provenance check confirmed the deployed NPZ values already carry the (1+z)^-4 Jacobian, so a rebuild from current source reproduces them directly and must NOT re-apply it). This isolates "new plumbing" from "new physics" before any nonlinear number is trusted. THEN rebuild tree on the densified grid and quantify the grid effect on FK xi(gamma) as its own result.
2. **Closure test.** Independent Limber projection of the SAME `B_model` (standalone script in `driver_field_emulators/code/`, the fastnc input-validation pattern) vs the pipeline zeta on the squeezed family, per shell. Order-of-magnitude/shape agreement gate as in the 0.70x fastnc closure. This catches unit/measure errors in the new `b_delta_fn` fold independently of canoes internals.
3. **h-invariance fold.** Physical-vs-h build fold must be 1.000000 (the decisive test that caught the h^6 bug). B_model works natively in h-units and plugs in BEFORE the h^4 conversion; this test is the contract check.
4. **FK kk baseline re-derivation.** End-to-end sft-wick FK with the nonlinear table -> the NEW baseline number, derived (not stamped); document tree -> A -> B shifts at gamma = 0.5' and at the wide-angle crossover (~2 deg).
5. **ell-band decomposition vs the archived table.** Re-run `ell_band_decomp.py` with the boosted HIGH; compare band-by-band against the archived z = 5.7-shell table in `ELL_BAND_DECOMPOSITION.md` (tree: high-ell share 100% at 0.5', crossover ~42'). Expected: high-ell bands boosted, low-ell [2,30] band unchanged, crossover gamma moves outward; any change in the LOW bands is a bug.
6. **Extreme-shell boundary tests.** Direct probes at the z ~ 0.1 shell (k = ell/chi reaches 10-100 h/Mpc, entirely in the blend/fallback region) and the z = 5.7 shell (fallback-to-tree region): finite, no catastrophic blowup, seam continuity — the canoes ell=5000 lesson says the fix must be re-probed at the extremes before landing.
7. **Promote probes to pytest**: parametrized boundary tests next to `test_equal_time_limber_helicity.py` (seam cells, extreme shells, pinned tree-regression values, pinned new-baseline value).

## Phase 4 - FK 2PCF + paper-figure regeneration (1-2 days)

- **sft-wick runs** (sft-wick env; chdir to the YAML's parent before `run_workflow`; callable-module paths YAML-relative, output/cache CWD-relative; `n_jobs <= 8` for fat NPZs): re-run the order-2 FK configs with `TABLE_PATH` pointed at the nonlinear NPZ — the production FK 2PCF (equal_time geometry, n_gauss 24/32) for kk, xi_+, xi_-, kappa-gamma. O0 and FF runs are UNTOUCHED (FK-only input change).
- **Affected figures** (check `SFT-lensing-paper-analyses/figure_manifest.md` for the authoritative map):
  - `figures/analysis3_NLO_FFFK.pdf` (Fig 11 family, O0/FF/FK decomposition) — regenerate.
  - `figures/analysis3_cl_full_vs_O0.pdf` (Fig 12, curved-sky Wigner-d C_ell) — regenerate; the FK low-ell excess claim must be re-checked (the boost mostly lifts small angles, but the wide-angle vertex carries the boosted hard pair too).
  - `figures/zeta_driving_field_slices_draft.pdf` — regenerate (zeta changed).
  - `figures/multiz_kappa_xi_cl.pdf` (2x5 multi-z) — the FK slice is PCHIP-in-z from the cached 20-z grid; either re-derive the cache with the new table or mark the figure tree-level explicitly. Decide with the user.
  - `figures/cl_EB_polarization.pdf` — reuses existing run outputs; re-run only if FK feeds its runs.
  - Appendix fastnc validation figure — its input model changed; regenerate or re-scope the caption (it validates open-triangle machinery, not the collapsed slice).
- Keep the tree-level outputs archived alongside (dated folder), since the paper may present tree vs nonlinear as a comparison rather than a replacement.
- **No `sections/*.tex` edits during the analysis runs.** Text updates (intro "tree-level is conservative" statements, discussion/forecast FK numbers, figure captions in timeless prose) are a separate pass after numbers stabilize.

## Risks and open questions

1. **canoes path drift** (confirmed): the hardwired `/Users/zzhang/projects/canoes` is dead; also verify the `_ARCHIVED_L2` fallback NPZ still exists after the 2026-06-17 cleanup before relying on it (prefer `--l2-npz` at the R_contracted table).
2. **cos-grid coverage at the collapsed family** `(1, cos g, cos g)`: the production FK fold queries this whole family with cos g in [0.116, 1] (topology probe 2026-08-25; the older cos = (1,1,1) colocation reading was a rounding artifact of a single 0.5' probe point). The deployed covgrid16 spans [-1, 1] including the (1,1,1) corner, but re-assess cKDTree nearest-neighbor interpolation error along the near-diagonal slice: the nonlinear zeta is steeper in the squeezed corner, so grid-resolution error grows exactly where FK reads. May need `n_cosine` densification near cos = 1.
3. **sigma8 / P(k) consistency with the 2-point side**: the FK baseline pins `PCAMBz0.txt`, the 2-point side uses `PCAMB_pyccl_stf_fid_z0.txt`. The boost and BiHalofit must be evaluated at the sigma8 of the P(k) actually used; quantify the residual inconsistency between the two tables and decide whether this rebuild is the moment to unify them (that would ALSO shift the tree baseline — keep as an explicit, separately-gated decision).
4. **baccoemu conventions**: cold-vs-total (massless-nu fiducial => cold = total, direct mapping, but pass `sigma8_cold` and `omega_cold` correctly); Omega_b is not pinned by canoes' FiducialCosmology; baryonification parameters for option B need a fiducial choice (and are a systematic axis, not a single number). Also TF-vs-jax coexistence in the canoes .venv (fallback: pre-tabulated boost).
5. **Response-model G_1 data availability** (option C): separate-universe G_1 measurements/fits (Barreira-Schmidt formalism) cover limited (k, z); confirm a usable public fitting function before promising option C as production. Recommended stance: C is the squeezed-corner cross-check of A/B, not the deployed model.
6. **BiHalofit domain**: calibrated z = 0-3, ~10-15% at k < 3-10 h/Mpc. In squeezed configurations it is the best-performing fitting formula (<= 10%, vs > 200% for SC01/GM12; Heydenreich et al. 2023, arXiv:2208.11686), but its ~10% model error may still dominate the FK shift's own error budget; the paper must state this (no unearned validation claims), and the response-function model (option C) is the squeezed-corner cross-check.
7. **ell_max extension cost**: if the sweep pushes `ell_high_max` to 2000-3000, HIGH runtime and the flat-sky approximation budget both grow; also confirm sft-wick's response-R application makes no assumption about the table's ell content (no double counting).
8. **Cross-shell sigma3 measure** remains OPEN from 2026-06-09 (does not feed the paper figures, but do not let the rebuild silently inherit it if that path is ever activated).
