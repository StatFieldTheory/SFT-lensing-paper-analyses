# Repairing the non-PSD two-ray `Sigma2`

*Synthesis of the four assessments, 2026-08-28. Every number below marked
**[V]** was re-measured by me from scratch in this session; the scripts are
`psd_fix/v1_verify.py` ... `psd_fix/v8_order0_patch.py`. Nothing outside
`psd_fix/` was modified; `git status --porcelain` over the repo is empty.*

---

## Verdict

**Adopt R1 (a spline broken at the corr_op table nodes), as a patch to the
`Sigma2Builder` *class* in `driver_stats.py`, in the same commit as a
quadrature fix to `order0_mc`.** The repair removes a covariance that is
returned with the *wrong sign* at two of the 1000 production nodes. It changes
nothing the paper reports: the exact FK expectation moves by -6.5e-6 relative
and the extrapolated marker ratios by at most +0.00013, against a quoted
headline uncertainty of +/- 0.006.

Separately, and at higher *scientific* priority but on a longer timescale, the
equal-time diagonal that corr_op serves is itself wrong by -5% to -24%
mid-cell. That is a bigger error than anything R1 fixes. The right lever for
it is a one-dimensional, ridge-aware reading of the table, **not** a denser
table: halving the cell width reduces the mid-cell error by only 1.50x **[V]**.

---

## 1. Root cause

`corr_op` serves its 21x21 lambda table through a **bilinear**
`RegularGridInterpolator` (`canoes/sachs/sft_input/corr_op/table.py:962`;
`interp_method` defaults to `"linear"` and `corr_op_C_callable.py:91` does not
override it), so on the equal-time diagonal `C_ab(t,t)` is *exactly* piecewise
quadratic (degree-2 fit residual 5.1e-14 to 1.8e-13 in every cell tested) and
only `C^0`: across the node at lambda = 1300 the one-sided derivative runs
+5.984e-13 on the left and -2.205e-13 on the right, a **sign flip**, while the
value itself is continuous to 4.1e-12 relative **[V]**.
`Sigma2Builder._density_splines` fits a single globally `C^2` `CubicSpline`
through 160 uniform samples (12.09 Mpc) of that function and returns its
analytic derivative, so it is forced to interpolate a function whose first
derivative jumps at 19 interior nodes, and it rings with alternating sign in
the cells on either side. The ringing is large enough that at lambda = 1308.93
the builder returns `Sigma2_00 = -1.895e-14` where the true value is
`+1.112e-13` **[V]**: the 6x6 becomes negative definite in all six
eigenvalues, `_cholesky_psd`'s eigen-floor returns `L = 0`, and the simulation
injects no driving noise at that node at all.

| quantity | value |
|---|---|
| interior table nodes in the builder's span [406, 2328] | 19 |
| non-PSD nodes, stock builder, 1500-point scan | 4 (0.5'), 4 (1'), 4 (17.3'), 3 (114.3') **[V]** |
| worst min-eig / max-eig, stock | -4.00e-02, -3.46e-02, -4.38e+00, -1.09e+01 **[V]** |
| same scan, repaired | 0 non-PSD at every gamma; worst ratio +3.13e-04 **[V]** |
| production failure window | lambda in [1307.7, 1311.1], 3.8 Mpc = 0.2% of the ray |

### Three corrections to the standing diagnosis

1. **The density is discontinuous, not kinked.** `DIAGNOSIS.md` says `C` is
   `C^1` and the density has a slope discontinuity. It is `C^0`; `dC/dlam`
   jumps (and flips sign at 1300), so `Sigma2` itself jumps at every interior
   node **[V]**. This strengthens the case for R1 (no globally smooth density
   can ever be right, which is why `R3Smoothed` was doomed a priori) and it is
   the reason a jump-aware quadrature is needed downstream.
2. **The max-error figures should be withdrawn.** `R0Reference` is not a valid
   yardstick within about 1 Mpc of a node: its dense per-cell finite
   differences smooth the jump. Measured against an analytic evaluator the
   stock builder's max error is ~960%, not 42.5%, and R1's is ~1.05%, not
   0.20%. Both *medians* (0.61% and ~1e-5) reproduce exactly, so every
   median-based conclusion in `DIAGNOSIS.md` stands.
3. **There are 19 ringing windows (17 in [410, 2310]), not eleven**, one
   starting at each interior node. Only one of them, [1307.5, 1311.2], drives
   the 6x6 indefinite, which is consistent with the reported 2/600.

---

## 2. The repair, and the exact minimal patch

**R1Knots.** Same interpolate-then-differentiate recipe as production, with the
spline broken at the corr_op table nodes so no spline cell ever straddles a
jump. It is not knot tuning: because `C(t,t)` is exactly piecewise quadratic,
the density is analytic with **zero free parameters**, and R1 reproduces it to
a median 8.4e-6 (reported by the adversarial pass; my transcription of the
patch is bitwise identical to `R1Knots`, max relative difference 0.00e+00 over
601 lambda at two gammas **[V]**).

Two things this patch must *not* be justified by, in a commit message or a code
comment:

* **"the negative eigenvalue disappeared."** Five of seven deliberately *wrong*
  break placements also give 0 non-PSD with an identical min-eig/max-eig.
  Breaking a global spline anywhere shortens its support and damps the ringing.
  PSD restoration is not diagnostic. Justify it by agreement with the analytic
  density and by the wrong-sign covariance it removes.
* **"the pointwise calibration failed because of this defect."**
  `PAPER_PROPOSAL.md:244-248` says the pointwise-versus-smeared wording is
  non-specific "because the reason is a defect in the tabulated covariance (it
  is not positive semi-definite at isolated lambda; see NOTES.md section 3.2)".
  Repairing the defect does not rescue the pointwise calibration: the
  nominal-calibration exact expectation goes 3.91x the fold -> 3.11x at
  gamma = 1' and 1.85x -> 2.48x at 17.3', i.e. worse at the second. The
  pointwise calibration fails for the physical reason the paper text already
  gives. The paper prose is unaffected; only that rationale note is wrong.

### Patch (do not apply from here; shown for review)

Three hunks in `Sigma2Builder` plus one in `order0_mc`. Nothing else in
`driver_stats.py` references `_nodes`, `_D4_nodes` or `_density_splines`, and
nothing outside the file does **[V]**.

```diff
--- a/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/driver_stats.py
+++ b/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt/driver_stats.py
@@ -124,6 +124,17 @@ ANCHOR_C0: Final[float] = 1.0 / ANCHOR_RATIO
 
+# corr_op interpolates its table BILINEARLY in (lam1, lam2), so on the equal-time
+# diagonal C_ab(t,t) is exactly piecewise quadratic and only C^0: dC/dlam jumps at
+# every table node (at lam = 1300 it flips sign, +5.98e-13 -> -2.21e-13), and the
+# density Sigma2 = dC/dlam + 4 (D'/D) C is genuinely DISCONTINUOUS there.  A single
+# globally C^2 spline through that function rings: it returned Sigma2_00 with the
+# wrong sign (-1.90e-14 against a true +1.11e-13) at lam = 1308.9, which made the
+# 6x6 negative definite and silenced the driving noise at that node.  Sample and
+# spline INSIDE each table cell instead, so no spline cell straddles a jump.
+_TABLE_LAM: Final[NDArray[np.float64]] = np.load(
+    _corr_op.TABLE_PATH, allow_pickle=True)["lambda_grid"].astype(float)
+
 # Reference rays for the two-ray reduced problem.
 _N1: Final[NDArray[np.float64]] = np.array([0.0, 0.0, 1.0])
 
@@ -162,13 +173,22 @@ class Sigma2Builder:
     def __init__(self, background: _bg.Background | None = None, n_nodes: int = 160,
                  lam_lo: float = 406.0, lam_hi: float = 2328.0,
                  apply_c0: bool = True) -> None:
         self.bg = background if background is not None else _bg.Background()
         self.lam_lo = float(lam_lo)
         self.lam_hi = float(lam_hi)
         self.apply_c0 = bool(apply_c0)
-        self._nodes = np.linspace(self.lam_lo, self.lam_hi, int(n_nodes))
-        self._D4_nodes = np.asarray(self.bg.D(self._nodes), dtype=float) ** 4
-        self._cache: dict[float, list[CubicSpline]] = {}
+        # Break the sampling at the corr_op table nodes; keep the nominal spacing
+        # that n_nodes implies over the whole span (160 -> 12.09 Mpc -> 210 samples).
+        inner = _TABLE_LAM[(_TABLE_LAM > self.lam_lo + 1e-9)
+                           & (_TABLE_LAM < self.lam_hi - 1e-9)]
+        self._breaks = np.concatenate([[self.lam_lo], inner, [self.lam_hi]])
+        dlam = (self.lam_hi - self.lam_lo) / max(int(n_nodes) - 1, 1)
+        self._cell_nodes = [
+            np.linspace(lo, hi, max(6, int(np.ceil((hi - lo) / dlam)) + 1))
+            for lo, hi in zip(self._breaks[:-1], self._breaks[1:])
+        ]
+        self._cache: dict[float, list[list[CubicSpline]]] = {}
 
-    # -- internal: per-cos spline of d/dlam[D^4 C_ab(t,t)] / D^4 ---------------
-    def _density_splines(self, cos_gamma: float) -> list[CubicSpline]:
+    # -- internal: per-cos, per-cell spline of d/dlam[D^4 C_ab(t,t)] / D^4 -----
+    def _density_splines(self, cos_gamma: float) -> list[list[CubicSpline]]:
         key = round(float(cos_gamma), 12)
         cached = self._cache.get(key)
         if cached is not None:
             return cached
         n2 = _n2_of_cos(cos_gamma)
-        # Sample the corr_op equal-time (3,3) diagonal on the lambda nodes.
-        c_diag = np.empty((self._nodes.size, N_COMP, N_COMP), dtype=float)
-        for i, t in enumerate(self._nodes):
-            c_diag[i] = _corr_op.C_fn(_N1, float(t), n2, float(t))
-        f = self._D4_nodes[:, None, None] * c_diag  # D^4 * C_ab(t,t)
-        splines: list[CubicSpline] = []
-        for a in range(N_COMP):
-            row: list[CubicSpline] = []
-            for b in range(N_COMP):
-                row.append(CubicSpline(self._nodes, f[:, a, b]))
-            splines.append(row)  # type: ignore[arg-type]
-        # Flatten to a 9-list for cache compactness.
-        flat = [splines[a][b] for a in range(N_COMP) for b in range(N_COMP)]
-        self._cache[key] = flat
-        return flat
+        # Sample the corr_op equal-time (3,3) diagonal cell by cell.
+        per_cell: list[list[CubicSpline]] = []
+        for xs in self._cell_nodes:
+            c_diag = np.empty((xs.size, N_COMP, N_COMP), dtype=float)
+            for i, t in enumerate(xs):
+                c_diag[i] = _corr_op.C_fn(_N1, float(t), n2, float(t))
+            f = (np.asarray(self.bg.D(xs), dtype=float) ** 4)[:, None, None] * c_diag
+            # 9-list per cell, for cache compactness.
+            per_cell.append([CubicSpline(xs, f[:, a, b])
+                             for a in range(N_COMP) for b in range(N_COMP)])
+        self._cache[key] = per_cell
+        return per_cell
 
     def matrix(self, cos_gamma: float, lam: float) -> NDArray[np.float64]:
         """Bare equal-time (3,3) covariance density at separation cos and shell lam."""
-        flat = self._density_splines(cos_gamma)
+        per_cell = self._density_splines(cos_gamma)
         lam_c = float(np.clip(lam, self.lam_lo, self.lam_hi))
+        # The density is two-valued at a table node; take the left limit there.
+        j = int(np.clip(np.searchsorted(self._breaks, lam_c) - 1,
+                        0, len(per_cell) - 1))
+        flat = per_cell[j]
         d4 = float(self.bg.D(lam_c)) ** 4
@@ -278,17 +298,26 @@ def order0_mc(cos_gamma: float, lam_f: float | None = None,
     b = builder if builder is not None else _builder()
     lf = float(lam_f) if lam_f is not None else _bg.LAM_SOURCE_BASELINE
+    # Sigma2 is discontinuous at the corr_op table nodes (see _TABLE_LAM), so a
+    # single global Gauss-Legendre rule does not converge: n_gauss = 48/256/512/2048
+    # give 0.8946/0.9041/0.8925/0.8900 for the same integrand.  Integrate panel by
+    # panel instead; n_gauss is now points PER PANEL, and 24 agrees with 96 to 3e-5.
+    edges = _TABLE_LAM[(_TABLE_LAM > b.lam_lo + 1e-9) & (_TABLE_LAM < lf - 1e-9)]
+    edges = np.concatenate([[b.lam_lo], edges, [lf]])
     nodes, weights = leggauss(int(n_gauss))
-    la = 0.5 * (lf - b.lam_lo) * nodes + 0.5 * (lf + b.lam_lo)
-    jw = 0.5 * (lf - b.lam_lo) * weights
+    la = np.concatenate([0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
+                         for lo, hi in zip(edges[:-1], edges[1:])])
+    jw = np.concatenate([0.5 * (hi - lo) * weights
+                         for lo, hi in zip(edges[:-1], edges[1:])])
     sig00 = np.array([b.matrix(cos_gamma, float(l))[0, 0] for l in la])
     G = order0_window(b.bg, la, lf, n_gauss=n_gauss)
     return float(np.sum(jw * sig00 * G ** 2))
```

**Patch the class, not an instance.** `sachs_mc_core._precompute` (line 119)
and `_field_precompute` (line 241) each construct their own
`_ds.Sigma2Builder`, and `driver_stats._builder()` caches a module-level
singleton behind `sigma2_matrix`, `sigma2_blocks`, `sigma2_6x6` and
`order0_mc`. Assigning `grid.builder` after `build_grid` (the pattern one of
the assessments used) repairs the `fk_expect_exact` path and nothing else:
`simulate`, `simulate_crn`, `simulate_ff_crn`, `simulate_fk_vr` and
`simulate_fk_pathA` would all keep the stock builder, and the figure would mix
two different `Sigma2`.

**Cost.** 210 `corr_op.C_fn` evaluations per `cos_gamma` instead of 160 (+31%),
paid once per gamma and cached. No change to any public signature.

**Housekeeping.** `psd_fix/sigma2_repaired.py`'s `_Base` and `R3Smoothed` read
`self._nodes`; they are diagnostic-only and would need a one-line update after
the patch lands.

---

## 3. Re-anchoring

`ANCHOR_RATIO = 0.881` is the median of
`order0_mc(apply_c0=False) / O0_analysis3`. It is traceable to `order0_mc`'s
**default `n_gauss = 48`** with the stock builder: I measure 0.880343 on the
ten-gamma grid **[V]** (an independent pass gets 0.88118 on the 26-gamma
analysis-3 grid). The three relevant numbers, all measured on the same grid
**[V]**:

| rule | stock | repaired | repaired/stock |
|---|---|---|---|
| GL 48 (today's default) | 0.880343 | 0.894567 | +1.6157% |
| GL 256 (what `simulate` prints) | 0.889792 | 0.904123 | +1.6106% |
| GL 512 | 0.889637 | 0.892509 | +0.3229% |
| GL 2048 | 0.889618 | 0.890027 | +0.0460% |
| **panelised, 24 / 48 / 96 per panel** | 0.889647 / 0.889626 / 0.889618 | **0.889357** (all three) | **-0.0303%** |
| uniform trapezoid, N = 4000 | 0.889618 | 0.889681 | +0.0071% |

**The famous "+1.61%" is a quadrature artefact, not a physics shift.** It is
one sample of an oscillating sequence: change only `n_gauss` and the same two
builders differ by -2.07% (128), +1.61% (256), +0.32% (512), +0.046% (2048).
Under a rule that breaks at the table nodes, self-converged to six digits, the
shift is **-0.030%**, fifty-five times smaller and of the opposite sign. The
stock builder only *looked* stable under Gauss-Legendre because its `C^2`
spline smooths the integrand into something a global polynomial rule can
integrate; it integrates the wrong density accurately.

**What to do.**

* **The repair alone requires no re-anchoring.** It moves the converged
  constant by -0.030%, which is thirty times smaller than the +0.95% by which
  the hardcoded 0.881 is *already* away from the converged stock value
  (0.889618). Do not re-derive `ANCHOR_RATIO` "because of the repair", and in
  particular never adopt 0.895 or 0.904: those are the artefact.
* **If you land the `order0_mc` panelisation in the same commit (recommended),
  re-derive it once**: `ANCHOR_RATIO = 0.8894`. That +0.95% change is driven by
  the quadrature fix, and it applies to the stock builder too.
* **What breaks if re-anchoring is skipped.** Nothing hard fails.
  `driver_stats.__main__`'s Part-A self-check tolerates 10% on the median and
  5% on the spread, so it keeps printing PASS while its median drifts from
  ~1.000 to ~1.010. The anchored `Sigma2` keeps a ~+0.95% normalisation error
  that has been there since the constant was fitted with the same unconverged
  rule that `order0_mc` still uses, so the analytic self-check hides it. The
  repair neither creates nor materially changes it (-0.03%). Skipping is
  therefore acceptable for this paper and wrong for the next one.
* **Do not fold the MC's own grid discretisation into the constant.** A
  trapezoid over the MC's production lambda grid (N = 1000, uniform) gives
  0.889617 for the stock builder **[V]**, i.e. the converged value, so the
  ~+0.8% by which the *realised* MC exceeds `order0_mc` comes from the AR(1)
  recursion's discretisation, not from the lambda rule. That is an N-dependent
  convergence property and must not be absorbed into an N-independent constant.

**One genuine cost of the repair, at the constant's own precision.** Under a
uniform rule the repaired integrand is discontinuous, so its error is O(h) with
a sign that depends on where the nodes fall: trapezoid N = 1000 gives 0.888327
and N = 4000 gives 0.889681 **[V]**, a 0.15% spread, where the stock builder
sits steady at 0.88962 for every N. The repaired builder is right on average
and converges; the stock one is smoothly and consistently wrong. Anyone quoting
a single MC anchor number after the repair will see a 0.1-0.3% grid dependence
that was not there before. Panelling removes it for the analytic reference.

---

## 4. Downstream

### Must be redone

| item | why |
|---|---|
| `driver_stats.__main__` Part-A self-check | its median moves; re-read the printed value rather than trusting PASS |
| `ANCHOR_RATIO` / `ANCHOR_C0` | only if the `order0_mc` panelisation lands (see section 3) |
| any cached `order0_mc` curve stored as a reference | recompute under the panelised rule |

### Need not be redone

Nothing in `PAPER_PROPOSAL.md` changes, and the figure does not need
regenerating.

* **`fk_reference` does not move at all**, at any gamma, to machine precision
  (agreement 6.7e-16 between stock and repaired). Note the *stated* structural
  reason is wrong as written: `cum3_from_Q(solve_Q(M, Z), M) = sym(Z)` holds
  only for positive-definite `M`, and it fails with 100% error at exactly the
  two indefinite nodes under the *nominal* calibration, where `fk_reference`
  does move by +0.125%. It is invariant on the production path because the
  AR(1) **smearing** restores PSD (calibration matrix `M`: 0/1000 non-PSD even
  though `V` is 2/1000). Right conclusion, wrong reason.
* **The exact FK expectation** at gamma = 1', sigma = 8, N = 1000, smeared:
  1.5547050e-05 -> 1.5546949e-05, a shift of -6.5e-6 relative.
* **The marker table.** Extrapolated ratio shifts +0.00013, +0.00012, +0.00011,
  +0.00010, +0.00008, +0.00008, +0.00005 across 0.5' to 17.3', against quoted
  errors of +/- 0.006 to +/- 0.069. The exact-expectation references 1.0005 /
  1.0010 / 1.0016 / 1.0058 become 1.0006 / 1.0011 / 1.0017 / 1.0059.
  "1.002 +/- 0.006" stays "1.002 +/- 0.006".
* **The gamma shape of every two-point observable.** Under a converged rule the
  repair changes kappa-kappa, xi_+, xi_- and kappa-gamma by the same -0.03%,
  in both the within-ray and the cross-ray block, from 0.1' to 600'; the
  spread of that change over gamma is 0.021-0.035%.
* **The p = 2 lever-arm anchor** cited in `PAPER_PROPOSAL.md` section 0. On the
  25-point analysis-3 gamma grid (0.633' to 183.262') the standard deviation of
  the normalised ratio is 0.0141 for the stock builder and 0.0140 for the
  repaired one, peak-to-peak 0.0705 -> 0.0701 **[V]**. The discriminant against
  p = 1 and p = 3 is untouched.
* **The FF channel**, bounded at -0.24% +/- 0.28% (gamma = 1', 32 paired seeds,
  n_real = 48000), with a noise-free linear-in-`Sigma2` companion predicting at
  most 0.4% at the production grid. Bounded, not resolved.
* **The eigen-floor's effect on what the MC realised**, +0.0029%, and upward:
  the floored `Sigma2_00` was negative, so flooring it to zero *raised* the
  realised Order-0. The repair removes a tiny upward bias rather than
  introducing one.

---

## 5. Should the corr_op table be rebuilt on a denser lambda grid?

**No. Do not rebuild.** But the defect the rebuild would address is real and is
larger than everything above.

**The measurement.** Against four dense rebuilds already on disk
(`psd_fix/dense/`), the production 21-node bilinear diagonal is wrong mid-cell
by **[V]**:

| cell | u = 0.25 | 0.50 | 0.75 | shared-node agreement |
|---|---|---|---|---|
| 550-700 (h = 150) | -20.20% | -24.00% | -16.12% | -0.003% / +0.128% |
| 1300-1447 (h = 147) | -12.78% | -15.80% | -10.60% | +0.370% / +0.164% |
| 1700-1900 (h = 200) | -15.40% | -15.01% | -8.17% | +0.019% / +0.010% |
| 2160-2250 (h = 90) | -5.57% | -1.49% | +2.19% | +0.222% / +0.074% |

The cause is structural: `C(lam1, lam2)` has a ridge on the diagonal (from the
`theta(lambda - chi)` window), so a product-grid bilinear interpolant
systematically *undershoots* the diagonal. `Sigma2Builder` queries precisely
that locus, where the error is maximal and one-signed; the analysis-3 folds
integrate over the whole cell, where the +4% far corners partly cancel the
-16% ridge (net 2-D cell-integral error -6.4% to +6.8%).

**Why refining the table is the wrong lever.** Halving the cell (21-node
h = 147 versus a 3-node rebuild h = 73.5, both against the 5-node truth at
lambda = 1336.75) reduces the mid-cell error only from -12.78% to -8.53%, a
factor **1.50x**, not the 4x a smooth function would give **[V]**. Reaching ~1%
therefore needs 16-60x more nodes, i.e. a 1-11 GB table (the stored size is
3 x 515 x N^2 x 8 bytes) and, at the measured ~N^1.8 build scaling from 172 s
(N=3) / 320 s (N=5) / 4264 s (N=21), weeks of compute. Switching
`RegularGridInterpolator` to `method="cubic"` would ring across the ridge
rather than undershoot it.

**The cheap lever, and evidence it is a real correction rather than a different
convention.** The table already stores `C_ab(lam_i, lam_i)` *exactly* at all 21
nodes. A one-dimensional sign-preserving PCHIP of `log|C|` through those values
(`psd_fix/t12_diag_aware_builder.py::R4DiagAware`) reproduces the dense truth to
<= 0.83% and is PSD everywhere tested. It moves the Order-0 anchor median by
**+6.04% (GL 256), +6.06% (panelised), +6.06% (trapezoid N = 4000)** -- the
same number under all three rules, because its density is smooth **[V]**. And,
decisively, it makes the MC track analysis-3's *gamma shape* much better: on the
25-point analysis-3 grid the normalised anchor ratio has

| builder | median | std(r / median) | peak-to-peak |
|---|---|---|---|
| stock | 0.89027 | 0.0141 | 0.0705 |
| R1 (repaired) | 0.89001 | 0.0140 | 0.0701 |
| R4 (diagonal-aware) | 0.94846 | **0.0052** | **0.0185** |

**[V]** -- a 2.7x tighter spread and a 3.8x tighter range, with the far tail
(144.7', 183.3') going from 0.9743 / 0.9321 to 1.0041 / 1.0120. Since
`O0_analysis3` is folded from the *same* table, a diagonal patch that were
merely "different" would degrade that agreement; it improves it sharply. This
is independent evidence that the bilinear diagonal, not the fold, is what is
wrong.

**Priority.**

* **P1, now, hours:** R1 + the `order0_mc` panelisation. Zero paper impact,
  removes a wrong-sign covariance. This is the commit to make.
* **P2, before any work that quotes the MC's Order-0 normalisation as physics
  rather than as an anchored quantity:** fix the equal-time diagonal. Do it
  **upstream in canoes**, not in `driver_stats`, by reparameterizing the
  interpolation in (mean, difference) coordinates so the ridge lies along a
  grid line (the trick already adopted for `chi_delta`), so that the analysis-3
  2-D folds benefit from the same fix. `R4DiagAware` is the validated
  driver-side stopgap. Requires re-anchoring by ~6% (`ANCHOR_RATIO` 0.889 ->
  0.948, `ANCHOR_C0` down by ~6%), which is why it is not a same-commit change.
* **P3, do not do:** rebuild the table on a denser lambda grid.

---

## 6. Residual risks and what is still unverified

1. **The whole exercise is about the interpolant, not the propagator.** Every
   number in sections 1-4 is an exact property of the bilinear interpolant that
   corr_op serves. Section 5 says that interpolant is itself wrong by -5% to
   -24% on the locus `Sigma2Builder` queries. R1 versus stock is a 0.03%
   argument inside a 6% error. That ordering should govern how much confidence
   is placed in the MC's absolute Order-0 normalisation, though not in the FK
   ratios the paper actually quotes, which are anchored.
2. **R1's own pointwise residual.** Reported as up to 1.05% just after a cell
   start, against an analytic evaluator, argmax at node + 0.65 Mpc; the likely
   cause is not-a-knot end conditions on the few samples inside a cell. I did
   not re-measure it and did not try a fix. A guard band of 1 Mpc around nodes
   brings it to 0.108%. It does not affect the integrated quantities.
3. **Node values are a convention.** The density is genuinely two-valued at a
   table node; the patch takes the left limit. Any single number quoted *at* a
   node inherits that choice.
4. **`corr_op`'s batch and scalar paths are not interchangeable in one
   process.** `C_fn_batch` caches an angular-channel table keyed on a rounded
   `(cos gamma, psi1, psi2)`, so whichever path touches a `cos` first fixes what
   both return; the difference is up to +0.165% at gamma = 1'. Every number here
   is from the scalar path (which is what both builders use). This is worth
   reporting upstream to canoes; I did not re-verify it.
5. **`psd_fix/sigma2_repaired.py::R0Reference` has an endpoint bug.** At exactly
   `lam == lam_lo` the cell lookup degenerates (`lo == hi`), `h` collapses to
   1e-3 and `np.clip` is called with `a_min > a_max`, so it evaluates outside
   the table and returns 3.117e-09 where the density is 2.375e-14, a factor
   1.3e5 **[V]**. Gauss-Legendre never lands on the endpoint, but any uniform
   grid that includes it does. `R0Reference` is diagnostic-only, but anyone
   reusing it must fix this first.
6. **The FF channel is bounded, not measured.** -0.24% +/- 0.28% at gamma = 1'
   with 32 paired seeds; a 4-seed first pass gave the opposite sign at 1.3
   sigma. The usual warning about small seed blocks applies.
7. **The 2-D off-diagonal bias is untested downstream.** The analysis-3
   O0/FF/FK products are folded from the same bilinear table over the full
   (lam1, lam2) domain, where the cell-integral error is -6.4% to +6.8%. Nobody
   measured whether those products carry a comparable bias. That needs its own
   experiment before section 5's P2 is acted on.
8. **The dense "truth" of section 5** is a rebuild with a coarser
   source-adaptive chi grid than production (5 sources rather than 21). Bounded
   three ways: shared-node agreement 0.003-0.37% **[V]**, two independent dense
   builds agreeing to +0.204% at a shared node **[V]**, and n_chi 256 -> 512
   moving a node value by -0.0004% **[V]** (0.50% mid-cell, reported). Treat the
   mid-cell errors as +/- 0.5%.
9. **Not tested at all:** `R2Density` and `R3Smoothed` downstream; source
   distances other than the baseline `lam_f`; `sigma_lambda` other than 8 (FK)
   and 0 (FF); the "full" nonlinear-arm FK flavour, as opposed to the T1+T2
   perturbative one the markers use.

---

## Appendix: what I re-measured, and how

```
export PYTHONPATH=/Users/zzhang/projects/angular_statistics/canoes/src
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
cd driver_field_emulators/code/mc_fk_complete

$PY psd_fix/v1_verify.py        # bilinear structure; derivative sign flip at 1300;
                                # wrong-sign Sigma2_00 and the negative-definite 6x6
$PY psd_fix/v2_anchor_rules.py  # d6 anchor reproduced (0.889792 / 0.904123), then
                                # re-integrated under 9 rules: the +1.61% is quadrature
$PY psd_fix/v3_dense.py         # production bilinear vs the four dense rebuilds
$PY psd_fix/v4_r4.py            # R4 anchor +6.04/+6.06/+6.06%; R0 endpoint bug
$PY psd_fix/v5_shape.py         # gamma-shape spread, 25-point analysis-3 grid
$PY psd_fix/v6_scaling.py       # halving h buys only 1.50x
$PY psd_fix/v7_patch_test.py    # the proposed patch == R1Knots bitwise; PSD scan
$PY psd_fix/v8_order0_patch.py  # the proposed order0_mc, converged at 24 pts/panel
```

Assessments I relied on without re-measuring: the downstream FK numbers
(exact expectation, markers, `fk_reference`, Monte Carlo), the FF bound, the
`cum3_from_Q` identity analysis, the knot-placement ablation, and R1's
pointwise residual against an analytic evaluator.
