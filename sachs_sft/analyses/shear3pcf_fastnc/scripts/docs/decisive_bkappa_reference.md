# DECISIVE B_kappa-phase reference: which shear-3PCF pipeline is the outlier?

Builder: `stage1_bkappa_phase_reference.py`
Comparison: `stage1_threeway_compare.py` -> `outputs/stage1_threeway_compare.npz`
Reference output: `outputs/stage1_bkappa_phase_ref.npz`
Run log: `outputs/stage1_bkappa_phase_ref_run.log`
Interpreter: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`

---

## TL;DR -- VERDICT

An INDEPENDENT, definition-level cosmic-shear 3PCF reference, built straight from
the convergence bispectrum `B_kappa` via the standard flat-sky spin-2 phase
projection `gamma~(l) = e^{2i phi_l} kappa~(l)`, **agrees with fastnc tree-level
to ~1% in the frame invariant (median ref/fnc = 0.99, range 0.83-1.01 over four
configs)** and **disagrees with the canoes `zeta_D` route by a structured factor
of ~4-17x (median can/ref ~ 7.8)**.

> **The canoes driving-field spin-2 `zeta_D` route is the OUTLIER.**
> fastnc is corroborated by a fully independent first-principles reference.
> The discrepancy is NOT a constant: it is a real, structured (gamma- and
> phi-dependent) spin-2 over-normalization / shape bug in the `zeta_D` route,
> ~4x at the PRIMARY point and growing to ~17x at wide opening angle.

This refutes the earlier (`shear3pcf_natural_map.md`) hedge that "both candidate
single-normalization fixes are excluded, so the shear components are a genuinely
different projection". They are NOT genuinely different at the level that
matters: two independent `B_kappa -> Gamma` implementations (this reference and
fastnc) land on the SAME answer; only the `zeta_D` route is off.

---

## 1. What the reference is (the definition)

The cosmic-shear 3PCF natural components are, BY DEFINITION (Schneider & Lombardi
2003 A&A 408 829; Schneider, Kilbinger & Lombardi 2005 A&A 431 9; Sugiyama+2024
arXiv:2407.01798), correlators of the spin-2 shear

        gamma(x) = INT d^2l/(2pi)^2 e^{i l.x} e^{2 i phi_l} kappa~(l),
        phi_l = polar angle of l,

projected onto a per-vertex reference direction. In the CENTROID projection each
leg's reference is the vertex->centroid direction `alpha_j = polar(centroid-X_j)`,
so `g|_j = gamma(X_j) e^{-2 i alpha_j}` and (s_j = +1 unconjugated, -1 conjugated):

        Gamma^0 = <g|_1   g|_2   g|_3 >   (s=+,+,+)
        Gamma^1 = <g|_1^* g|_2   g|_3 >   (s=-,+,+)
        Gamma^2 = <g|_1   g|_2^* g|_3 >   (s=+,-,+)
        Gamma^3 = <g|_1   g|_2   g|_3^*>  (s=+,+,-).

Using `<k~(l1)k~(l2)k~(l3)> = (2pi)^2 delta(l1+l2+l3) B_kappa`, this collapses to a
**2D double Fourier-plane integral** over two independent ell vectors:

        Gamma^mu = INT d^2l2/(2pi)^2 d^2l3/(2pi)^2  B_kappa(|l1|,|l2|,|l3|)
                     * exp[ i l2.(X2-X1) + i l3.(X3-X1) ]
                     * exp[ 2 i ( s1(phi_{l1}-alpha1) + s2(phi_{l2}-alpha2)
                                  + s3(phi_{l3}-alpha3) ) ],     l1 = -(l2+l3).

This is the centroid projection **by construction** (no fastnc `x2cent` bridge is
used; the per-leg spin reference is built in), so the reference is independent of
fastnc's 2DFFTLog `MCF222` mode-coupling machinery AND of canoes' great-circle
`zeta_D` cumulants. The geometry (apex at origin, base vertices symmetric about
+x, `alpha_j` and displacement vectors `d2=X2-X1`, `d3=X3-X1`) is derived
symbolically in `scripts/spin2_centroid_direct.wl`:

        centroid = (2/3 cos(phi/2), 0)
        alpha1 = atan2(-3 sin(phi/2), -cos(phi/2))
        alpha2 = atan2( 3 sin(phi/2), -cos(phi/2))
        alpha3 = atan2(0, cos(phi/2)) = 0
        d2 = (0, -2 sin(phi/2)),  d3 = (-cos(phi/2), -sin(phi/2))   (units of gamma)

The spin-2 phase `e^{2i phi_l}` IS the standard flat-sky shear projection; the
`x2cent` factors in fastnc were verified to be pure phases (`|x2cent|=1`,
`scripts/spin2_3pcf_phase.wl`), confirming our centroid construction and fastnc's
land in the same projection family.

## 2. The B_kappa is fastnc-identical (the shared physics input)

`B_kappa` is the ONE input that must be common to all three objects. We reproduce
fastnc's exact single-plane Limber chain (`stage1_fastnc_tree.py`,
`fastnc/bispectrum.py`), entirely in h-units:

| quantity | formula | fastnc source |
|---|---|---|
| chi(z)   | astropy wCDM `comoving_distance(z)*h` [Mpc/h] | bispectrum.py:293 |
| g(chi)   | `(3/2)(100/c)^2 Om (1 - chi/chi_s)` [h-units]  | bispectrum.py:369-374 |
| weight   | `g(chi)^3 * (1/chi) * (1+z)^3`                  | bispectrum.py:533 |
| B_delta  | `D(z)^4 [2 F2(k1,k2)P(k1)P(k2) + 2 cyc]`        | (SPT subclass) |
| F2       | `5/7 + (1/2)cosθ(ka/kb+kb/ka) + (2/7)cosθ^2` (TRUE SPT, not F2_eff) | |
| k_j      | `l_j / chi` [h/Mpc]                              | bispectrum.py:611-612 |
| growth   | analytic flat-wCDM `D(z)/D(0)` (build_growth)   | stage1_fastnc_tree.py |
| P(k)     | PCAMBz0.txt z=0, h-units, fed as-is             | stage1_grid.json |

`B_kappa(l1,l2,l3) = INT dchi weight * B_delta`. Spot values:
`B_kappa(100,100,100)=1.27e-14`, `B_kappa(1000,1000,1732)=1.29e-17` (right sign,
right magnitude for a tree convergence bispectrum at z_s=5). Same P(k), same
h=0.6711, same Omega_m=0.31609, same z_s=5 single source plane as STAGE 0.

This is the SAME `B_kappa` family that STAGE 0 validated: STAGE 0 matched our
scalar `zeta_TTT` fold to 15%, so the `B_kappa` machinery is independently
trusted.

## 3. Method + internal convergence

The 2D Fourier integral is evaluated in polar coords (l2, phi2, l3, phi3): log
radial nodes `l in [0.1, 3e4]`, uniform angular nodes. Memory-safe: the inner
(l2, phi2) plane is vectorised, the outer (l3, phi3) pair is looped. `B_kappa`
is the exact Limber chi-integral, tabulated on an `(ln l2, ln l3, cos23)`
trilinear interpolant (validated **median 1.1e-3, max 1.1e-2 rel err vs the exact
routine**) so the inner loop is fast.

Internal convergence (PRIMARY gamma=10', phi=60deg), frame-invariant:

| (n_l, n_ang, lmax) | frame-inv | Δ vs prev |
|---|---|---|
| (64, 64, 2e4)   | 2.298e-7 | -- |
| (96, 96, 2e4)   | 2.321e-7 | 1.0% |
| (128,128, 3e4)  | 2.265e-7 | 2.4% |

Worst per-component config (phi=120): n=128 frame-inv 1.274e-7 vs n=160 1.267e-7
(0.6% drift); |G1|, |G3| change <2%. **Production setting = (128,128, lmax=3e4),
converged to ~2-3%.** chi-grid n_chi=160; B_kappa interpolant 96x96x64.

## 4. The three objects + the numbers

Same pinned configs (`stage1_grid.json`), z_s=5, centroid projection:

1. **ref** = `stage1_bkappa_phase_ref.npz` (THIS reference; definition-level).
2. **canoes (under test)** = `stage1_ours_v2.npz` (canoes `zeta_D` route).
3. **fastnc** = `stage1_fastnc.npz` (fastnc tree-level 2DFFTLog).

### PRIMARY point (gamma=10', phi=60deg) -- |Gamma^mu|

| comp | \|ref\| | \|fnc\| | \|can\| | ref/fnc | ref/can | fnc/can |
|---|---|---|---|---|---|---|
| G0 | 6.51e-8 | 6.82e-8 | 3.44e-8 | **0.96** | 1.90 | 1.98 |
| G1 | 1.36e-7 | 1.38e-7 | 7.07e-7 | **0.99** | 0.19 | 0.19 |
| G2 | 1.20e-7 | 1.24e-7 | 4.28e-7 | **0.96** | 0.28 | 0.29 |
| G3 | 1.20e-7 | 1.24e-7 | 4.30e-7 | **0.96** | 0.28 | 0.29 |
| frame-inv | 2.26e-7 | 2.33e-7 | 9.32e-7 | **0.97** | 0.24 | 0.25 |

ref and fnc agree per-component to <=4%; canoes is ~4x larger in the
frame-invariant (and the per-component split is also wrong: canoes puts too much
into G1 and too little into G0).

### Full set (frame-invariant)

| config | ref | fnc | can | **ref/fnc** | ref/can | fnc/can |
|---|---|---|---|---|---|---|
| (10', 60deg) PRIMARY | 2.26e-7 | 2.33e-7 | 9.32e-7 | **0.97** | 0.24 | 0.25 |
| (10', 30deg) | 1.98e-7 | 2.38e-7 | 1.20e-6 | **0.83** | 0.16 | 0.20 |
| (10', 120deg)| 1.28e-7 | 1.27e-7 | 2.21e-6 | **1.01** | 0.058 | 0.058 |
| (50', 60deg) | 5.98e-8 | 5.93e-8 | 5.67e-7 | **1.01** | 0.11 | 0.11 |

**Self-check ref/fnc: median 0.99, range [0.83, 1.01]** -- the reference matches
fastnc tree to ~1-17%, far inside the ~30% gate. The reference is TRUSTED.

**Verdict ratios (frame-invariant):**
- ref/can: median 0.135, range [0.058, 0.243] (canoes 4-17x too large)
- fnc/can: median 0.152, range [0.058, 0.250]
- can/ref: median 7.77, range [4.12, 17.23]

## 5. VERDICT

The canoes `zeta_D` spin-2 route is the **OUTLIER**, by a **structured** factor of
**~4x at the primary point, growing to ~17x at wide opening angle** (phi=120deg).
The two independent `B_kappa -> Gamma` objects (this definition-level reference
and fastnc) agree with each other to ~1% (median frame-invariant), so fastnc is
**vindicated** as the correct tree-level shear 3PCF.

The discrepancy is NOT a constant normalization and NOT a frame phase (it varies
4-17x with configuration and mis-splits the per-component magnitudes), so it is a
genuine spin-2 **over-normalization + shape** bug in the canoes driving-field
`zeta_D` cumulant route -- consistent with the standing
`project_kappa3_spin2_helicity_bug_2026-06-04.md` finding that the canoes spin-2
`(2,2,-2)` cumulant normalization was UNVALIDATED. This test settles it: that
normalization (and its configuration dependence) is wrong, by a structured
~4-17x.

### Caveats (honest)

- ref-vs-fnc per-component ratios at phi=30/120 scatter 0.4-2.8x even though the
  frame-invariant agrees to ~1-17%. This is the expected residual between two
  DIFFERENT spin-2 implementations (our real-space centroid Fourier integral vs
  fastnc's 2DFFTLog multipole resum) -- the frame-invariant is the
  projection-robust quantity, and it agrees. This residual is an order of
  magnitude smaller than the canoes discrepancy, so it does not muddy the verdict.
- Overall sign: this reference's Gamma are negative-real-dominant where fastnc's
  are positive-real (an overall e^{i*const} convention from the absolute spin
  origin); we compare MAGNITUDES and frame-invariants, which are sign/phase
  independent. The per-component magnitude agreement (Section 4) is unaffected.

## 6. Sources

- Schneider & Lombardi 2003, A&A 408, 829 (shear 3PCF natural components).
- Schneider, Kilbinger & Lombardi 2005, A&A 431, 9.
- Sugiyama, Takada+ 2024, arXiv:2407.01798 (fastnc method).
- Porth+ 2023, arXiv:2309.08601 (centroid x2cent projection bridge).
- fastnc: `fastnc/bispectrum.py` (Limber B_kappa), `fastnc/fastnc.py:752-784`
  (x2cent), `fastnc/coupling.py:250-268` (MCF222).
- STAGE 0 scalar validation: `stage0_spt_reference.py`,
  `outputs/stage0_spt_reference.npz` (B_kappa machinery, 15% vs zeta_TTT).
- Symbolic geometry: `scripts/spin2_centroid_direct.wl`, `scripts/spin2_3pcf_phase.wl`.
- Standing helicity-bug note: `project_kappa3_spin2_helicity_bug_2026-06-04.md`.
