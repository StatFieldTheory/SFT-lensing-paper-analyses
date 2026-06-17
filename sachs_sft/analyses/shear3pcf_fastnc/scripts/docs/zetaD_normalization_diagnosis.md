# zeta_D spin-2 modulus normalization: root-cause diagnosis

Read-only investigation. Probe: `scripts/probe_zetaD_normalization.py`
(env `sft-wick`). Symbolic anchor: `scripts/derive_spin2_gaunt_norm.wl`.
Outputs: `outputs/zetaD_normalization_probe.npz`.

## TL;DR verdict

**The sharp hypothesis is REFUTED.** The deployed `zeta_D` is NOT
over-normalized by a spin-independent `(0,0,0)`-anchor error. Two independent
facts kill that hypothesis:

1. The `(0,0,0)` anchor `w0 = wigner_3j(L1,L2,L3;0,0,0)` is applied
   **identically** to the spin-2 channels `(2,2,2)`/`(0,2,2)` and to
   `(2,2,-2)`/`(0,2,-2)`. The `(2,2,2)`/`(0,2,2)` path **reproduces the
   validated scalar `zeta_PPP`/`zeta_TPP` to rel-err < 1e-6** (49/49
   `test_kappa3_mod_normalization.py` pass, re-run this session). So `w0` is the
   deliberate, validated reduced->full convention for **all** channels, not a
   spin-2-specific bug.

2. There is **no** spin-weighted Gaunt symbol that could "correctly" replace
   `w0` for this channel. The modulus channel `(2,2,-2)` has total spin
   `+2 != 0`; it is not a rotational scalar. The naive spin anchor
   `wigner_3j(L1,L2,L3; -2,-2,+2)` **vanishes identically** (bottom-row sum
   `= -2 != 0`) -- confirmed symbolically (`derive_spin2_gaunt_norm.wl`) and
   numerically (the brute-force "spin anchor" gave exactly 0 at every config).

**What the data actually shows:** the deployed `zeta_D` is, if anything,
*under* the appendix anchor, not over it. At near-squeezed triples
`|zeta_D|/|zeta_TTT| = 0.15` and `0.006` -- it is *smaller* than `zeta_TTT`,
never a "few-x" larger. (At the literal `cos=(1,1,1)` row in the real-builder
data, `theta_ij = 0`, i.e. all three lines of sight coincide -- a degenerate
point where `zeta_D` is exactly 0; at a *finite* equilateral triangle, e.g.
g=10', `zeta_D (2,2,-2) = -1.3e-7` is nonzero. The companion modulus
`zeta_B (0,2,-2) = +8.1e-3` there is comparable to `zeta_TTT ~ 0.02`,
consistent with the appendix `zeta_B -> zeta_TTT`. It is specifically the
`(2,2,-2)` `D`-channel that comes out small, ~1e-7, not large.)

Therefore the fastnc ~3.5-5x (frame-invariant median ~8x) **over-LARGE** shear
natural components `Gamma^1..3` are **not** caused by an over-normalized
`zeta_D` cumulant. The residual is the genuine, configuration-dependent
projection-convention difference already documented in
`docs/shear3pcf_natural_map.md` (canoes great-circle real-space spin-2 cumulants
vs fastnc fixed-frame 2DFFTLog `MCF222` projection of the scalar convergence
bispectrum). No single normalization factor fixes it.

## The traced deployed code path (file:line)

The deployed `equal_time_limber` `zeta_D` is built by:

- **LOW** `compute_kappa3_mod_zeta_equal_time`
  (`canoes/src/canoes/sachs/kappa3.py:3737`), spins `_spin_weights_dmod=(2,2,-2)`
  (`kappa3.py:3755`), projected by `SpinAwareZetaEvaluator(..., spin_weights=(2,2,-2))`
  (`kappa3.py:3813-3816`).
- **HIGH** `compute_kappa3_mod_sigma3_high` (`kappa3.py:3912`), spins at
  `kappa3.py:3927`, projected by `_kappa3_limber_alpha_kernel`
  (`kappa3.py:1747`) via `i^S exp(iS arg(Q)) J_S(|Q|)` -- a correctly
  spin-threaded flat-sky orientation integral (no `w0` factor; spin-correct).

The over-normalization candidate is the LOW angular synthesis. Its kernel build
is `SpinAwareZetaEvaluator._build_kernel_numpy`
(`canoes/src/canoes/nuell/correlation/_spin_aware_three_pt.py:565-693`). The
per-triple coefficient (lines 598-604) is:

```
h     = sqrt[(2L1+1)(2L2+1)(2L3+1)/(4 pi)]      # line 598
sq_l1 = sqrt[(2L1+1)/(4 pi)]                      # line 602
w0    = wigner_3j(L1,L2,L3; 0,0,0)               # via _filter_active_triples (line 587)
coef  = h * w0 * sq_l1 * (-1)^s1                  # line 604
```

times the spin m-sum (lines 618-691):

```
M_s = Re Sum_{m2} wigner_3j(L1,L2,L3; -s1, m2, s1-m2)
            * {s2}Y_{L2 m2}(g12, 0) * {s3}Y_{L3, s1-m2}(g13, phi3)
```

**Yes, it uses the `(0,0,0)` anchor `w0` for the spin-2 channel.** But this is
the *same* `coef` form the scalar `ZetaEvaluator` uses
(`three_pt.py:567-572`: `coef = h * w0 * sq_l1`), differing only by the
spin-1 pole phase `(-1)^s1`. The docstring states this explicitly
(`_spin_aware_three_pt.py:25-28`): "the W3j0 anchor is kept as a parity selector
for all channels". The m-sum with `m1 = -s1`, `m1+m2+m3=0` is the correct
leg-1-at-pole reduction `{s1}Y_{L1,m1}(pole) ∝ delta_{m1,-s1}` and **is**
spin-threaded.

## Why the `(0,0,0)` anchor is correct here (not the bug)

The canoes reduced bispectrum `b_{L1L2L3}` from `compute_bispectra_deri` is the
Komatsu-Spergel **isotropic reduced bispectrum** (the full bispectrum is
`B = h * (L1 L2 L3; 0 0 0) * b`, see `nuell/utils/wigner.py:86`
`full_bispectrum`). The reduced->full anchor `(0,0,0)` belongs to the *scalar
shape function*; the spin structure of the cumulant lives entirely in the
**configuration synthesis** -- the spin-weighted `{s}Y_{Lm}` tables and the
`wigner_3j(...; -s1, m2, s1-m2)` m-sum. This is exactly how the paper's own
derivation frames it (`scripts/mathematica/derive_kappa3_cross_cumulants.wl:66-71`
keeps `(L1 L2 L3; 0 0 0)` as the parity selector for all channels and puts the
spin in the `{s}Y` projection). The `(2,2,2)==PPP` and `(0,2,2)==TPP`
self-consistency tests are a direct proof that this synthesis is normalized
correctly for spin-2 legs.

## What the probe measured (with actual numbers)

### Per-triple W3j (`derive_spin2_gaunt_norm.wl`)

Along the diagonal `(L,L,L)` the naive spin Gaunt `(L L L; -2,-2,+2)` is **0**
for every `L` (Mathematica: "not physical", bottom row sums to `-2`). The
`B`-channel naive anchor `(L L L; 0,-2,2)` is finite with ratio
`(000)/(0,-2,2) -> -2` as `L` grows (so any "fix" that swapped anchors would
*rescale* `B`, not `D`, and by a sign-flipping factor -> -2, not +few). This
already shows the "swap `w0` for a spin anchor" idea is ill-posed for `D`.

### Deployed evaluator vs the appendix small-angle anchor

(`probe_zetaD_normalization.py`, bare reduced `b`, ell_max=16, phi=90):

| gamma  | \|zeta_D\| | \|zeta_TTT\| | \|D\|/\|TTT\| |
|-------:|-----------:|-------------:|--------------:|
| 0.5'   | 1.00e-10   | 2.06e-02     | ~5e-9         |
| 1.0'   | 4.01e-10   | 2.06e-02     | ~2e-8         |
| 2.0'   | 1.60e-09   | 2.06e-02     | ~8e-8         |

`zeta_D ~ gamma^2` and is ~8-9 orders of magnitude **below** `zeta_TTT` at small
angle. (`zeta_PPP (2,2,2) ~ 1e-15`, the `d^L_{2,-2} ~ gamma^4` suppression.)
This is the OPPOSITE of an over-normalization. The appendix
`zeta_D -> zeta_TTT` limit requires the `COEFS_PSI0` per-leg eth^2 operator
`sqrt(L(L+1)(L-1)(L+2))/chi^2`, which the deployed builders DO apply via
`compute_bispectra_deri_psi0`.

### Real deployed builders (the apples-to-apples cumulants)

(`consistency_deployed_vs_fresh.npz`; `reldiff_Dmod = reldiff_TTT = 0`, i.e.
deployed == fresh callable, confirming this is the deployed `zeta_D`):

| triple `cos`              | sum\|zeta_TTT\| | sum\|zeta_D\| | \|D\|/\|T\| |
|---------------------------|----------------:|--------------:|------------:|
| (1, 0.857, 0.857)         | 2.15e-13        | 3.18e-14      | 0.148       |
| **(1, 1, 1)** = theta=0   | 5.57e-11        | **0.0**       | **0.0**     |
| (0.857, 0.857, 0.714)     | 1.56e-13        | 9.23e-16      | 0.006       |

The `(1,1,1)` row is the degenerate coincident-point limit (`theta_ij=0`) where
`zeta_D` is identically 0; the near-squeezed finite triples have `zeta_D`
*below* `zeta_TTT`. There is no configuration where the deployed `zeta_D`
cumulant is "3.5-5x too large".

## Where the fastnc over-LARGE Gamma actually comes from

The downstream map (`derive_shear3pcf_natural_components.wl`,
`stage1_ours_v2.py`) is `Gamma^mu = Z_channel(mu) * pure_phase`, so
`|Gamma^mu| = |fold(zeta_channel)|` is set entirely by the cumulant magnitude
and CANNOT be inflated by the projection phase. Since the deployed `zeta_D` is
not over-large, the 3.5-5x is structural:

- fastnc builds **all four** `Gamma^mu` from the **same scalar convergence
  bispectrum** `b_kappa` via the fixed-frame 2DFFTLog spin-2 mode coupling
  `MCF222` (`fastnc/coupling.py:250`, `bispectrum.py:746`).
- canoes builds them from **great-circle-frame** spin-2 driving-field cumulants
  with the exact-Wigner-3j (ell<=60) + Limber (ell<=1000) split.

These are genuinely different projections of the same physics. The doc
`shear3pcf_natural_map.md` already established (with a decisive cross-check) that
NEITHER a constant spin-2 factor times the scalar convergence 3PCF NOR a
constant times the great-circle spin-2 cumulants matches fastnc: the
frame-invariant ratio ranges 3-17x (median ~8), and the
`fnc / |fold(zeta_TTT)|` ratio ranges 0.02-25.6 (std/mean = 3.35). Both
"single-normalization" fixes are excluded by the data. My independent
cumulant-level probe corroborates this: a zeta_D renormalization is not the
lever, because zeta_D is not the thing that is too big.

## Honest residual / caveats

- The appendix statement `zeta_B -> zeta_TTT, zeta_D` finite-at-small-angle
  describes the **modulus** combination `<Phi|Psi0|^2>` / `<Psi0|Psi0|^2>` whose
  small-angle behaviour is carried by `d^L_{2,2} -> 1`. The deployed
  `(2,2,-2)`/`(0,2,-2)` SpinAware channels instead **vanish on the squeezed
  diagonal**. This is a real tension worth a second look: it suggests the
  great-circle `(2,2,-2)` SpinAware reduction does **not** realize the
  `d^L_{2,2}`-dominated modulus the appendix has in mind at small separation
  (it behaves like a `d^L_{2,-2}`-suppressed, total-spin-2 object on the
  diagonal). If a bug exists, it is here -- a possible spin-index / m-sum
  convention in the `(2,2,-2)` reduction making it the wrong helicity
  combination -- NOT an over-normalization. But this would make `zeta_D` too
  SMALL, which cannot explain a too-LARGE `Gamma`. So it is orthogonal to the
  fastnc magnitude discrepancy.
- I did not (read-only) attempt to fix anything. No factor was fit to fastnc.

## One-line "fix" descriptions (NOT applied)

- For the fastnc magnitude discrepancy: **no scalar/normalization fix exists**;
  closing it requires implementing fastnc's fixed-frame 2DFFTLog `MCF222`
  spin-2 projection on our side (a reimplementation of Sugiyama+2024), out of
  scope. The honest report stands as in `shear3pcf_natural_map.md`.
- For the separate (orthogonal) squeezed-diagonal tension: investigate whether
  the `(2,2,-2)` SpinAware m-sum/helicity assignment realizes the
  `d^L_{2,2}`-modulus the appendix expects (a convention check in
  `_spin_aware_three_pt.py:606-691`), but note this would make `zeta_D`
  *larger/finite*, not smaller, and still would not reduce `Gamma`.
```
