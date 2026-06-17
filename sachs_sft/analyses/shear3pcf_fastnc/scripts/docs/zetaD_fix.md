# zeta_D spin-2 kernel fix: pinned mechanism, minimal diff, verification gates

Branch (canoes): `fix/spin2-zetaD-2_2_-2-norm` (uncommitted, for human review).
Oracle: `scripts/probe_zetaD_kernel_oracle.py` -> `outputs/zetaD_kernel_oracle.npz`.
Env: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`.

## PHASE 1 -- the EXACT mechanism (pinned by an independent oracle)

The bug is NOT in the negative-spin `{s}Y` convention (the original narrowed
hypothesis). The oracle verified cell-by-cell that canoes' `_build_spin_y_table`
spin-weighted harmonic d-part is **exactly** `d^l_{m,-s}` (matches sympy
`Rotation.d` to ~1e-15 for spin 2), and that the per-`m2` Wigner-3j m-sum
(`_build_spin_wigner3j_m_sum`, value `(l1 l2 l3; -s1, m2, s1-m2)`) is correct.

The bug is in the **group-batched m-sum REDUCTION** inside
`SpinAwareZetaEvaluator._build_kernel_numpy`
(`canoes/src/canoes/nuell/correlation/_spin_aware_three_pt.py`).

Triples are grouped by `(l2, l3)` (sorted with `l1` ascending). The loop read the
per-group `m`-count and the W3j reshape from the **group-FIRST triple `g0`**:

```python
count_g = int(m_count[g0])          # line ~665, WRONG
if count_g == 0:
    K[:, g0:g1] = 0.0               # zeroes the ENTIRE (l2,l3) group
    continue
m2_start = int(m2_min[g0])
...
W_group = wm_flat[o_start:o_end].reshape(n_l1_in_group, count_g)   # line ~688
```

`m_count` is uniform across an `(l2,l3)` group EXCEPT that triples with
`l1 < |s1|` are masked to `m_count = 0` (the `|m1| = |-s1| <= l1` selection,
`_build_kernel_numpy` line ~614 `bad_m1 = np.abs(s1) > triples_l1`). Because the
group is sorted with `l1` ascending, the masked triples are a PREFIX, and the
group-first triple often has `l1 < |s1|`.

For `s1 = +/-2` channels this is catastrophic. Example (ell_max=6, group
`l2=l3=2`, triples `l1 in {0,2,4}`): the first triple is `l1=0 < |s1|=2`, so
`m_count[g0]=0`, and the loop hits `K[:,g0:g1]=0.0` -- **zeroing the entire group
including the valid `l1=2,4` triples**. Even when `count_g != 0`, the
`wm_flat[...].reshape(n_l1_in_group, count_g)` mis-shapes, because the flat slice
holds only the unmasked rows (masked triples contribute 0 flat entries), so
`o_end - o_start = n_valid * count_g != n_l1_in_group * count_g`.

### Oracle numbers (ell_max=6, two explicit triangles)

Independent full-`m`-sum reference (trusted sympy Wigner-d, NO group batching) vs
canoes production K, per active triple:

| channel  | cells the reference says nonzero but canoes ZEROES |
|----------|----------------------------------------------------|
| (2,2,2)  | 41 / 83                                            |
| (2,2,-2) | 60 / 120                                            |

So canoes was dropping ~half of every `s1!=0` channel's angular kernel -> the
LOW `zeta_D` came out far too SMALL (consistent with the prior read-only
diagnosis `zetaD_normalization_diagnosis.md`, which found `zeta_D ~ 1e-7`, BELOW
`zeta_TTT`, NOT over-large).

Affected channels: any with `s1 != 0`, i.e. **PPP `(2,2,2)`** and the modulus
**D-channels `(2,2,-2)`, `(-2,2,2)`, `(2,-2,2)`**. Channels with `s1 = 0` (TTT
`(0,0,0)`, TTP `(0,0,2)`, TPP `(0,2,2)`, zeta_B `(0,2,-2)`) are UNAFFECTED.

Why the existing `(2,2,2)` self-consistency test still passed before the fix:
the canonical PPP `zeta_PPP` is built with the SAME buggy `SpinAwareZetaEvaluator
(2,2,2)`, so the mod-path `(2,2,2)` and canonical PPP are IDENTICALLY wrong and
the self-comparison cancels the bug.

## PHASE 4 -- the minimal fix

File: `canoes/src/canoes/nuell/correlation/_spin_aware_three_pt.py`, group loop in
`_build_kernel_numpy` (was lines ~661-691).

BEFORE:
```python
count_g = int(m_count[g0])
if count_g == 0:
    K[:, g0:g1] = 0.0
    continue
m2_start = int(m2_min[g0])
...
o_start = int(offsets[g0]); o_end = int(offsets[g1])
n_l1_in_group = g1 - g0
W_group = wm_flat[o_start:o_end].reshape(n_l1_in_group, count_g)
```

AFTER (minimal):
```python
n_l1_in_group = g1 - g0
grp_counts = m_count[g0:g1]
count_g = int(grp_counts.max())          # natural count from an UNMASKED triple
if count_g == 0:
    K[:, g0:g1] = 0.0
    continue
m2_start = int(max(-l2, -l3 + s1))       # natural m-range, not group-first
...
# per-triple fill from each triple's own offset slice; masked rows stay 0.
W_group = np.zeros((n_l1_in_group, count_g), dtype=np.float64)
for i in range(n_l1_in_group):
    c_i = int(grp_counts[i])
    if c_i == 0:
        continue
    oi = int(offsets[g0 + i])
    W_group[i, :c_i] = wm_flat[oi:oi + c_i]
```

Properties:
- For `s1 = 0` (no masked triples) every `c_i == count_g` and the per-row fill
  reproduces the old contiguous `reshape` **bit-for-bit** -> all scalar/`s1=0`
  channels are an exact no-op (verified: `(0,0,0)` matches `ZetaEvaluator` to
  2.4e-15; `(0,0,2)/(0,2,2)/(0,2,-2)` bit-identical pre/post fix to 12 digits).
- For `s1 != 0` the masked prefix triples stay 0 and the valid `l1 >= |s1|`
  triples now contribute their full `m`-sum.

## VERIFICATION GATES (5)

### Gate 1 -- ORACLE (canoes K == independent full-m-sum reference)
After fix (`outputs/zetaD_kernel_oracle.npz`):

| channel  | max rel-err | abs-diff | zeroed cells |
|----------|-------------|----------|--------------|
| (2,2,2)  | 3.4e-4 (cancellation; max abs-diff 6.8e-18 on ~1e-12 values) | 6.8e-18 | 0/83 |
| (2,2,-2) | 1.5e-9      | 7.9e-18  | 0/120 |

PASS. The `(2,2,2)` "rel-err" is pure float64 cancellation: max ABSOLUTE
difference 6.8e-18 on kernels of magnitude ~1e-12 (the `d^l_{2,-2}~gamma^4`
suppression). Both channels are bit-level agreement.

### Gate 2 -- existing canoes tests STAY GREEN
`pytest tests/unit/sachs/test_kappa3_mod_normalization.py
tests/unit/correlation/test_spin_aware_three_pt.py` -> **71 passed** (same as
baseline 71). The +2 PPP/TPP self-consistency tests still pass (the canonical PPP
and mod-path PPP use the same now-fixed evaluator, so they change together). PASS.

### Gate 3 -- zeta_B (0,2,-2) still honors the appendix small-angle limit
`(0,2,-2)` has `s1 = 0` -> provably unchanged by the fix. Verified bit-identical
pre/post fix (12 sig figs) at 3 configs:
`zeta_B(0,2,-2) @ (10',60) = 8.081056745780e-03` (pre = post), comparable to
`zeta_TTT ~ 2e-2`, i.e. the appendix `zeta_B -> zeta_TTT` behaviour is preserved
exactly. PASS.

### Gate 4 -- SHEAR CLOSURE  (the decisive end-to-end gate)
KEY STRUCTURAL FINDING: the deployed `zeta_D` is built as LOW + HIGH, where
- LOW = `compute_kappa3_mod_zeta_equal_time` -> uses the (now-fixed)
  `SpinAwareZetaEvaluator` angular synthesis.
- HIGH = `compute_kappa3_mod_sigma3_high` -> uses `_kappa3_limber_alpha_kernel`,
  an independent flat-sky orientation integral that DOES NOT use the buggy
  evaluator (and was NOT changed by this fix).

At the worst config (gamma=10', phi=120deg, channel (2,2,-2)), deployed
ell_max=60 (LOW) + ell_cut=60/ell_high_max=1000 (HIGH):

| piece | sum\|.\| pre-fix | sum\|.\| post-fix |
|-------|------------------|-------------------|
| LOW   | 2.94e-21         | 3.19e-21 (+8.5%)  |
| HIGH  | 7.08e-13         | 7.08e-13 (unchanged) |
| total | 7.08e-13         | 7.08e-13 (change ~3e-22, ~1 part in 1e9) |

Two compounding reasons the fix is negligible for the deployed shear:

1. HIGH dominates LOW by ~8 orders of magnitude (the LOW spin-2 channels are
   `d^l_{2,-2}~gamma^4`-suppressed). The shipped canoes shear at this config is
   `|Gamma2|=|Gamma3|=1.56e-6 = |fold(zeta_D)|`, driven ENTIRELY by HIGH, which
   uses the bug-free `_kappa3_limber_alpha_kernel` (not the fixed evaluator).

2. At the deployed ell_max=60 the group-zeroing bug only kills 4156/58156
   cells = **7.1%** of triples carrying **7.2%** of the kernel |K| (vs ~half at
   ell_max=6 in the oracle). The masked-prefix problem only bites `(l2,l3)`
   groups whose smallest valid l1 is < |s1|=2, which is a small fraction at high
   ell_max. So the LOW change is only ~7-8%, not ~10^4x.

Net: the LOW-path fix moves the total `zeta_D` (hence the shear) by ~1 part in
1e9. The cosmic-shear discrepancy is UNCHANGED by this fix.

Three-way frame-invariant `can/ref` from a FRESH FULL REBUILD with the fixed
code (`stage1_ours_v2_FIXED.py` -> `stage1_ours_v2_FIXED.npz`;
`stage1_threeway_compare_FIXED.py` -> `outputs/stage1_threeway_compare_FIXED.npz`):

| config        | can/ref BEFORE | can/ref AFTER | frame-inv can_FIXED/can_old |
|---------------|----------------|---------------|-----------------------------|
| (10', 60deg)  | 4.118          | 4.118         | 1.0000 |
| (10', 30deg)  | 6.074          | 6.074         | 1.0000 |
| (10', 120deg) | 17.230         | 17.230        | 1.0000 |
| (50', 60deg)  | 9.473          | 9.473         | 1.0000 |
| median        | 7.773          | 7.773         | 1.0000 |

ref/fnc self-check median = 0.989 (reference still trusted).

VERDICT (Gate 4): the shear does NOT close. `can/ref` is UNCHANGED (bit-identical
frame-invariants, can_FIXED/can_old = 1.0000 at every config). The 4-17x
over-largeness lives in the HIGH/Limber spin-2 projection
(`_kappa3_limber_alpha_kernel`), NOT in the fixed LOW evaluator. This is the
honest, expected outcome given the LOW/HIGH magnitude hierarchy above.

### Gate 5 -- EXTREMES (boundary-validation re-probe)
Re-probed on the full grid (gamma in {10,50,150}', phi in {10..120}deg),
including phi=120 (worst, 17x), squeezed phi=10, and gamma=150' (`can_FIXED` vs
`can_old` frame-invariant, vs fastnc):

| config (extreme)        | frame-inv FIX/old | FIX/fnc | finite |
|-------------------------|-------------------|---------|--------|
| (10', 120) WORST        | 1.00000           | 17.36   | yes    |
| (10', 10) squeezed      | 1.00000           | 3.06    | yes    |
| (150', 120) wide        | 1.00000           | 9.61    | yes    |
| (150', 60)              | 1.00000           | 5.46    | yes    |
| (50', 90) mid           | 1.00000           | 8.52    | yes    |

All FIXED Gamma finite; max|Gamma| = 1.56e-6 (NO new blow-up). Every
frame-invariant is unchanged to 5 decimals.

Caveat (investigated, benign): at two SQUEEZED cells -- (50',10) and (150',90) --
an individual component ratio FIX/old reaches ~407-456x. This is NOT a blow-up: it
is a per-component LABEL PERMUTATION among the near-degenerate D-channels. E.g. at
(50',10) Gamma1 (2.0e-9 -> 8.3e-7) and Gamma3 (8.3e-7 -> 2.0e-9) simply SWAP; the
frame-invariant is identical (FIX/old = 1.0000). The tiny LOW change (the
gamma^4-suppressed channels are mutually comparable at squeezed configs) tips
which component holds which value, with zero effect on the physical
(frame-invariant) magnitude or the closure. PASS.

## HONEST VERDICT
The fix is a real, minimal, oracle-verified correction of a genuine bug (the
group-batched reduction dropped ~half of every `s1!=0` angular kernel, making the
LOW `zeta_D`/`zeta_PPP` too small). It does NOT close the cosmic-shear 3PCF
discrepancy, because the shear `zeta_D` is dominated by the HIGH (Limber)
`_kappa3_limber_alpha_kernel` path, which is independent of the fixed evaluator.
The 4-17x over-largeness vs the trusted `B_kappa`-phase reference therefore lives
in the HIGH/Limber spin-2 projection (or in a genuine great-circle-vs-flat-sky
projection-convention difference), NOT in the negative-spin LOW kernel. Per the
boundary-validation rule, no further fix was stacked to force closure.
