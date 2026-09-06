# The kappa3 vertex rebuilt at n_phi = 512, with the pair-phase fix

**Date.** 2026-09-06.
**Status.** Products built and checked. **Nothing in the manuscript has been changed,
and the deployed table has not been replaced.** This document is the record so the
paper side can decide whether to adopt.

**What it answers.** The two open defects recorded in
`project_deployed_vertex_defects_2026-09-06.md` are both closed. The handoff that
canoes wrote for this job on 2026-09-02
(`canoes/docs/handoff/2026_09_02_parent_handoff_nphi_rebuild.md`, never delivered)
was followed as written; its predicted acceptance behaviour reproduced.

---

## 1. What was wrong, and what was done about it

| defect | fix |
|---|---|
| `zeta_Bmod` computed with the residual pair phase dropped. `_kappa3_limber_alpha_kernel` took a fast path `if spin_sum == 0: return 2*pi*J_0`, correct for (0,0,0) but wrong at the modulus channels' (0,2,-2), where S = 0 is reached by cancellation and the residual phase is exp(-2i phi). | Already fixed upstream in canoes `a95c0d6` (2026-09-02). The deployed table predated it. Picked up for free by rebuilding against current canoes. |
| `n_phi = 64` under-integrates the FK placement. In slot 1, `(c, 1, c)`, the phi integrand oscillates as exp(i ell theta cos phi), so at ell_max = 15360 a 12 arcmin separation carries about 17 oscillations against 64 Gauss nodes. | Rebuilt with `--n-phi 512` on every band job. |

A third, opportunistic change went in at the same time, as the canoes handoff
suggested: the LOW piece was built at `--Nmax 384` rather than the default 128,
which is what canoes decision record E2 assumes production callers use.

## 2. Acceptance test

`zeta_TTT` has spins (0,0,0) and must be invariant under relabelling the legs, so
the slot 0 `(1, c, c)` and slot 1 `(c, 1, c)` rows at the same (gamma, shell) must
agree. The residual is a calibrated error bar with known true value zero, and the
test needs no external reference.

```
canoes/.venv/bin/python canoes/scripts/check_vertex_leg_symmetry.py TABLE.npz
```

| table | result |
|---|---|
| deployed `equal_time_limber_cut15360_permaware/table_permclosed.npz` | 864 of 864 cells outside tolerance, **FAIL** |
| rebuilt `pieces_nphi512/table_permclosed_np512.npz` | 0 of 864 outside tolerance, **PASS** |

Median |slot0/slot1 - 1| over 171 collapsed placements: 1.652e-1 deployed,
1.512e-2 rebuilt, an 11x improvement. Run independently twice, by the session that
built the table and by a separate verification agent.

This is an internal consistency check against a known true value, not a comparison
with an external ground truth. It says the quadrature now integrates what it claims
to integrate. It does not by itself establish that the vertex is right.

## 3. What moved

All six channels moved. **The note's statement that the phase fix leaves TTT, TTP,
TPP, PPP and Dmod bit-identical is about the phase fix in isolation; it does not
describe this rebuild**, where the n_phi change touches everything.

Median |new/old| per channel: TTT 0.9922, TTP 0.6653, TPP 0.6373, PPP 0.7118,
Bmod 0.4112, Dmod 0.7191. Bmod/TTT median falls from 0.98 to 0.35, which is the
pair-phase fix arriving: the modulus channel was the scalar channel.

### FK against Order-0, z_s = 5

The two defects push in opposite directions. The net is downward, and it is
scale dependent: the suppression is strongest at small separation and reverses
sign near 30 arcmin, where the n_phi gain overtakes the phase loss.

| gamma | deployed | rebuilt | new/old |
|---|---|---|---|
| 0.5' | 2.307 % | 1.742 % | 0.755 |
| 2' | 1.737 % | 1.339 % | 0.771 |
| 12' | 1.307 % | 1.148 % | 0.878 |
| 35' | 0.970 % | 1.058 % | 1.091 |

Over the 2' to 12' band: 1.33 to 1.73 % becomes 1.15 to 1.34 %.
FK/FF at 2.06': 6.94 becomes 5.35. At 10.77': 3.53 becomes 3.06.

**Callable caveat, important for like-for-like quoting.** The deployed FK sweep was
folded with the perm-aware callable; `run_fk_variant.py` defaults to the production
callable. Holding the callable fixed, the table-only change at z_s = 5 and 0.5' is
2.307 % to 1.761 % (perm-aware both sides) or 2.295 % to 1.742 % (production
callable both sides). Quote one of those pairs, not a mixture. The callable effect
on the harmonic B/E statistic is null (+0.0002 in band median); it is not null for
the real-space |xi_- / xi_+| at 0.5', where it is a factor 1.44.

### Cutoff ladder, tree arm, FK/Order-0 at 0.5'

| ell_max | deployed | rebuilt |
|---|---|---|
| 960 | 0.483 % | 0.339 % |
| 1920 | 0.967 % | 0.695 % |
| 3840 | 1.531 % | 1.126 % |
| 7680 | 2.023 % | 1.514 % |
| 15360 | 2.310 % | 1.742 % |
| geometric extrapolation | 2.710 % | 2.064 % |

The convergence behaviour is essentially unchanged and marginally slower. The final
doubling still adds 14 to 15 % at 0.5' (14.17 % deployed, 15.04 % rebuilt), and the
extrapolation overshoot grows from 17.4 % to 18.5 %. **The statement that the quoted
amplitude is a floor still climbing survives the rebuild.** What changed is the
level, not the story.

The extrapolation recipe is the one in `reproduce/check_numbers.py` lines 126-129,
unchanged: `inc = diff(v); r = inc[-1]/inc[-2]; extrap = v[-1] + inc[-1]*r/(1-r)`.
Fitting the last three increment ratios instead of the last two gives 2.349 %, so
2.064 % should not be read as carrying three significant figures.

### Multi-redshift, FK/Order-0 at 0.5'

| z_s | 1 | 1.7 | 2.5 | 3.2 | 4 | 5 |
|---|---|---|---|---|---|---|
| deployed | 1.364 % | 1.772 % | 2.037 % | 2.138 % | 2.245 % | 2.307 % |
| rebuilt | 1.056 % | 1.364 % | 1.562 % | 1.636 % | 1.716 % | 1.761 % |

The suppression factor is essentially independent of redshift (0.7739 down to
0.7634), so the monotone growth with source redshift is unchanged in shape.

t_final grid used, taken from the deployed run record rather than reconstructed:
1822.7212585, 2094.89417743, 2216.79339499, 2266.58182201, 2297.28925685,
2313.0288751857356.

## 4. The B-mode result, and a correction to the project note

**The note's prediction is wrong, and the error is specific.** It states that fixing
the phase makes the FK B-mode to E-mode ratio rise by up to 2x. The harmonic ratio
falls.

| quantity | deployed | rebuilt |
|---|---|---|
| C_BB / C_EE, band median over 50 <= ell <= 1500 | 0.6287 | 0.2889 |
| C_BB / C_EE at ell = 60 | 0.9511 | 0.6784 |
| C_BB / C_EE at ell = 1500 | 0.4203 | 0.0266 |
| real-space \|xi_- / xi_+\| at 0.5' | 8.576e-4 | 2.390e-3 |

The note's *measurement* is correct and reproduces to four digits: |xi_- / xi_+|
does rise. What is wrong is the identification of that ratio as "the B/E ratio",
and the direction of the inference drawn from it. With
EE = (Cl[xi_+] + Cl[xi_-])/2 and BB = (Cl[xi_+] - Cl[xi_-])/2,

    B/E = (1 - q) / (1 + q),   q = Cl[xi_-] / Cl[xi_+]

which is strictly decreasing in q. xi_- = 0 is the *maximal* B-mode limit, BB = EE
exactly, not the minimal one. A rising |xi_- / xi_+| therefore implies a falling
B/E. Checked as a counterfactual on the note's own premise: scaling the deployed
xi_+ with xi_- held fixed, which is the phase-fix-alone scenario, gives band-median
B/E of 0.629, 0.596, 0.557, 0.529, 0.450, 0.374 while |xi_- / xi_+| rises. Even
restricted to the phase fix alone the note's conclusion is backwards.

The note's remaining claim, that the B-mode does not go to zero, stands: B/E is
still 0.68 at ell = 60 and 0.25 to 0.47 through ell of about 400. It collapses only
above ell of about 1000.

### What can and cannot be quoted about the B-mode

Robust:
* FK remains the dominant Order-2 B-mode source. Its band-median B/E of 0.289 is
  about 9x FF's 0.032, and FF is unchanged by the rebuild because it contains no
  n-point cumulant vertex.
* Delta C_ell^EB is exactly 0.000e+00 on every fold, deployed and rebuilt, as is the
  raw real-space <g+ gx>. This is the parity theorem, not a numerical result.
* C_EE > 0 and C_BB > 0 at all 71 dense multipoles in the band, on all four folds.
* The fall with ell, and its rough size.

Not robust:
* **The value at high ell.** C_BB there is a small difference of two nearly equal
  transforms; the cancellation factor (|Cl[xi_+]| + |Cl[xi_-]|) / (2 |C_BB|) reaches
  40 at ell = 1500 on the rebuilt fold. Propagating the table's own median
  leg-symmetry residual of 1.5 % gives sigma(B/E) = 0.015, constant across the band.
  So B/E at ell = 1500 is 0.027 +/- 0.015, an error bar of order the value.
* **The exact band median.** 0.289 is on the paper's 16 log-spaced multipoles. On
  all 1451 integer multipoles the pair is (0.520, 0.138) and the ratio is 0.27
  rather than 0.46. State the grid with the number.
* **Smoothness.** C_BB is not a power law in ell in either table
  (d ln|C_BB| / d ln ell swings from -4.45 to +1.22 deployed, -8.13 to +3.54
  rebuilt). This is pre-existing and unchanged in character, so it is not evidence
  against the rebuild, but smoothness cannot be used as an argument for physicality.

## 5. How to reproduce

```bash
CAN=/Users/zzhang/projects/angular_statistics/canoes
cd sachs_sft/callables/kappa3_vertex/rebuild
export PYTHONPATH=$CAN/src CANOES_SUPPRESS_METAL_WARNING=1

# 1. the 17 build jobs: 8 tree bands on r4, 8 extra bands, 1 LOW at Nmax 384.
#    ~4.5 h wall at 8 concurrent, 8 of 28 cores, about 23 GB peak.
python run_queue.py --jobs-json jobs_np512_stage1.json --max-jobs 8 --min-free-gb 30

# 2. merge r4 (1828 rows) with extra (314) into the perm-closed 2142-row bands
for w in 00060_00120 00120_00240 00240_00480 00480_00960 \
         00960_01920 01920_03840 03840_07680 07680_15360; do
  $CAN/.venv/bin/python merge_bands.py \
    --parts products/pieces_nphi512/band_tree_r4_${w}_np512.npz \
            products/pieces_nphi512/band_tree_extra_${w}_np512.npz \
    --triples products/triples_permclosed_r4.npz \
    --out products/pieces_nphi512/band_tree_permclosed_${w}_np512.npz
done

# 3. assemble, then check
$CAN/.venv/bin/python assemble.py \
  --low products/pieces_nphi512/low_permclosed_r4_np512.npz \
  --bands 'products/pieces_nphi512/band_tree_permclosed_*_np512.npz' \
  --cutoff 15360 --out products/pieces_nphi512/table_permclosed_np512.npz
$CAN/.venv/bin/python $CAN/scripts/check_vertex_leg_symmetry.py \
  products/pieces_nphi512/table_permclosed_np512.npz

# 4. fold (97 s). Add --callable for a perm-aware fold; see the caveat in section 3.
<sft-wick python> run_fk_variant.py products/pieces_nphi512/table_permclosed_np512.npz --n-jobs 8
```

The build is single-threaded per job at about 2.8 GB, so concurrency is limited by
cores rather than memory on a 28-core machine. `--chunk-rows 24` at n_phi = 512
holds the (rows, n_ell, n_ell, n_phi) working array at the same size that
`--chunk-rows 192` held it at n_phi = 64.

**Serialise the folds.** The project rule that no two folds run at once was violated
during this work by concurrent agent sessions. Nothing downstream appears to have
been corrupted, but gate on
`pgrep -f 'python.*run_FF_single|python.*run_fk_variant'`. A bare `run_FF_single`
pattern also matches waiting shells and never clears.

## 6. Products

Tracked in this repo:

| file | what |
|---|---|
| `pieces_nphi512/table_permclosed_np512.npz` | the rebuilt vertex, 2142 rows, 16 shells, cut 15360 |
| `pieces_nphi512/low_permclosed_r4_np512.npz` | the LOW piece at Nmax 384 |
| `pieces_nphi512/table_permclosed_np512_xi.npz` | FK fold, production callable |
| `pieces_nphi512/table_permclosed_np512_permaware_xi.npz` | FK fold, perm-aware callable |
| `pieces_nphi512/table_deployed_plaincall_xi.npz` | deployed table, production callable, the 2x2 control |
| `pieces_nphi512/multiz_fk_np512_6planes.npz` | the six-plane FK series |
| `pieces_nphi512/ladder/table_tree_np512_cut*_xi.npz` | the four lower ladder rungs |
| `verify_bmode_np512.py` | the independent B-mode verification |
| `run_fk_variant_multiz_np512.py`, `report_multiz_np512.py` | the multi-z driver and report |
| `pieces_nphi512/ladder/measure_ladder_np512.py` | the ladder measurement |
| `jobs_np512_stage1.json` | the recipe above |

Excluded, following the same convention as `products/pieces/`: the 24 octave-band
intermediates and the 5 ladder tables, about 26 MB of regenerable canoes output.
Every one is recorded with its md5 in
`reorg_2026-09/nphi512_manifest_2026-09-06.txt`.

## 7. What was NOT done

* **The manuscript is untouched.** No .tex file, no figure in `figures/`, no entry in
  `reproduce/numbers.json`. Adopting this rebuild means re-running
  `check_numbers.py` and `regen_figure.py --all` and rewriting the affected claims.
* **The deployed table is still deployed.** `equal_time_limber_cut15360_permaware/`
  is untouched, so every existing figure still reproduces byte-for-byte.
* **BiHalofit was not rebuilt.** The nine jobs sit unbuilt in
  `jobs_np512_stage2.json`, about 4 h. Until they are run the BiHalofit ladder is
  still the n_phi = 64 product and **must not be plotted on the same axes as the new
  tree curve**.
* **The exact-versus-Limber question (D-CN-05) is untouched.** D-CN-06 bounds it:
  inside 0.5 to 12 arcmin the exact branch is at most 0.44 % of the vertex.
* **The spin-2 zeta_D over-normalisation is untouched.** Diagnosed on the
  uncommitted branch `fix/spin2-zetaD-2_2_-2-norm`; the response derivation has not
  been done.
* **The FK callable still sorts the query triple and symmetrises the coupling
  tensor**, which is wrong for the spin channels, and its 8-decimal deduplication
  still cannot represent separations below gamma = 0.344 arcmin. See
  `project_fk_vertex_sampling_artifact_2026-08-26.md`.

## 8. One pre-existing inconsistency found in passing

In the **deployed** products the cutoff-ladder top rung and the FK sweep agree at
0.5' (2.310 % against 2.307 %) but not across the band: the ladder's 15360 rung
gives 1.836 / 1.753 / 1.656 / 1.596 % at 2 / 5 / 8 / 12', where the sweep gives the
1.33 to 1.73 % the paper quotes. They are different products. The ladder arm is
built on the 1828-row r4 triple set with LOW Nmax = 128, the sweep on the 2142-row
perm-closed set. This is not caused by the rebuild, and the new ladder does not have
it: its rungs and the rebuilt table share a row set, and its top rung is literally
the same fold.

## 9. Where the talk stands

The v2 conference deck at `~/Documents/MyTalks/SFT_WL_talk/v2/` (branch
`talk/v2-sft-wick-2026`) has been updated to the rebuilt numbers throughout, and its
four lensing figures are regenerated from these products. Mixing the rebuilt FK with
the existing Order-0 and FF sweeps is correct there, since neither contains an
n-point cumulant vertex. If the paper adopts the rebuild the two will agree; if it
does not, the deck and the manuscript will quote different numbers, which is worth
knowing before the talk is given.
