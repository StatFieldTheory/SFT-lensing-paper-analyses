# Corrected kappa3 vertex callable

> **Table archived on 2026-10-02.** `table_permclosed.npz`, the August table (n_phi = 64, pair phase
> dropped) that `perm_aware_kappa3_callable.py` opens by default, was moved to
> `_archive/cleanup_2026-10-02/` (index: `reproduce/archived_inputs.json`). The callable is unchanged and
> is still the production callable; set `TABLE_PATH` before the first call. The manuscript's table is
> `../rebuild/products/pieces_nphi512/table_permclosed_np512.npz`.

Fixes the two defects the note reports in the production callable,
`sachs_sft/callables/kappa3_vertex/equal_time_limber/equal_time_limber_kappa3_callable.py`.
Nothing here overwrites the production file; the replacement is exercised
through the production runner with `run_fk_variant.py --callable`.

## Defect 1: leg permutations

The vertex is a function of three (field, direction) pairs, and the
tabulated channels fix which leg carries which field: `zeta_TTP` has spin
weights `(0, 0, 2)`, so the spin-2 leg is the third, and `zeta_TPP` has
`(0, 2, 2)`, so the scalar leg is the first.

The production callable sorts the query triple and then writes one value
into every index permutation of the coupling tensor. Asking for the spin
field on a different leg is a question about a **permuted geometry**, not
about the same geometry with a relabelled channel, so the symmetrised fill
is exact only for `zeta_TTT` and for geometries that happen to be
permutation symmetric.

Measured cost on the deployed table, over asymmetric rows: the corrected
`K001` and `K010` differ by a median **115%** of their magnitude, where the
production callable forces them equal. At the collapsed configurations the
channel-level ambiguity is a median 23% on `zeta_TPP` and 70% on
`zeta_Dmod`.

The correction evaluates each entry at the geometry its own field
assignment requires:

| entry | rule |
|---|---|
| no spin leg | `zeta_TTT`, any ordering |
| one spin leg at `p` | permute `p` into slot 2, read `zeta_TTP` |
| two spin legs, scalar at `q` | permute `q` into slot 0, read `(zeta_TPP, zeta_Bmod)` |
| three spin legs, all `Re` | `zeta_PPP` and `zeta_Dmod`, any ordering |
| three spin legs, one `Re` two `Im` | three cyclic orderings, combined |

The last row is the case the table cannot answer with one lookup.
`(zeta_Dmod - zeta_PPP)/4` is the **average** of two of the three placements
of the `Re` leg, so three cyclically permuted evaluations are combined as
`mean(A) - mean(B) + mean(C)` to isolate the wanted placement. On a
symmetric geometry the three are equal and the combination collapses to the
production value, which is one of the tests.

Because the corrected lookup asks for permuted geometries, the table must
contain them. The production mesh already does, being a symmetric meshgrid;
the densified collapsed family is made so by
`code/rebuild/make_perm_closed_triples.py`. The callable checks this at
load and warns if the table is not closed, rather than letting a missing
row be answered silently by a different configuration.

## Defect 2: the mesh key

The production callable deduplicates triples after rounding to 8 decimals,
so two collapsed configurations closer than `1 - cos g = 5e-9` become one
node however finely the table samples. That is a hard floor at
**g = 0.344 arcmin**, and the draft plots down to 0.5.

The key here is 12 decimals, a floor at 0.0034 arcmin. The production mesh
deduplicates to the same 351 rows at either precision, so nothing else
changes.

## Cost

About 1.9 times the production callable. The naive implementation was 21
times, because it queried the tree once per (channel, ordering) pair; the
tree query depends only on the ordering, so it is done once per ordering
and the six channels are read off it.

## Files

| file | what it is |
|---|---|
| `perm_aware_kappa3_callable.py` | the replacement callable |
| `tests/test_perm_aware.py` | 7 tests, including one that fails if the fix becomes a no-op |

## Running

```sh
CAN=/Users/zzhang/projects/angular_statistics/canoes
PYTHONPATH=$CAN/src $CAN/.venv/bin/python -m pytest tests -q

SFT=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
$SFT ../run_fk_variant.py <permutation-closed table>.npz \
     --callable perm_aware_kappa3_callable.py --tag permfix --n-jobs 6
```
