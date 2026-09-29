# Corrected n_phi=512 inputs and minimal ladder rebuild

## Status

The author approved the corrected kappa3 inputs. No vertex build or framework
fold has been launched by this preparation. Historical artifacts and the
controlled old-input framework replay are preserved.

The 31 existing corrected files match the September 6 checksum manifest.
`corrected_input_inventory.json` records their hashes, geometry and timing
evidence. The main permutation-closed table has 2,142 rows, 16 affine shells,
HIGH quadrature `(n_ell, n_phi) = (96, 512)`, LOW `Nmax=384`, and
`ell_max=15360`. Its permutation check passed all 864 tested comparisons.
Corrected tree inputs exist at all five ladder cutoffs. Corrected BiHalofit
inputs still require a build.

## Prepared main run

After the root-controlled FF probe finishes, the prepared main batch can run
with four workers. It retains all 40 angles, all six requested observable
pairs and three internal field components. Its endpoint is the aligned true
source redshift 5, `lambda=2317.9695958203943` Mpc.

Run from the reproduction-package root:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -u \
  sachs_sft/analyses/r1_sft061/replay_suite.py run \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24_corrected_main
```

The batch uses the corrected main table and the existing permutation-aware
callable. `replay_manifest_corrected_main.json` records its inputs. Main and
multi-redshift products continue to require 40 angles.

## Exact restricted geometry for the cutoff ladder

Figure 17 needs the convergence result at 0.5 arcmin. The associated angular
claims require interpolation at 2, 5, 8 and 12 arcmin. Original angle indices
`[0,5,6,7,8,9,10,11,12,13,14]` retain 0.5 arcmin, the relevant interpolation
neighbors and every original grid point inside 2 to 12 arcmin.

The actual cached sft-wick 0.6.1 FK diagrams have no integrated free angular
direction. Their vertex legs are assigned to the two external directions.
`corrected_kappa3_ladder52/query_contract.json` records the expansion-cache
hash, both spatial maps and an enumeration of all eight three-leg placements
of the external directions. This enumeration is a superset of the diagrams'
collapsed configurations.

The historical ladder callable sorts the cosine triple, groups rows rounded
to eight decimals, averages all raw duplicate members of each group and
interpolates between four nearest groups. For the retained queries, only
44 canonical groups are needed. Preserving all their raw members gives
52 of the original 1,828 rows. Both the four neighbor coordinates and their
distances are bit-identical between the full and restricted geometries.

The actual callable was then compared using the corrected tree data on the
full historical r4 geometry and the restricted geometry. All 2,728 sampled
full 3 by 3 by 3 tensors were bit-identical, with maximum absolute difference
zero. The samples include 11 angles, eight leg placements, all 16 affine
shells and all 15 inter-shell midpoints. This verifies the table restriction
for the specified queries. It does not establish quadrature convergence.

The restriction retains all six vertex channels, all 16 radial shells,
`n_ell=96`, `n_phi=512`, and LOW `Nmax=384`. No missing component is filled
with zero. The restricted LOW file is an exact row selection from the
existing corrected LOW piece, with its fingerprint recomputed. No LOW
integration is needed.

Prepared files:

- `corrected_kappa3_ladder52/triples_ladder52.npz`
- `corrected_kappa3_ladder52/low_ladder52_np512.npz`
- `corrected_kappa3_ladder52/query_contract.json`
- `corrected_kappa3_ladder52/bihalofit_jobs.json`

## Serial BiHalofit build

The eight jobs cover the original disjoint HIGH windows from 60 to 15360.
Each uses 24-row chunks and full quadrature. Scaling the measured historical
serial BiHalofit work by the quadrature and row counts gives approximately
0.89 hours (54 minutes), plus startup and chunk overhead. This is an estimate,
to be replaced by measured times. The earlier 31-hour estimate applied to
all 1,828 rows and is not the proposed execution plan.

After root approval and only while no other vertex build is active:

```sh
mkdir -p sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/pieces \
  sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/logs
/Users/zzhang/projects/angular_statistics/canoes/.venv/bin/python -u \
  sachs_sft/callables/kappa3_vertex/rebuild/run_queue.py \
  --jobs-json sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/bihalofit_jobs.json \
  --max-jobs 1 --min-free-gb 30
```

The queue skips existing output paths without validating their contents.
Use fresh paths for the first run and validate existing products before any
resume. Its current free-memory check does not block the first job when no
job is active, so independently check available memory before dispatch.
The older `bihalofit_np512_jobs.json` describes full-r4 work and is retained
as provenance only. Do not launch that full-r4 queue.

## Assemble five cutoff tables for each model

Use cutoffs `960, 1920, 3840, 7680, 15360`, replacing `CUTOFF` below with one
of these values. Commands below must only target paths that do not exist.
`assemble.py` itself does not refuse overwrites. It validates matching LOW
and HIGH fingerprints, contiguous band coverage and a common model, and
retains all six channels.

BiHalofit, after all required HIGH bands finish:

```sh
/Users/zzhang/projects/angular_statistics/canoes/.venv/bin/python \
  sachs_sft/callables/kappa3_vertex/rebuild/assemble.py \
  --low sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/low_ladder52_np512.npz \
  --bands 'sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/pieces/band_bihalofit_r4_*_np512.npz' \
  --cutoff CUTOFF \
  --out sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/table_bihalofit_cutCUTOFF_np512.npz
```

Tree uses existing corrected permutation-closed pieces and selects the same
52 raw rows during assembly. No new tree quadrature is required:

```sh
/Users/zzhang/projects/angular_statistics/canoes/.venv/bin/python \
  sachs_sft/callables/kappa3_vertex/rebuild/assemble.py \
  --low sachs_sft/callables/kappa3_vertex/rebuild/products/pieces_nphi512/low_permclosed_r4_np512.npz \
  --bands 'sachs_sft/callables/kappa3_vertex/rebuild/products/pieces_nphi512/band_tree_permclosed_*_np512.npz' \
  --subset-npz sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/triples_ladder52.npz \
  --cutoff CUTOFF \
  --out sachs_sft/analyses/r1_sft061/corrected_kappa3_ladder52/table_tree_cutCUTOFF_np512.npz
```

## Implemented integration and prepared batches

`replay_suite.py --ladder-geometry contract52` explicitly selects the
restricted geometry with `--input-mode corrected-nphi512`. The default
historical-geometry option remains available. Preparation verifies every raw
member of the retained duplicate groups, bit-identical scalar and batch KD
neighborhoods, input hashes, all six finite channels, radial shells, cutoff,
model and quadrature metadata.

`restricted_kappa3_guard.py` wraps both actual callable entry points. It
requires exact identity with the external config directions before any
lookup and requires every non-x leg in one sample to use the same selected
y direction. An unknown direction, even a one-ULP perturbation, or a mixture
of different angular pairs fails. The wrapper imports and hashes the explicit
query contract and table, and uses the historical sorting implementation for
the numerical work. It does not approximate an unsupported query.

The five tree cutoff tables have been assembled without new quadrature.
The fresh batch `true_redshift_fk_gl24_corrected_tree52` is prepared with
four workers, 11 angles, three internal fields and only pair `(0,0)`.
`true_redshift_fk_gl24_corrected_multiz` is independently prepared with four
workers, five aligned source planes, 40 angles, the corrected permutation-
closed table and only pair `(0,0)`. The existing prepared main inputs and
generated files retain their recorded hashes.

`verify_restricted_ladder.py` exercised each of the five generated tree
wrappers against its full-r4 corrected counterpart at 2,728 scalar and
2,728 batched samples. All 27 tensor entries are bit-identical at every
sample. It also checked rejection of unsupported, mixed-angle and nonfinite
directions in both entry points and rejection of a missing modulus channel.
`corrected_kappa3_ladder52/guard_verification.json` records the PASS result.
No fold is part of that verification.

Run contract validation and restricted folds with the pinned `sft-wick`
interpreter shown above. Although the PyCCL environment reports the same NumPy
version, its scalar normalization produces a different exact query set. The
contract deliberately fails closed in that environment. The actual framework
environment passed both scalar and batch verification.

After code review, each prepared batch can be run separately through the
shared lock, substituting its directory below:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python -u \
  sachs_sft/analyses/r1_sft061/replay_suite.py run \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24_corrected_tree52
```

Once BiHalofit assembly finishes, prepare it in a fresh directory:

```sh
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
  sachs_sft/analyses/r1_sft061/replay_suite.py prepare \
  --source-mode true-redshift --input-mode corrected-nphi512 \
  --ladder-geometry contract52 \
  --source-map sachs_sft/analyses/r1_sft061/source_map.json \
  --targets bihalofit_960 bihalofit_1920 bihalofit_3840 bihalofit_7680 bihalofit_15360 \
  --n-jobs 4 \
  --work-dir sachs_sft/analyses/r1_sft061/true_redshift_fk_gl24_corrected_bihalofit52
```

Plot and number readers must consume the available `(0,0)` output directly.
The generic `load_observables` expects all six observable pairs and must not
be satisfied by invented zero arrays. The full source-folded main computation
is independent of this restricted-ladder integration.

## Figure 7 table slices

The corrected main table includes collapsed raw rows `(1,c,c)` for all
22 main-grid angles at or below 85 arcmin (the largest selected grid angle
is approximately 71.255 arcmin), with all 16 radial shells. The largest
nearest-row angular offset caused by rounded config directions is
`1.657566599533311e-7` arcmin. Use the table's actual angular coordinate as
the plotted abscissa, or record nearest-row selection with an explicit
`1e-6` arcmin tolerance. New vertex integrations are unnecessary for these
slices.
