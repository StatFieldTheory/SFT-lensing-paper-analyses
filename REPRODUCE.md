# Reproducing the SFT weak-lensing results

**AI-generated documentation:** this guide was written with AI assistance and
checked against the repository's entry points and saved records. The recorded
configurations, source and input hashes determine the calculation.

There are two levels of reproduction: regenerating figures and numerical
summaries from saved products, and recomputing those products from the driving
fields. The first was revalidated in an isolated package copy on 2026-09-30. The second requires the
external runtimes, potentially hours of computation, and care when relocating
recorded paths. A clean-machine installation has not been validated.

## 1. Current results and inputs

[reproduce/active_products.json](reproduce/active_products.json) selects the R1
revision results. The main source is `z_s = 5`, with
`lambda_s = 2317.9695958203943 Mpc`; the main FK calculation uses
`ell_max = 15360`, corrected pair phases, `n_phi = 512`, and the
permutation-aware three-point-vertex callable. Lower source planes use their
corresponding affine endpoints.

The manifest contains 12 selections and 22 hashed file members, including the
cutoff-ladder directory. Readers check the selected products before use. A
missing or modified member causes an error. The manifest is required: without
it the readers stop with an error, because the historical inputs they once
fell back to were archived on 2026-10-02. A candidate manifest is not an
activation record.

[reproduce/archived_inputs.json](reproduce/archived_inputs.json) lists the 53
archived files with their SHA256 hashes and the commit that still holds them.
The copies live under the local, untracked `_archive/cleanup_2026-10-02/`. The
referee-revision run records pin 31 of them by path and hash, so the
recomputation drivers read them through `reproduce/archived_inputs.py`, which
checks each archived copy against its recorded hash and keeps the recorded
path. On a machine without the archive those drivers stop and name the commit
to recover from. Figure regeneration and the number check do not need it.

Useful records:

- [Cleanup validation](reproduce/provenance/cleanup_20260930_validation.json):
  17 figures regenerated without local notes or archived analyses, matching
  figure content and 33 unchanged numerical summaries. Existing dependency
  environments were reused.
- [Archive validation](reproduce/provenance/archive_20261002_validation.json):
  the same 17 figures and the 32 remaining numerical summaries after the
  superseded inputs were archived, with the driver checks that were run.
- [Framework provenance](sachs_sft/analyses/r1_sft061/framework_provenance.json):
  external commits, working-state hashes and selected input hashes.
- [Source map](sachs_sft/analyses/r1_sft061/source_map.json): source redshifts,
  affine endpoints and response convention.
- [Closeout validation](sachs_sft/analyses/r1_sft061/closeout_validation.json) and
  [figure comparison](sachs_sft/analyses/r1_sft061/closeout_figure_check.txt):
  17 regenerated figures with identical content and zero raster differences;
  the number checker completed 33 entries.
- [Order-0 convergence checks](sachs_sft/analyses/analysis1/r1_aligned/source_folded_final_validation.json)
  and [cutoff assembly](sachs_sft/analyses/r1_sft061/figure_products/cutoff_ladder/assembly_manifest.json):
  checks specific to the accepted calculations.

## 2. Environments and layout

Run commands from this repository's root. Figure regeneration launches the
interpreters selected by `STF_PYCCL` and `STF_SFTWICK`; set these and
`CANOES_ROOT` to your installations. The original workstation defaults are:

```bash
export STF_PYCCL=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
export STF_SFTWICK=/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
export CANOES_ROOT="$HOME/projects/angular_statistics/canoes"
```

These are examples of the validated layout, not portable installation commands.
The PyCCL environment supplies NumPy, SciPy, Matplotlib, PyCCL, CAMB and Pillow;
numerical run preparation also uses PyYAML. PDF comparison requires Poppler's
`pdftoppm` executable. SFT-Wick supplies the expansion engine and diagram
rendering. Some generators and checks import canoes or read its power-spectrum
tables even when using saved products.

The revision runtime pins are:

| Dependency | Recorded version or state |
|---|---|
| SFT-Wick | 0.6.1, commit `6945815ad9e6efa2b9e768891d6bf92b1214670d`, clean tracked source |
| canoes | Commit `077a64c86788add650ece263b054e7803ca0c644` plus the [recorded patch](sachs_sft/analyses/r1_sft061/canoes_worktree.patch) |
| PyCCL | 3.3.1 in the aligned comparison |

Apply the canoes patch to the recorded commit in a separate checkout if
reconstructing the numerical runtime. Its SHA256 is recorded in
`framework_provenance.json`. The package lists in
`reorg_2026-09/env/pip_freeze_PyCCL.txt` and
`reorg_2026-09/env/pip_freeze_sft-wick.txt` are historical inventories, not
complete replacement specifications for the later revision.

**Portability:** interpreter overrides do not relocate every path. In particular,
several numerical scripts and saved configurations contain absolute paths to
inputs or external repositories. The reproduction entry points resolve the
package location, but their numerical dependencies may still require rebasing. Rebase these paths before running elsewhere. Preserve the
original records and create new run manifests; the runtime/input guards reject
unreviewed changes. Saved `*.source.py` files are immutable run snapshots.

The author's layout is:

```text
STF_lensing/                       separate manuscript project
  main.tex, sections/, figures/
  SFT-lensing-paper-analyses/       this repository
```

The manuscript source and deployed figure PDFs are not supplied by this analysis
repository. `check_figures.py` expects those PDFs in `../figures/`.

## 3. Regenerate figures and numerical summaries

With the environments and paths configured:

```bash
"$STF_PYCCL" reproduce/regen_figure.py --list
"$STF_PYCCL" reproduce/regen_figure.py --all
"$STF_PYCCL" reproduce/check_figures.py
"$STF_PYCCL" reproduce/check_numbers.py
```

The figure command runs generators sequentially, reusing the accepted arrays and
Monte Carlo cache. It writes each generator's output directory and execution
logs. It does not run new vertex builds or SFT folds. The number command prints
32 summaries and writes `reproduce/numbers.json`; it is a claim-versus-product
check, not an automatic comparison with the current manuscript text or a proof
of numerical convergence.

To rebuild and compare selected figures:

```bash
"$STF_PYCCL" reproduce/regen_figure.py --figure 1 2 3
"$STF_PYCCL" reproduce/check_figures.py --figure 1 2 3
```

The integers are reproduction IDs, not manuscript figure numbers. IDs 8/9 and
15/16 share generators. `--list` prints the active paths, arguments and interpreter choices
selected by the current code. In particular, ID 1 uses the aligned
`analysis1/r1_aligned/plot_comparison.py` under the active manifest.

`IDENTICAL(raw)` means byte equality; `IDENTICAL(content)` permits differences
in PDF date metadata. `EQUIVALENT(bbox drift)` indicates equal raster content at
the check's resolution despite a bounding-box difference. `DIFFERS` or `MISSING`
requires investigation. Regeneration does not replace `../figures/`;
`reproduce/deploy.py` is the separate deployment step and updates
`reproduce/deployed_md5.txt`.

## 4. Recompute the numerical products

Full recomputation is a collection of guarded workflows, not a single turnkey
command. Inspect the accepted manifests and `--help` for each driver, use fresh
output paths and validate each replacement before changing the active manifest.
The accepted products remain available for comparison.

| Calculation | Drivers and accepted records |
|---|---|
| Order-0 and PyCCL reference | `sachs_sft/analyses/analysis1/r1_aligned/source_folded_reference.py`, `ccl_reference.py`; saved `*.json` configurations and `source_folded_final_validation.json` |
| Main and lower-redshift FK | `sachs_sft/analyses/r1_sft061/replay_suite.py`; `true_redshift_fk_gl24_corrected_main/prepared.json` and `true_redshift_fk_gl24_corrected_multiz/prepared.json` |
| FF source planes | `sachs_sft/analyses/r1_sft061/companion_suite.py`; `true_redshift_ff_gl24_all/prepared.json` |
| Permutation-aware cutoff ladder | `sachs_sft/analyses/r1_sft061/permaware_ladder/pieces.py`, `suite.py`, `guard.py`; `contract.json` and `folds_{tree,bihalofit}/prepared.json` |
| Monte Carlo controls | `sachs_sft/analyses/r1_sft061/mc_update/` drivers and per-run manifests |
| Figure-product assembly | `sachs_sft/analyses/r1_sft061/assemble_figure_products.py`, `prepare_plot_inputs.py` and saved assembly manifests |

For example, the FK replay driver separates `plan`, `prepare` and `run`.
Preparing a fresh main calculation requires `--source-mode true-redshift`,
`--input-mode corrected-nphi512`, the reviewed `--source-map`, and
`--targets main`. The defaults include historical input modes, so do not use
an unqualified command as a reproduction of the current paper. The corrected
cutoff ladder has its own permutation-aware workflow; the older restricted
sorting-callable batches are retained as rejected controls.

Run SFT folds serially and avoid overlapping vertex builds. Worker count affects
memory consumption substantially. FF sweeps, vertex integration and Monte Carlo
controls can take hours. Some regeneration-only intermediates, including
historical octave-band pieces, are excluded from Git; reconstruct those inputs
with the corresponding build workflow when needed. Accepted assembled products
and their checks remain tracked.

Validation has limits. The FF outer-quadrature probes do not establish convergence
of its underlying 21-source covariance table. The FF Monte Carlo control uses a
matched local covariance model. Recovered FK BB in this scalar demonstration is
a numerical residual; the pure-E transform control is input-specific. Finite
cutoff increments do not establish an infinite-cutoff limit. Keep these limits
when comparing new products with the recorded results.

## 5. Build the separate manuscript

If the manuscript checkout is available one level above this repository, run
from that manuscript directory:

```bash
cd ..
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

This compiles the manuscript's own sources and deployed figures. It does not
build a document from the analysis repository. Review the current compiler
output rather than assuming that an earlier manuscript build still applies.

Local working notes, rollback snapshots, generated previews, caches and logs
are excluded where specified in `.gitignore`. They are not prerequisites for
redrawing the figures from the accepted products. The scientific inputs and
machine-readable provenance needed for that path are retained in the package.
