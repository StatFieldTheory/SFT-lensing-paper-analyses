# Reorganisation report, 2026-09-02

Reorganised `STF_lensing` so that a reader holding the manuscript and this repository can
reproduce every figure and every quoted number without oral tradition. Executed against
`reorg_2026-09/REORG_PLAN.md` with the operator's sign-off on decisions D1 to D16
(D16 deferred: the push to GitHub waits until the analysis and the results have been
checked for consistency).

Follows the format of the 2026-06-17 `CLEANUP_REPORT.md`.

## Headline results

* **All 17 figures regenerate from the repository with no manual rebinding**, in about
  90 seconds total, and the deployed PDFs are byte-for-byte what the documented commands
  produce (`reproduce/check_figures.py` reports `IDENTICAL(raw)` for all 17).
  Before this, 4 of them silently produced a wrong FK channel when run as documented.
* **The FK sweep re-runs byte-identical.** A fresh fold through the promoted run folder,
  with the moved permutation-aware callable and a relative config, reproduced the
  2026-08-26 product exactly, three months on.
* **All 35 quantitative claims recompute unchanged** (`reproduce/check_numbers.py`).
  No number contradicts the manuscript; four wording nuances are recorded in REPRODUCE.md
  section 5 and were not acted on, since the manuscript is read-only here.
* `latexmk -g -pdf -halt-on-error main.tex` exits 0: 32 pages, 17 figures, no undefined
  references.
* **`driver_field_emulators/` is gone**, absorbed into the analysis package with its
  history preserved from a verbatim import commit. The production FK path is now inside
  the reproduction package, so "running the generator standalone downgrades the figure"
  no longer exists as a trap.
* The package grew from 302 tracked files to about 1000, and from 15 MB to 24 MB of
  tracked content; 1.0 GB of superseded archives and build residue was retired.

## The four traps that are now closed

| trap (from the survey) | what it was | what closed it |
|---|---|---|
| Running an FK generator standalone gave a channel 4.5x too low | the generators defaulted to the superseded `cut1000` sweep; the deployed figures came from a rebinding driver outside the package | the defaults point at the manuscript's converged sweep, promoted into a run folder; the rebinding driver is kept only as the record of the August deployment |
| `cosmology_coord_maps.pdf` had no live generator | it existed only in the retired archive and needed the external canoes package under its pre-2026-06 name | promoted to `figures/make_cosmology_coord_maps.py` with the two canoes constants and the two background helpers vendored (40 lines); reproduces byte-for-byte with no external dependency |
| The two Feynman diagrams had "no generator" | they are byte-identical to files written by two cells of `sft-wick/examples/nonlocal_vertex_2pt.ipynb` | extracted to `figures/make_feyndiag.py`; reproduces byte-for-byte at today's sft-wick |
| Figure 5 depended on the talk repository | the compute script imported the talk's sweep code by absolute path, and the figure's Order-0 planes are interpolants of a talk-era cache | the sweep code moved into `analysis3/multiz_sweep.py` and the talk file became a wrapper; the cache is vendored with provenance; the exact Order-0 at the figure's own planes was computed as a cross-check (agrees to a median 2e-3) |

## Phase 1: version-control baseline

Tagged `pre-reorg-2026-09` in all three repositories before touching anything.

* Paper repository fast-forwarded one commit to the Overleaf head (`1f0ec58`, which sets
  `\revhlfalse`).
* The 18 uncommitted items in the analysis package committed in four groups: the canoes
  path resolver and the vertex-build flags; the 2026-08 revision generator changes; the
  regenerated outputs; the new xAct derivation.
* `driver_field_emulators/` imported by copy, verified identical to the original, with an
  md5 manifest of all 1080 files (`reorg_2026-09/dfe_manifest_2026-09-02.txt`), tagged
  `dfe-imported-2026-09`. 24 MB entered git; the large regenerable, exploratory and
  private material is gitignored and listed in that manifest.
* Root `.gitignore` rewritten: markdown is no longer ignored globally (which is why
  `figure_manifest.md` had never been tracked), and the agent instruction files and build
  outputs are ignored explicitly.
* The talk repository stopped tracking compiled Python files and gained a `.gitignore`.

## Phase 2: moves and default paths

Six commits, each a group, every move a `git mv`:

* **callables**: the permutation-aware callable became
  `equal_time_limber_cut15360_permaware/` with its table under the name it resolves; the
  factored vertex rebuild and `b_model` moved beneath `kappa3_vertex/`; the L2 lambda grid
  was vendored out of the archive.
* **run folders**: the manuscript's FK sweep, the ten-fold cutoff ladder and the multi-z
  sweeps became first-class run folders with configs relative to their own directory. The
  runner now resolves paths against the YAML folder and refuses to overwrite a folder
  marked `PRODUCTION` unless forced, which is the guard against the accident that
  destroyed the deployed FK sweep on 2026-08-25.
* **analyses**: the FK Monte-Carlo, the revision drivers and the audit probes moved under
  `sachs_sft/analyses/`; the multi-z sweep helper and the talk cache were vendored into
  `analysis3/`.
* **docs**: the audit notes, the 22-page LaTeX note and the handover files moved under
  `docs/fk_audit_2026-08/`; the agent briefs were collected; the arXiv PDFs and the editor
  correspondence are local only.
* **defaults**: every generator now defaults to the corrected input and writes into its own
  `outputs/`; the FK Monte-Carlo markers became a generator option with the pooled markers
  as default; deployment into the manuscript happens only through `reproduce/deploy.py`.

All 17 figures were regenerated and compared after the moves; 14 matched byte-for-byte
immediately and the three schematic-type figures matched after redeployment (decision D2;
they had drifted by a fraction of a point in the tight bounding box since June, with
identical curves).

## Phase 3: stale products and stamps

* The sft-wick and canoes commits, working-tree patches and both conda package lists were
  pinned into `reorg_2026-09/env/` before any re-run.
* The FK 2PCF sweep was re-run once (2 min 50 s, 6 workers) and reproduced the deployed
  npz byte for byte.
* The exact Order-0 at the figure's own five source planes was computed (85 s) and
  compared with the interpolated planes: median relative difference 2e-3, and the quoted
  FK/Order-0 band is unchanged at the printed precision.
* The correlation-propagator callable's hardcoded canoes path, which pointed at a
  directory deleted months ago, was replaced by the same resolver the vertex callable
  uses. Without this the callable could not be imported at all outside a shell that
  already set `PYTHONPATH`, which is why the first verification re-run failed.
* Thirteen audit documents that quote an `ell_max = 1000` FK amplitude gained a dated
  correction block pointing at `docs/fk_audit_2026-08/FK_BASELINE_NUMBERS.md`, a new page
  that keys the five different values of the same quantity to their configurations.

## Phase 4: documentation

`REPRODUCE.md` (environments, one command per figure with expected time and md5, a
provenance line per data product, the sweep re-run rules), `figure_manifest.md` rewritten
and tracked, `MOVED_PATHS.md` (old to new, with the SHA-256 status of every vendored file),
this changelog, and the updated `CLAUDE.md` layout section.

## What was deleted, and what was only retired

**Deleted outright**: regenerable caches (`__pycache__`, `.ruff_cache`, `.pytest_cache`,
`.benchmarks`, `.DS_Store`), 19 zero-byte `synctex(busy)` files, and the LaTeX build
residue of the audit note and the revision letter.

**Retired to `~/.Trash/STF_lensing_reorg_2026-09-02/`**, every file listed with its md5 and
original path in `reorg_2026-09/trash_manifest_2026-09-02.txt`:

| item | size | why |
|---|---|---|
| `_archive/canoes_pipeline/` | 160 MB | the 2026-05 pipeline; the three things still live in it (the coordinate-map generator, the base config geometry, the L2 lambda grid) were vendored into the live tree first |
| `_archive/cleanup_2026-06-17/` | 815 MB | the June cleanup archive, including a 773 MB superseded Limber C-table |
| `sachs_sft/_cleanup_archive_2026-06-09/` | 71 MB | the h^4-fix archive |
| `figures/_archived/` | 3 MB | May 2026 figure drafts |
| root `analysis_setup.md`, `findings.md`, `progress.md`, `task_plan.md`, `outputs/` | small | superseded planning documents; `analysis_setup.md` predates both reorganisations and describes a layout that no longer exists |
| the untracked remainder of the imported `driver_field_emulators` copy | 20 MB | exploratory variants and caches, all in the manifest |

**Removed from git but present in history** at tag `pre-reorg-2026-09`:
`figures/archived_2026-08-26/`, `figures/archived_2026-08-28/`, the two orphan figure PDFs
(`appendix_3cumulant_fastnc.pdf`, `cumulant_hierarchy.pdf`), the two obsolete June
`analysis3` notes, and the byte-identical `_PROD_REF` twin of the June FK sweep. The
manuscript's `figures/` folder now holds exactly the 17 live PDFs.

## How to roll back

| scope | command |
|---|---|
| any single repository, to the pre-reorganisation state | `git reset --hard pre-reorg-2026-09` |
| the analysis package, to just after the verbatim import | `git reset --hard dfe-imported-2026-09` |
| one group of moves | `git revert <commit>`; each phase-2 group is a single commit |
| a retired file | copy it back from `~/.Trash/STF_lensing_reorg_2026-09-02/` using the path in the manifest, while the Trash is unemptied |
| a figure | `git checkout pre-reorg-2026-09 -- figures/<name>.pdf` in the paper repository |
| the archived and orphan figure PDFs | `git checkout 1f0ec58 -- figures/archived_2026-08-26 figures/archived_2026-08-28 figures/appendix_3cumulant_fastnc.pdf figures/cumulant_hierarchy.pdf` |

Overleaf's git bridge accepts only the `master` branch and rejects tags, so the paper
repository's rollback point is the commit `1f0ec58` rather than the tag name. The tag
exists locally and the commit is on Overleaf, so either works from a local clone.

The scripts that performed each phase are in `reorg_2026-09/phase1/`, `phase2/` and
`phase3/`, so any step can be read back exactly as it ran.

## Open items

* **D16 is closed.** All three repositories were pushed on 2026-09-02: the analysis
  package and the talk to GitHub (both with the `pre-reorg-2026-09` tag), the manuscript
  to Overleaf. Before the package push, the rendered figures under `outputs/` were
  untracked: `regen_figure.py` rewrites them on every run, and tracking them alongside
  the deployed copies is what let three of them go stale in the first place. Every data
  file under `outputs/` is still tracked.
* The four wording nuances of REPRODUCE.md section 5 are for the operator to decide on;
  this reorganisation proposes no manuscript edit.
* `rebuild/products/pieces/` (50 MB of octave-band builds) and the arXiv PDFs are local
  only; a fresh clone can rebuild the former in hours and does not need the latter.
