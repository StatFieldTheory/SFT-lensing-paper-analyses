# SFT-lensing-paper-analyses

Code and saved numerical products for *Statistical Field Theory for Weak
Gravitational Lensing*. The package maps Ricci/Weyl driving-field statistics
through the Sachs dynamics and the SFT-Wick expansion to weak-lensing
correlation functions. It includes figure generators, numerical checks,
callables, run configurations and symbolic-verification scripts.

**AI-generated documentation:** this README and the reproduction guide were
written with AI assistance and checked against the repository. Consult the
source, recorded inputs and validation results when interpreting a calculation.

## Reproduce the paper's figures and numerical summaries

Start with [REPRODUCE.md](REPRODUCE.md) for environment setup and path requirements.
Obtain the package with:

```bash
git clone https://github.com/StatFieldTheory/SFT-lensing-paper-analyses.git
cd SFT-lensing-paper-analyses
```

Run the following from this repository's root in the configured Python environment:

```bash
python reproduce/regen_figure.py --list     # inspect current generators and interpreters
python reproduce/regen_figure.py --all      # regenerate 17 figure PDFs from saved products
python reproduce/check_figures.py          # compare with the manuscript's ../figures/
python reproduce/check_numbers.py          # print 33 numerical summaries and update numbers.json
```

These commands reuse the saved scientific products. They do not rerun the full
Sachs/SFT calculations, vertex builds or Monte Carlo experiments. Runtime depends
on the machine and environment; the operator figures also evaluate their stored
callables. The comparison step requires the separate manuscript's deployed PDFs.

The current R1 results are selected by
[reproduce/active_products.json](reproduce/active_products.json). It pins 12
product selections, including their SHA256 hashes. Keep this manifest in place:
removing it selects the historical inputs. Missing or altered selected products
cause an error rather than a silent fallback.

The [cleanup validation](reproduce/provenance/cleanup_20260930_validation.json)
on 2026-09-30 regenerated all 17 figures in an isolated copy without local notes
or archived analyses. All figure contents matched, with zero raster differences;
all 33 numerical summaries were unchanged. This test reused the workstation's
dependency environments, so it is not a clean-machine installation test or a
new numerical convergence proof.

## Where to look

| Path | Contents |
|---|---|
| `reproduce/` | Figure/number entry points, active manifest and deployment checksums |
| `sachs_sft/analyses/r1_sft061/` | Current revision products, input provenance, fold and Monte Carlo drivers |
| `sachs_sft/analyses/analysis1/r1_aligned/` | Source-folded Order-0 and aligned PyCCL comparison |
| `sachs_sft/analyses/analysis3/` | Correlation and angular-spectrum figure generators |
| `sachs_sft/callables/` | Correlation-propagator and three-point-vertex implementations and tables |
| `sachs_sft/sftwick_outputs/2PCF/` | Historical calculations retained for comparison |
| `figures/` | Schematic and operator-figure generators |
| `mathematica/` | Symbolic derivations and checks |
| `reorg_2026-09/env/` | Historical environment inventories; current runtime pins are in the revision provenance |

The authoritative manuscript source is maintained separately. In the author's
workspace it is one directory above this repository. Generators write local
outputs; only `reproduce/deploy.py` copies figures into the manuscript.

The current figure-to-input mapping is listed in [figure_manifest.md](figure_manifest.md).
For full numerical recomputation, follow the workflow and limitations in
[REPRODUCE.md](REPRODUCE.md). Use fresh output directories, retain provenance,
and run folds and vertex builds sequentially to control memory use.

Repository cleanup and numerical revision history are recorded in
[CHANGELOG.md](CHANGELOG.md). Local research notes and agent working files are
excluded from the published repository.
