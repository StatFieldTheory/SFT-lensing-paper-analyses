# SFT-lensing-paper-analyses

The reproduction package of *Statistical Field Theory for Weak Gravitational Lensing*:
every generator, callable and data product the manuscript's figures and numbers come from.
The manuscript itself lives one level up (`main.tex`, `sections/`, `figures/`, `biblio.bib`)
and is read-only from here.

**Start with [`REPRODUCE.md`](REPRODUCE.md).** Three commands rebuild all 17 figures,
check them against the deployed PDFs, and recompute every quantitative claim:

```bash
python reproduce/regen_figure.py --all      # about 90 seconds
python reproduce/check_figures.py           # IDENTICAL(raw) for all 17
python reproduce/check_numbers.py           # 35 claims, each next to its products
```

Then:

* [`figure_manifest.md`](figure_manifest.md): figure by figure, generator, inputs,
  environment, runtime.
* [`MOVED_PATHS.md`](MOVED_PATHS.md): old path to new path for the 2026-09 reorganisation,
  for anything holding a reference into this tree.
* [`CHANGELOG_REORG_2026-09.md`](CHANGELOG_REORG_2026-09.md): what changed, why, and how to
  roll it back.
* [`docs/`](docs/README.md): the 2026-08 FK audit that drove the revision, and the design
  records.

## Map

```
reproduce/            regen_figure.py, check_figures.py, check_numbers.py, deploy.py,
                      pdfcmp.py, deployed_md5.txt, numbers.json, logs/
figures/              generators for the schematic and operator figures (+ the coordinate
                      maps and the Feynman diagrams); outputs/ holds their renders
sachs_sft/
  callables/          what sft-wick is driven by: two kinds of callable object
    C_propagator/corr_op/                    correlation propagator C = <Phi Phi> + table
    kappa3_vertex/
      equal_time_limber_cut15360_permaware/  PRODUCTION vertex (permutation-aware, ell_max 15360)
      equal_time_limber/                     the superseded cut1000 vertex (DEPRECATED.md),
                                             the build script, the ell-band decomposition,
                                             the figure-7 generator and its outputs
      rebuild/                               factored vertex rebuild: LOW branch, octave bands,
                                             assembly, the fold runner, probes, products/, logs/
      b_model/                               matter-bispectrum models (tree, BiHalofit)
  sftwick_outputs/2PCF/                      one folder per sft-wick run, each with its config
    C_corr_op_O0/                              Order-0                       (PRODUCTION)
    C_corr_op_K_limber_FF/                     nonlinear propagation         (PRODUCTION)
    C_corr_op_K_limber_FK_cut15360_permfix/    the paper's FK                (PRODUCTION)
    C_corr_op_K_limber_FK/                     the superseded June FK        (SUPERSEDED)
    cutoff_ladder/                             the ten folds behind figure 17
    multiz/                                    the source-redshift sweeps behind figure 5
  analyses/           analysis1, analysis3, mc_sachs_2pt, mc_fk_complete,
                      revision_2026-08 (the drivers that built the 2026-08 figures),
                      fk_audit_2026-08 (the audit probes), shear3pcf_fastnc (backup)
  scripts/            shared code and the sweep runner
docs/                 fk_audit_2026-08 (notes + the 22-page audit note), agent_briefs,
                      the (1+z)^4 design records; papers/ and private/ are local only
mathematica/          symbolic proofs backing the manuscript's equations
reorg_2026-09/        the 2026-09 reorganisation: plan, evidence, scripts, environment pins
```

## Rules that keep this reproducible

* **One folder per callable implementation**, and configs name a callable folder by path.
  This is what makes a result traceable to one file; it replaced an environment-variable
  scheme that lost provenance twice.
* **A run folder marked `PRODUCTION` is protected.** `scripts/run_FF_single.py` refuses to
  overwrite its output unless `SFT_WICK_FORCE_OVERWRITE=1` is set. A variant config that
  inherited a production output path destroyed the deployed FK sweep once, on 2026-08-25.
  Use `rebuild/run_fk_variant.py` to fold anything experimental.
* **Generators never write into the manuscript.** Each writes into its own `outputs/`;
  `reproduce/deploy.py` is the single path into `../figures/`.
* **Never quote an FK amplitude without its multipole cutoff.** The same quantity appears
  in this repository as five different numbers; `docs/fk_audit_2026-08/FK_BASELINE_NUMBERS.md`
  is the key.
* **`--n-jobs` is a memory control**, not a speed knob: sft-wick forks with loky and every
  worker reloads the inputs. Never run two folds at once, and never fan out the canoes
  vertex builds; that cost this machine two out-of-memory reboots in August.
* English in code, comments and commits. Archive rather than delete. Do not edit
  `sections/*.tex`, `main.tex` or `biblio.bib` from here.
