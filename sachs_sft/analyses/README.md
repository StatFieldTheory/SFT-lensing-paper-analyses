# analyses/

One folder per analysis (the role the old `analysis_1/1b/2/3/4` played). Each
analysis folder:
- references callables BY PATH from `../callables/<kind>/<impl>/` (explicit, no env).
- holds its own config(s), driver, outputs, and figures — self-contained.
- states which C_propagator impl and which kappa3_vertex impl it uses.

The old analyses live in `scripts/_archive/canoes_pipeline/analysis_*` for
reference; port forward only what is needed, re-wired to the new callable folders.
