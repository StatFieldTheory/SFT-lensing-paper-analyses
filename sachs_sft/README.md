# sachs_sft — clean sft-wick driving pipeline (since 2026-05-30)

Fresh start replacing the tangled `scripts/canoes_pipeline` (archived at
`scripts/_archive/canoes_pipeline`). The organizing principle is **traceability by
directory structure**: each callable implementation is a single self-contained
folder, so there is no ambiguity about which products/conventions belong to which
implementation.

## What sft-wick consumes

sft-wick's config file is driven by exactly **two kinds of callable objects**:

1. **Correlation propagator** `C` — the 2-point cumulant `⟨Φ Φ⟩`.
2. **Non-local κ³ vertex** `K` — the 3-point driving-field cumulant.

Each kind has **two validated implementations** (different folders under
`callables/`). The config picks one folder per kind.

## Layout

```
sachs_sft/
├── callables/
│   ├── C_propagator/
│   │   ├── corr_op/          # impl 1: canoes corr_op C_ell(λ',λ'') table → propagator
│   │   └── limber/           # impl 2: canoes Limber-version C → propagator
│   └── kappa3_vertex/
│       ├── R_contracted/     # impl 1: R-contracted-already κ³ (windowed)
│       └── equal_time_limber/# impl 2: (non-contracted) equal-time Limber κ³
├── sftwick_outputs/          # one folder per concrete sft-wick RUN (L2 config + output)
│   ├── 2PCF/                 # observable=[phi(x),phi(y)] 2-point; C-only O0 + F/K post-Born O[0,2]
│   │   {C_corr_op_O0, C_limber_O0, C_corr_op_K_contracted_FF/FK, C_corr_op_K_limber_FK}/
│   └── 3PCF/                 # (reserved) genuine 3-point observable [phi(x),phi(y),phi(z)]
├── analyses/                 # one folder per analysis (like the old analysis_1/2/3/4)
└── scripts/                  # ONLY necessary shared/common code (D_callable, kappa2_callable,
                              # run_FF_single, gen_sftwick_configs, F_tensor + kappa2 NPZ)
```

## Rules

- **One folder = one implementation.** ALL files an implementation needs (the
  callable module, its product NPZ(s) or build recipe, its convention notes, its
  validation) live inside that folder. Nothing shared implicitly across folders.
- **No env-var product resolution.** A config points at a specific callable folder
  by path; products are explicit, not resolved from `SFT_WICK_*` env defaults.
- **Each callable folder carries its own README** stating: the sft-wick contract it
  satisfies (coordinate, measure, units, already-R-contracted?), its provenance
  (which archived build produced it), and its validation status.
- `scripts/` holds only genuinely-shared utilities — not analysis logic, not
  callable implementations.

## Status (2026-05-30)

Skeleton only. The four callable implementations are to be (re)built together,
pulling validated pieces from `scripts/_archive/canoes_pipeline` as needed:
- C_propagator: corr_op + canoes-Limber (both previously validated).
- kappa3_vertex: R-contracted (windowed) + equal-time-Limber (both previously
  provided/validated).
