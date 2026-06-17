# Prompt — run 3 sachs_sft 2PCF sft-wick jobs and report results

> Hand this whole file to a FRESH session (clean context). Self-contained.
> Conversation in 中文; code/comments English. End each turn with a 中文小结.

## Goal

Run THREE already-authored sft-wick L2 configs (the wiring is done + reviewed;
you only EXECUTE them and report the resulting 2-point correlation `ξ(γ)`):

1. `2PCF/C_corr_op_O0`      — corr_op C propagator, Order-0 (Born). No κ³.
2. `2PCF/C_limber_O0`       — Born-Limber C propagator, Order-0. No κ³.
3. `2PCF/C_corr_op_K_contracted_FF` — corr_op C + R-contracted κ³ vertex, Order-[0,2],
   FF sweep (only the F local-cubic vertex active in the O2 diagrams; the K vertex is
   declared but excluded by `vertex_types: ["F"]`).

All three compute the **lensing 2-point function** `ξ(γ)` (observable
`[phi_a(x), phi_b(y)]`); "FF" is a sft-wick vertex type, not the observable point count.

## Locations (all absolute)

- Runs dir: `/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/`
  Each run folder has a ready `config_L2.yaml` (do NOT edit it) + a `README.md`.
- Runner: `/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/scripts/run_FF_single.py`
- Shared kappa2 NPZ (required arg, even for the closed-form C path):
  `/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/scripts/inputs/kappa2_grid_lcut0_low_empty.npz`
- Python: `/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python`
  (prefix `CANOES_SUPPRESS_METAL_WARNING=1`).

## How the runner consumes a config (verified)

`run_FF_single.py <kappa2_npz> <output_npz> [sigma2_npz]`:
- reads the config from env `SFT_WICK_CONFIG_YAML` (accepts an ABSOLUTE path);
- `os.chdir`s to the `scripts/` dir, but our configs use ABSOLUTE module/output/cache
  paths, so chdir does not misresolve anything;
- clears the config's expand+propagator caches at start (clean re-run);
- sets `SFT_WICK_KAPPA2_NPZ` from argv[1]; the configs' κ³/C callables use their OWN
  bundled tables (no `SFT_WICK_KAPPA3_*` env needed);
- writes to `cfg.output[0].path` (already absolute, inside the run folder) then
  `shutil.move`s to argv[2]. Pass the SAME absolute output path as argv[2] (move = no-op).
- `SFT_WICK_EXPAND_ORDERS` / `SFT_WICK_SWEEP_N_GAUSS` env can override the YAML; do NOT
  set them — let each config's own `orders` and `n_gauss=24` stand.

## Run recipe (per run)

For run 1 (`C_corr_op_O0`), from any CWD:

```bash
cd /Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft
R=sftwick_outputs/2PCF/C_corr_op_O0
CANOES_SUPPRESS_METAL_WARNING=1 \
SFT_WICK_CONFIG_YAML="$(pwd)/$R/config_L2.yaml" \
  /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python scripts/run_FF_single.py \
    "$(pwd)/scripts/inputs/kappa2_grid_lcut0_low_empty.npz" \
    "$(pwd)/$R/xi_C_corr_op_O0.npz"
```

Repeat for run 2 (`R=sftwick_outputs/2PCF/C_limber_O0`, output `xi_C_limber_O0.npz`)
and run 3 (`R=sftwick_outputs/2PCF/C_corr_op_K_contracted_FF`, output
`xi_C_corr_op_K_contracted_FF.npz`). Each config already encodes its own `orders`
and `output[0].path`; do not override.

- These can be slow (n_gauss=24, 40 γ-points, 6 component pairs; run 3 is Order-[0,2]
  so heavier). Run them in the BACKGROUND and poll; do not block.
- If a run errors, capture the FULL traceback and diagnose (do not silently retry).
  Likely first-run issues: a callable import error (check the callable's own
  `python <callable>.py` self-check runs), a table-not-found (the bundled NPZ must
  exist in the callable folder), or a sft-wick contract mismatch (the callable's
  fail-loud guard — report its message verbatim).

## What to report (per run + a comparison)

Read each output NPZ. The schema follows the archived sft-wick FF output
(`x`, `y`, `t_final`, `a`, `b`, `order`, `value` arrays — inspect keys). Extract the
**κκ component** (`a==0 & b==0`, highest `order`) as a function of γ (compute
`γ = arccos(n̂_x·n̂_y)` in arcmin from the stored `x`,`y` unit vectors), mirroring the
archived `analysis_3/plot_fk_windowed_vs_limber.py::load_kk` logic.

Report:
1. Per run: wall time; κκ(γ) at a few γ (e.g. 0.5′, 10′, 100′, 1000′); finite + sign.
2. **corr_op vs limber (runs 1 vs 2)**: the 2-point κκ A/B — same geometry, only the
   C propagator differs (full non-Limber vs Born-Limber). Quantify the difference.
3. Run 3 (O[0,2] FF): the Order-0 part should match run 1's κκ (same corr_op C, same
   Born observable); the Order-2 (F-vertex post-Born) is the new piece — report its
   size relative to Order-0.
4. A figure: κκ(γ) for all three on one log-x plot (View the PNG and critique it; do
   not stop at "the script ran"). Save under the respective run folders or a shared
   `2PCF/` figure.
5. Any anomaly (non-finite, wrong sign, a run that is suspiciously identical to
   another) — flag it; do NOT fit-to-target.

Do NOT edit the configs or callables. If a config looks wrong, REPORT it rather than
patch it (the wiring was reviewed; a real issue is worth surfacing, not silently fixing).
Commit nothing. Output artifacts (NPZ, PNG) into the run folders.

中文小结要点：三个 run 的 κκ(γ) 数值、corr_op-vs-limber 的 2 点 A/B、run3 的 Order-2 相对 Order-0 大小、以及图的批判性检查结论。
