# analysis3 — Order-2 NLO decomposition: O0 vs FF vs FK

**Question.** How do the two Order-2 post-Born channels (**FF** and **FK**) compare to
the Order-0 Born 2-point signal, across the lensing observables? This is the *internal*
structure of the L2 NLO expansion (no external reference — analysis1 already validated
O0 against PyCCL).

```
Full = Order-0 + Order-2(FF) + Order-2(FK)
  O0 : Born 2-point
  FF : two F-vertex insertions      -> nonlinear-propagation correction
  FK : one F + one K(kappa3) vertex -> 3-point-cumulant contribution
```
KK vanishes at Order-2 (six tilde-φ legs cannot pair in MSR), so FF + FK are the
*only* Order-2 channels.

> **Why FK appears in κκ but not in shear** (the apparent κ↔γ "same information"
> puzzle): see [`WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md`](WHY_FK_BREAKS_KAPPA_GAMMA_EQUIVALENCE.md).
> Short version: the κ=γ equivalence is a *linear-theory* statement; the κ³ FK term is
> non-Gaussian + post-Born, and the Wigner-3j selection rule confines its leakage to the
> spin-0 Ricci (convergence) channel, suppressing it in the spin-2 Weyl (shear) channels.
>
> **Why FK_κκ is flat while FK in ξ±/κγ_t rises toward large scales** (the per-panel
> γ-trends): see [`FK_GAMMA_TRENDS_ACROSS_PANELS.md`](FK_GAMMA_TRENDS_ACROSS_PANELS.md).
> Short version: κκ is a real squeezed-κ³ signal (separate-universe long-mode modulation,
> flat within the long-mode coherence angle); the spin-2 channels have no signal, only a
> low-ℓ residual that the large-γ Wigner-d kernels lift off the floor.

## Inputs (three vertex-isolated sweeps, all corr_op C, identical geometry)

| Curve | Source NPZ | order rows used |
|---|---|---|
| **O0** | `../../sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz` | order 0 |
| **FF** | `../../sftwick_outputs/2PCF/C_corr_op_K_limber_FF/xi_C_corr_op_K_limber_FF.npz` | order 2 |
| **FK** | `../../sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz` | order 2 |

Each vertex-type-isolated run carries only its own order rows (the `vertex_types`
filter drops the zero-vertex Born term, so FF/FK runs have **order 2 only**). They
share the 40-pt γ grid + t_final=2313.029 + n_gauss=24, so `Full = O0 + FF + FK`
sums coherently.

## Observables + sign convention
Four panels: κκ `[(0,0)]`, ξ_+ `[(1,1)+(2,2)]`, ξ_- `[(1,1)−(2,2)]`,
κγ_t `[(0,1)]`. **κγ_t uses weight +1 on (0,1)** (corr_op convention: the corr_op
(0,1) path already carries the flipped sign), vs −1 for the archived Born-Limber path.

## Run
```
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python plot_analysis3_nlo_decomposition.py
```
Writes `outputs/analysis3_nlo_O0_FF_FK.{png,pdf}` and prints per-observable O0/FF/FK
ranges + FF/O0, FK/O0 ratios at 0.5′.

## Style / provenance
Ported from `scripts/_archive/canoes_pipeline/analysis_3/plot_paper_analysis3.py` via
the self-contained `_plot_style.py` (Full = thick gray sign-aware background line;
O0/FF/FK = filled/hollow sign-aware markers — filled +, hollow −).

## Files
- `plot_analysis3_nlo_decomposition.py` — driver/plotter.
- `_plot_style.py` — paper-figure style (mpl + numpy only).
- `outputs/` — figures.
