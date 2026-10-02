# analysis1 — Order-0 lensing 2PCFs: SFT corr_op (non-Limber) vs PyCCL FKEM

> **Superseded and archived on 2026-10-02.** `plot_analysis1_O0_vs_pyccl_fkem.py`,
> `build_pyccl_fkem_reference.py` and `pyccl_xi_shear_reference_fkem.npz` were moved to
> `_archive/cleanup_2026-10-02/` (index: `reproduce/archived_inputs.json`). The manuscript's Order-0
> comparison is in `r1_aligned/`. The text below describes the earlier comparison and is kept as a record.

**Question.** Does the sft-wick L2 expansion's **Order-0** lensing 2-point
functions reproduce PyCCL when *both sides are non-Limber*? Earlier the PyCCL
reference was pure Limber, so the corr_op (full non-Limber) curve could only
match it at small gamma (high ell, where Limber is exact). Switching the
reference to the **FKEM non-Limber** integrator makes this an apples-to-apples
validation across the whole gamma range.

## Production figure — 4 panels, no ratio panels (analysis3 style)

`plot_analysis1_O0_vs_pyccl_fkem.py` → `outputs/analysis1_O0_vs_pyccl_fkem.{png,pdf}`

2x2 grid of the four Order-0 observables, styled exactly like analysis3 via the
shared `_plot_style.py` (thick gray sign-aware background line = reference;
coloured sign-aware markers = SFT; filled = +, hollow = −; shared top legend).

| Panel | Observable | SFT component combo | PyCCL ref key |
|---|---|---|---|
| κκ | `xi_kappa` | `(0,0)·+1` | `xi_kappa` (`P_L`) |
| ξ₊ | `xi_plus` | `(1,1)·+1 + (2,2)·+1` | `xi_plus` (`d^L_{2,2}`) |
| ξ₋ | `xi_minus` | `(1,1)·+1 + (2,2)·−1` | `xi_minus` (`d^L_{2,-2}`) |
| κγ_t | `xi_kappa_gamma` | `(0,1)·−1` | `xi_kappa_gamma` (`d^L_{0,2}`) |

Two curves per panel, both **non-Limber**:

| Curve | Source | C propagator |
|---|---|---|
| **PyCCL FKEM (non-Limber)** | `pyccl_xi_shear_reference_fkem.npz` | `C_L^{κκ}` from a z_s=5 `CMBLensingTracer`, FKEM non-Limber for ℓ<250 + Limber for ℓ≥250, summed full-sky with `P_L` / Wigner-d |
| **SFT O0 (corr_op, non-Limber)** | `../../sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz` | corr_op — full **non-Limber** |

Shared cosmology (z_s=5, Ω_m=0.3161, h=0.6711, n_s=0.97, σ₈=0.81) and γ grid
(`logspace(0.5, 5000, 40)` arcmin).

**Result (2026-06-04).** SFT corr_op Order-0 tracks PyCCL **FKEM** to sub-percent
where each observable is large: SFT/PyCCL @ small γ = 0.998 (κκ), 0.995 (ξ₊),
1.008 (ξ₋), 1.000 (κγ_t). Residual departures appear only at large γ where each
observable crosses zero (ratio ill-defined) and on the tiny negative tail. With
both sides non-Limber the agreement now extends across the full γ range, not
just the small-γ / high-ℓ regime the old Limber reference matched.

## Why FKEM instead of Limber?

The old reference (`pyccl_xi_shear_reference_extended.npz`, built by the
archived `build_pyccl_shear_reference.py`) called `ccl.angular_cl` with the
default `l_limber=-1` → **pure Limber for every ℓ**. The non-Limber correction
is a *low-ℓ* effect (FKEM/Limber = +11.3% at ℓ=2, +5.6% at ℓ=3, <0.1% by ℓ≈30
for this z_s=5 kernel), which in real space lives at *large* γ. So FKEM only
changes the large-γ tail of the reference; the sub-arcminute regime (high ℓ,
Limber exact) is untouched. Running FKEM up to ℓ=4500 would be pointless and
slow, so the build uses FKEM below `--l-limber` (default 250) and Limber above;
the crossover ratio (printed at build time) is smooth (0.9999 → 1.0000).

## Build the FKEM reference

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python build_pyccl_fkem_reference.py
```

Writes `pyccl_xi_shear_reference_fkem.npz` (ells, Cl_kappa, gamma_arcmin,
xi_kappa, xi_plus, xi_minus, xi_kappa_gamma, cosmology, `l_limber`,
`non_limber_method="FKEM"`).

## Plot

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python plot_analysis1_O0_vs_pyccl_fkem.py
```

Consumes the corr_op O0 sweep + the FKEM reference; writes the 2x2 figure and
prints a per-observable SFT/PyCCL summary.

## Provenance / style

- The PyCCL reference uses a full-sky Legendre/Wigner-d sum, NOT
  `ccl.correlation` (FFTlog Hankel issues at γ ≲ 1/ℓ_max); see the archived
  `build_pyccl_shear_reference.py` docstring.
- `_plot_style.py` is identical to analysis3's (shared paper-figure visual
  language; mpl + numpy only).

## Files

- `plot_analysis1_O0_vs_pyccl_fkem.py` — **production** 4-panel plotter (FKEM ref).
- `build_pyccl_fkem_reference.py` — builds the FKEM non-Limber reference npz.
- `pyccl_xi_shear_reference_fkem.npz` — PyCCL FKEM reference (current).
- `_plot_style.py` — self-contained paper-figure style (mpl + numpy only).
- `outputs/` — figures.

### Superseded (kept for provenance)

- `plot_analysis1_kk_vs_pyccl.py` — old **single-panel, κκ-only** plotter that
  compared corr_op against the **Limber** reference. Superseded by the 4-panel
  FKEM figure above.
- `pyccl_xi_shear_reference_extended.npz` — old **Limber** reference (default
  `l_limber=-1`). Kept; the FKEM npz is the production reference.
