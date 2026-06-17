# equal_time_limber kappa3 vertex -- diagnostic slices

Multi-faceted slices of the `equal_time_limber` non-local kappa3 vertex `K`,
with the **angular separation gamma** (`cos gamma_ij = n_i . n_j`) as the primary
axis. Reproduce with:

```
/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python slice_equal_time_limber.py
```

The script imports the sibling callable directly (so it slices exactly what
sft-wick consumes: the k=4 inverse-distance kNN surface in cosine-space + linear
interp in lambda), and also reads the raw NPZ grid nodes for ground truth.

## What the table is

15-node cosine grid (`cos = k/7`, k=-7..7, i.e. gamma 0..180 deg) x 16 lambda
shells (397..2327 Mpc <-> z 0.1..5.7). 1671 cosine-triples stored; canonicalized
to 351 sorted-descending unique triples. Four driving-field channels:
`zeta_TTT` (phi^3), `zeta_TTP`, `zeta_TPP`, `zeta_PPP` (psi^3). All signed.

## gamma families swept (both on the native grid)

- **Squeezed / isoceles** `(1, cos g, cos g)`: apex pair coincident, third at
  separation g. 15 nodes, g in [0, 180] deg. This is the family adjacent to the
  FK 2-point apex `cos=(1,1,1)`.
- **Equilateral** `(cos g, cos g, cos g)`: all three pairwise separated by g.
  11 nodes, g in [0, 120] deg (equilateral spherical triangles need cos g >= -1/2).

## Figures

| file | slice |
|------|-------|
| `slice_gamma_TTT_contact.png` | **headline**: log\|zeta_TTT\| vs gamma, squeezed, all 16 z-shells (viridis). |
| `slice_gamma_channels_zhi.png` | all 4 channels vs gamma at z=5.7; squeezed (solid) vs equilateral (dashed). |
| `slice_redshift_squeezed.png` | zeta vs z at fixed gamma (20/60/90/135 deg), all 4 channels. |
| `heatmap_gamma_z_squeezed.png` | 2D zeta(gamma, z) per channel (RdBu symlog). |
| `interp_fidelity_TTT.png` | kNN interp vs exact grid nodes near the apex, z=5.7. |
| `slice_figdata.npz` | raw arrays behind the curves. |

## Key findings (gamma dependence)

1. **The vertex is a near-contact object.** `zeta_TTT` at the fully-collapsed apex
   `cos=(1,1,1)` (g=0) is `2.45e-11` -- **112x larger** than the next-largest
   triple `(1, 0.714, 0.714)` at `2.18e-13`. Only **1 of 351 triples** exceeds 1%
   of the global max. Off-apex, all configs sit on a flat `~1e-13` plateau (median
   peak `7.5e-14`) with sign-oscillating, near-noise structure. So "vs gamma" =
   a sharp spike at g=0 plus a flat low floor; not a smooth fall-off.

2. **Only TTT spikes at the apex.** The cross/psi channels (TTP, TPP, PPP) have no
   contact enhancement -- they stay at the `~1e-14` plateau even at g=0. The strong
   contact term is specific to the phi^3 auto-cumulant.

3. **High-redshift object.** Every channel is `~0` for z<2 and grows steeply past
   z~3, peaking at the farthest shell z=5.7. The equal-time vertex accumulates with
   source distance.

4. **This explains the FK kappa-kappa baseline.** The FK 2-point topology colocates
   all three vertex legs at `cos=(1,1,1)` (memory: kappa3 colocation finding), i.e.
   it samples exactly the contact spike -- which is why `+3.086e-5` is a clean,
   well-defined baseline. The 3-point fan-out to wide triangles samples the ~100x
   weaker plateau.

## Interpolation caveat (boundary-validation flag)

The cosine grid is **coarse near the apex**: the first off-apex node is at g=31 deg.
Between g=0 (`-2.45e-11`) and g=31 deg (`~+2e-13`) the production kNN returns a
roughly linear ramp of a 100x drop (see `interp_fidelity_TTT.png`). Production is
safe *because* the FK 2-point query lands exactly on the g=0 node (`dist~0` ->
exact path, no interpolation). But any future consumer that queries `0 < g < 31 deg`
would ride an under-resolved interpolant of the contact spike. If 3-point work ever
probes small-but-nonzero gamma, densify the cosine grid near cos=1 first.
