# STAGE 1 fastnc-side notes (reference half)

These notes document the fastnc REFERENCE half of STAGE 1 (cosmic-shear 3PCF
natural components Gamma^0..3 on the pinned isoceles family at z_s=5). The
our-side Gamma^i is produced independently by a parallel agent; this side did
NOT touch the our-side pipeline. README.md is owned by the orchestrator and was
NOT edited (per task constraints); all fastnc-side notes live here.

## Files produced

| file | what |
|------|------|
| `stage1_fastnc_tree.py` | STEP 2 driver: tree-level SPT subclass + full pinned grid (the VALIDATION run). Interpreter `/Users/zzhang/projects/fastnc_venv/bin/python`. |
| `stage1_fastnc_bihalofit_plumb.py` | STEP 1 driver: BiHalofit large-gamma plumbing. |
| `outputs/stage1_fastnc.npz` | STEP 2 result, npz_schema-compliant. Gamma0..3 complex128 (3 gamma x 7 phi), projection='cent'. |
| `outputs/stage1_fastnc_bihalofit_plumb.npz` | STEP 1 result, schema-compliant. Gamma0..3 complex128 (4 gamma x 2 phi). |
| `outputs/stage1_fastnc_selfcheck.png` | self-diagnostic (Gamma^0 vs phi, |Gamma^0| log, Gamma^2 Re/Im). |
| `outputs/stage1_fastnc_run.log`, `outputs/stage1_bihalofit_plumb.log` | full run logs. |

## Resolved fastnc commit + deps

- fastnc commit `d94fd2a` (`d94fd2a027ceecb4f483ff1ba39b62b77d9bc651`, "added baryon at init"), version string `1.2.3` (no git tag).
- venv `/Users/zzhang/projects/fastnc_venv`; key deps: numpy 1.26.4, scipy 1.13.1,
  astropy 6.0.1, mpi4py 4.1.2, pandas 2.3.3, matplotlib 3.9.4 (installed by this
  agent for the self-check PNG only; not needed for the physics).

## PRIMARY POINT (gamma=10', phi=60deg) tree-level Gamma^0..3 (cent)

Converged setting (Lmax=40, Lmax_diag=80, Mmax=40, epmu=1e-7):

```
Gamma0 = (+6.536657e-08) + (+6.86e-24)j     # real; Im is FFT noise
Gamma1 = (+1.341805e-07) + (-1.04e-23)j     # real
Gamma2 = (+1.261109e-07) + (+4.967576e-09)j
Gamma3 = (+1.261109e-07) + (-4.967576e-09)j # = conj(Gamma2) for isoceles t1=t2
```

Production-grid value (Lmax=24, Lmax_diag=48, Mmax=24; what is stored at
`stage1_fastnc.npz[gamma=10', phi=60deg]`):

```
Gamma0 = +6.81455e-08
Gamma1 = +1.37522e-07
Gamma2 = (+1.242222e-07) + (+6.227785e-09)j
Gamma3 = (+1.242222e-07) + (-6.227785e-09)j
```

The production setting sits ~4% above the highest-resolution run for Gamma^0
(see convergence below); the converged numbers above are the reference for the
single-point gate. `npz['primary_point_gamma0..3']` stores the CONVERGED values.

## Convergence at the primary point

Reference = (Lmax=40, Lmax_diag=80, Mmax=40, epmu=1e-7). reldiff vs reference:

| Lmax | Lmax_diag | Mmax | epmu | Re Gamma0 | reldiff(G0) | reldiff(G2) |
|------|-----------|------|------|-----------|-------------|-------------|
| 20 | 20 | 20 | 1e-7 | 6.9497e-08 | 6.3e-2 | 3.8e-2 |
| 20 | 40 | 20 | 1e-7 | 6.9497e-08 | 6.3e-2 | 3.8e-2 |
| 30 | 60 | 30 | 1e-7 | 6.6775e-08 | 2.2e-2 | 8.7e-3 |
| 40 | 80 | 40 | 1e-7 | 6.5367e-08 | 0 (ref)| 0 (ref) |
| 30 | 60 | 30 | 1e-6 | 6.6691e-08 | 2.0e-2 | 8.7e-3 |
| 30 | 60 | 30 | 1e-8 | 6.6870e-08 | 2.3e-2 | 8.7e-3 |

Reading:
- **Lmax / Mmax** (raised together) is the dominant knob: Gamma^0 converges
  monotonically from below, 6.3% (Lmax=20) -> 2.2% (Lmax=30) -> 0 (Lmax=40).
  The trend is smooth and the increments shrink, so the Lmax=40 value is good to
  ~1-2%.
- **Lmax_diag** at fixed Lmax=20 made NO difference (20 vs 40 identical to all
  printed digits): the squeezed-diagonal multipole tail is not activated at this
  near-equilateral config (t1=t2=t3=10' at phi=60deg). It only matters for very
  squeezed phi (small phi / large phi).
- **epmu** (1e-6 / 1e-7 / 1e-8) moves Gamma^0 by <0.3% and Gamma^2 by <0.1%: the
  collinear regularization is not limiting at the primary point.

Production grid uses (24, 48, 24, 1e-7) as a cost/accuracy compromise (~3-4% of
the Lmax=40 reference). The orchestrator may want to bump to Lmax>=30 if a
sub-percent reference is needed; the script's `conv_settings` makes that a
one-line change.

## STEP 1 BiHalofit plumbing numbers (large gamma)

Pinned plumbing point gamma=150', phi=60deg (BiHalofit, NONLINEAR, cent):

```
Gamma0 = (+6.036363e-09)
Gamma1 = (+2.891700e-09)
Gamma2 = (+1.668943e-09) + (-5.759318e-10)j
Gamma3 = (+1.668943e-09) + (+5.759318e-10)j
```

Nearby large-gamma points at phi=60deg: G0(100')=1.441e-8, G0(200')=3.048e-9,
G0(300')=1.049e-9. The G0(100',60deg)=1.4414e-8 reproduces the prior smoke-test
value 1.441e-8 EXACTLY (geometry/units/projection plumbing confirmed stable).

**Plumbing cross-check (tree vs nonlinear at large gamma):** at gamma=150',
phi=60deg the tree-level Gamma^0 = 6.079e-9 vs BiHalofit Gamma^0 = 6.036e-9, i.e.
ratio tree/BiHalofit = 1.007 (<1% apart). This is the expected behavior: at large
gamma tree-level ~ nonlinear. This confirms the units/geometry are right end to
end. (At arcmin gamma the two diverge by ~10x -- the smoke value at 10',60deg was
7.557e-7 for BiHalofit vs 6.5e-8 tree-level, an ~11.6x nonlinear boost -- so do
NOT compare BiHalofit at arcmin.)

## h/Mpc units assertion (CONFIRMED)

- PCAMBz0.txt is fed DIRECTLY to `set_pklin` with NO h-conversion: k in column 0
  is [h/Mpc], P in column 1 is [(Mpc/h)^3] at z=0.
- Inside `BispectrumSPT.matter_bispectrum_no_baryon`, k1,k2,k3 arrive in h/Mpc
  (verified from source: `bispectrum.py:293` `chi = comoving_distance(z)*cosmo.h`
  so chi is Mpc/h; `bispectrum.py:611-612` `K = ELL/CHI` so k = ell/(Mpc/h) =
  h/Mpc; the `matter_bispectrum_no_baryon` docstring at `bispectrum.py:488-491`
  states "k in h/Mpc unit"). The h-unit P(k) spline is evaluated at the h-unit k
  with no rescaling. The subclass asserts finiteness/positivity of the handed k
  and documents the convention in its docstring.
- sigma8 derived from the h-unit P(k) (top-hat R=8 Mpc/h) = **0.808988**, matching
  the smoke-test value 0.809 -- an independent confirmation the P(k) is in the
  expected h-units (a stray factor of h would have shifted sigma8 by ~0.67^1.5).
- Equilateral bispectrum sanity: B_delta(k=0.1 h/Mpc, z=0) = 4.982e7 (Mpc/h)^6
  matches the analytic `6 * F2_equilat(=0.2857) * P(0.1)^2` = 4.982e7 to 0.005%.

## Sanity: tree-level Gamma^0 BELOW BiHalofit smoke at arcmin

Confirmed. At gamma=10', phi=60deg: tree-level Gamma^0 = 6.5e-8 (converged) sits
BELOW the BiHalofit smoke value 7.557e-7 (a factor ~11.6). Nonlinear power boosts
the arcmin 3PCF well above tree-level, as expected. Tree-level is the correct,
clean comparison object for the our-side zeta (tree_phi).

## The exact fastnc Gamma^i definition + centroid projection (source refs)

Method: Sugiyama+2024 (arXiv:2407.01798), Schneider-Lombardi natural components
Gamma^0..3 (4 complex, flat-sky), SAS isoceles t1=t2=gamma FORCED, opening phi,
t3 = 2 gamma sin(phi/2).

- mu -> Gamma map and multipole indices (`fastnc/fastnc.py` `compute()`,
  line ~363): for each angular mode M and mu in {0,1,2,3},
  `(m, n) = [(M-3,-M-3), (-M-1,M-1), (M+1,-M-3), (M-3,-M+1)][mu]`, with
  mu=0->Gamma0, 1->Gamma1, 2->Gamma2, 3->Gamma3. The (-1)^m,(-1)^n signs come
  from `J_m = (-1)^m J_{-m}`; normalization `/(2*pi)^3`. The Gamma are resummed
  over phi as `Gamma = sum_M GammaM exp(i M phi) / (2 pi)`
  (`fastnc.py:401-405`).
- Projection (`fastnc.py` `_change_shear_projection`, lines 411-432, and `x2cent`,
  lines 752-789): we use `projection='cent'` (centroid), applied from the
  FFT-native x-projection by multiplying each Gamma^mu by the phase factor
  `x2cent(mu, T1, T2, PHI)`. x2cent implements the equations between Eq. (15) and
  (16) of arXiv:2309.08601:

  ```
  v = t1 + t2 exp(-i phi);     q1 = v / conj(v)
  v = -2 t1 + t2 exp(-i phi);  q2 = v / conj(v)
  v = t1 - 2 t2 exp(-i phi);   q3 = v / conj(v)
  Gamma^0 *= q1 q2 q3 exp(+3i phi)
  Gamma^1 *= conj(q1) q2 q3 exp(+1i phi)
  Gamma^2 *= q1 conj(q2) q3 exp(+3i phi)
  Gamma^3 *= q1 q2 conj(q3) exp(-1i phi)
  ```

  We took `cent` DIRECTLY from x-native (the validated x2cent path); we did NOT
  use the `ortho2cent` round-trip, which is flagged unvalidated in-source
  (`fastnc.py:728`, bare `NotImplementedError(...)` expression).
- Matter bispectrum: STANDARD SPT TREE-LEVEL (NOT Gil-Marin F2_eff,
  `bispectrum.py:1038-1051`, NOT BiHalofit). Implemented in a ~40-line subclass
  `BispectrumSPT(BispectrumBase)` in `stage1_fastnc_tree.py`:
  `B_delta = 2 F2(k1,k2) P(k1) P(k2) + 2 cyclic` with the TRUE SPT kernel
  `F2 = 5/7 + (1/2) cos_theta (k1/k2 + k2/k1) + (2/7) cos_theta^2`,
  `cos_theta = (k3^2 - k1^2 - k2^2)/(2 k1 k2)`. z-dependence via `D(z)^4` (each
  linear P carries D(z)^2; growth spline fed by `set_lgr`). fastnc applies its own
  built-in single-plane convergence Limber kernel and 2DFFTLog.

## Geometry / FFT-grid gotcha encountered + fix

fastnc's `get_tuned_fftgrid` (`fastnc.py:885-898`) requires the requested t1 array
to be (a) >= 2 points (it uses `np.diff(t)[0]`) and (b) UNIFORM-IN-LOG (it places
the real-space grid at `start + nskip*arange(t.size)` from a log-uniform stencil).
The pinned gamma grid [10,50,150]' is NOT log-uniform (ratios 5, 3), so passing
all three at once trips an internal AssertionError (`t1_fft` lands on [10,50,250]).
Fix: evaluate EACH gamma on its own 2-point LOG-UNIFORM stencil [gamma, gamma*3]
and take the diagonal isoceles element. Verified (BiHalofit, gamma=10', phi=60deg)
that the diagonal Gamma is companion-stencil-independent to ~0.15% (10' paired with
100' vs 150' gave 7.558e-7 vs 7.547e-7), and that it reproduces the smoke-test
value 7.557e-7. This is the exact, robust way to evaluate the non-log-uniform
pinned grid.

## z_s=5 finiteness

Confirmed finite near the top of the los grid (no patch needed, matching the prior
agent's finding): kernel g(z=4.99)=3.56e-11, g(z=5.0)=0.0 exactly (g ~ (1-chi_l/chi_s)
vanishes at the source plane). All Gamma^i are finite across the full grid.

## phi -> 0 (collapsed slice) handling

fastnc was NOT evaluated at phi=0 (degenerate t3=0 triangle, delicate spin phases).
The finite-phi curve is provided down to phi=10deg (the smallest pinned phi) for the
orchestrator to extrapolate to phi->0. The Gamma^0 curve is smooth and monotone
toward small phi (gamma=10': 6.81e-8 at 60deg rising to 1.22e-7 at 10deg) and
finite, so the phi->0 extrapolation is well-posed. Convergence in Mmax / Lmax_diag /
epmu was checked at the primary point; the orchestrator should re-check at the
smallest phi if it pushes the extrapolation aggressively (small phi activates the
Lmax_diag squeezed tail).
```
