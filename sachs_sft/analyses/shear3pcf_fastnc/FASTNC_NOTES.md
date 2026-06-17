# FASTNC_NOTES.md

API surface, install log, and z_s=5 handling for the fastnc v1.2.3 installation.
Audience: the later agent writing the tree-level SPT bispectrum subclass (STAGE 1, STEP 2).

---

## 1. Install Steps

```bash
# 1. Create fresh venv (Python 3.9 from macOS default python3)
python3 -m venv /Users/zzhang/projects/fastnc_venv

# 2. Upgrade pip
/Users/zzhang/projects/fastnc_venv/bin/pip install --upgrade pip

# 3. Core deps (numpy, scipy, astropy, mpi4py, pandas)
#    mpi4py builds against /opt/homebrew/bin/mpicc (Apple open-mpi from brew; clang-1700)
#    Downloaded a pre-built wheel (macosx_11_0_arm64) -- no compilation needed.
MPICC=/opt/homebrew/bin/mpicc \
  /Users/zzhang/projects/fastnc_venv/bin/pip install numpy scipy astropy mpi4py pandas

# 4. Clone and checkout exact commit
git clone https://github.com/git-sunao/fastnc /Users/zzhang/projects/fastnc
cd /Users/zzhang/projects/fastnc
git checkout d94fd2a

# 5. Install fastnc from local clone
/Users/zzhang/projects/fastnc_venv/bin/pip install /Users/zzhang/projects/fastnc
```

### Resolved commit hash

```
d94fd2a  added baryon at init
```

This is the HEAD of the main branch as of the checkout date. There is no v1.2.3 tag in the
repo; the version string `1.2.3` is set in `fastnc/__init__.py`.

### Full pip freeze

```
astropy==6.0.1
astropy-iers-data==0.2026.3.16.0.53.33
fastnc @ file:///Users/zzhang/projects/fastnc
mpi4py==4.1.2
numpy==1.26.4
packaging==26.2
pandas==2.3.3
pyerfa==2.0.1.5
python-dateutil==2.9.0.post0
pytz==2026.2
PyYAML==6.0.3
scipy==1.13.1
six==1.17.0
tzdata==2026.2
```

---

## 2. z_s = 5 Spline-Edge Handling

`BispectrumBase.set_cosmology` (bispectrum.py line 292) builds the chi<->z spline on:

```python
z   = np.linspace(0, 5, 100)
chi = self.cosmo.comoving_distance(z).value * self.cosmo.h  # Mpc/h
self.z2chi = ius(z, chi)       # scipy IUS, ext default (extrapolate)
self.chi2z = ius(chi, z)
```

With `z_s = 5` the source is at the last spline knot.

CHOICE: NO PATCH. We run with the unmodified clone and `z_s = 5.0`.

REASON: `z2chi` queries are only made at `zl in [zmin_losint, z_s]` (the lens redshift array)
which is within the grid. The lensing kernel `g -> 0` smoothly as `zl -> z_s` (confirmed in the
smoke test: g(chi(4.95)) = 1.79e-10, g(chi(5.0)) = 0.0 exactly, consistent with
`g = prefactor*(1 - chil/chis)` vanishing when `chil = chis`). No extrapolation artifacts.

The `chi2z` inverse spline is NOT used in the standard `compute_kernel` / `kappa_bispectrum_direct`
path (only `z2chi` is used there). `chi2z` is unused in the pipeline we invoke.

The `chi2g` spline is built with `ext=1` (returns 0 outside bounds), so the kernel correctly
drops to zero at and beyond z_s.

ALTERNATIVE (in case issues arise): patch line 292 to `linspace(0, 5.5, 110)` in the local clone
and rebuild. This would give 5 extra knots above z_s=5, making evaluation at z_s an interior
point. Not needed based on current tests.

---

## 3. API Surface for the SPT Bispectrum Subclass

### 3.1  BispectrumBase (fastnc/bispectrum.py)

The method to implement:

```python
def matter_bispectrum_no_baryon(self, k1, k2, k3, z):
    """
    k1, k2, k3 : np.ndarray  [h/Mpc]   <- UNITS ARE h/Mpc
    z          : np.ndarray  (redshift)
    Returns    : np.ndarray, same shape as k1/k2/k3/z
    """
    raise NotImplementedError
```

The k UNITS are confirmed **h/Mpc** (not 1/Mpc). Evidence:
- `kappa_bispectrum_direct` (bispectrum.py lines 611-612) computes
  `K1, K2, K3 = ELL1/CHI, ELL2/CHI, ELL3/CHI`
  where `CHI = self.z2chi(Z)` is in `Mpc/h` (line 293: `chi = ... * self.cosmo.h`).
  So `K = ell / (Mpc/h) = ell * h/Mpc`, confirming k is in h/Mpc.

The call signature inside `kappa_bispectrum_direct` is:

```python
bm = self.matter_bispectrum(K1, K2, K3, Z, **args)
```

where `K1, K2, K3` are 2D arrays of shape `(n_ell, n_z)` and `Z` has the same shape.
`matter_bispectrum` is the public wrapper that calls `matter_bispectrum_no_baryon` then
multiplies by the baryon suppression factor (default = 1).

### 3.2  set_pklin / set_lgr (on BispectrumHalofit; same interface on the SPT subclass)

```python
bs.set_pklin(k, pklin)
# k     : np.ndarray [h/Mpc]  -- matches PCAMBz0.txt column 0
# pklin : np.ndarray [(Mpc/h)^3] at z=0 -- matches PCAMBz0.txt column 1

bs.set_lgr(z, lgr)
# z   : np.ndarray of redshifts
# lgr : np.ndarray linear growth rate D(z)/D(0) (normalized to 1 at z=0)
# Internally stored as scipy IUS self.z2lgr
```

For the SPT tree-level subclass, you do NOT need to delegate to `self.halofit` (that is only
used by BiHalofit/GilMarin). Instead, store the splines directly:

```python
from scipy.interpolate import InterpolatedUnivariateSpline as ius

class BispectrumSPT(BispectrumBase):
    config_scale = dict(ell1min=1e-1, ell1max=1e5, epmu=1e-7)  # copy from BiHalofit

    def set_pklin(self, k, pklin):
        self._pk_spline = ius(np.log(k), np.log(pklin), ext=0)
        self.has_changed = True

    def set_lgr(self, z, lgr):
        self.z2lgr = ius(z, lgr, ext=1)
        self.has_changed = True

    def _pklin_at(self, k):
        return np.exp(self._pk_spline(np.log(k)))

    def matter_bispectrum_no_baryon(self, k1, k2, k3, z):
        # See section 3.5 for the F2 formula
        ...
```

### 3.3  Full pipeline call sequence

```python
from fastnc.bispectrum import BispectrumHalofit   # or your SPT subclass
from fastnc.fastnc import FastNaturalComponents
from astropy.cosmology import wCDM
import numpy as np

# 1. Instantiate with Lmax and Lmax_diag set at construction time
#    (required: without Lmax, set_multipole_grid does nothing and decompose() fails)
bs = BispectrumHalofit(Lmax=20, Lmax_diag=20)

# 2. Feed cosmology, P(k), growth rate
bs.set_cosmology(cosmo, ns=n_s, sigma8=sigma8)
bs.set_pklin(k_h, pk_h)       # h-units
bs.set_lgr(z_arr, D_arr)

# 3. Feed source distribution (single delta plane)
bs.set_source_distribution([np.array([z_s])], [np.array([1.0])])

# 4. Compute lensing kernel
bs.compute_kernel()

# 5. Build interpolation grid over (r,u,v) and decompose into Legendre multipoles
bs.interpolate()
bs.decompose()                 # stores bL_multipole[scomb][L, i_ell, i_psi]

# 6. Set up FastNaturalComponents
t1_rad = t1_arcmin * np.pi / 10800.0
phi_arr = np.array([phi1, phi2, ...])  # opening angles in radians

fnc = FastNaturalComponents(
    Lmax=20,
    Mmax=20,
    projection='cent',   # centroid projection (recommended; see section 3.6)
    t1=t1_rad,
    phi=phi_arr,
    verbose=True,
)
fnc.set_bispectrum(bs)   # also calls set_fftgrid(), sets t2 = t1

# 7. Compute
fnc.compute()

# 8. Access results
# Gamma0, Gamma1, Gamma2, Gamma3  shape: (n_t1, n_t1, n_phi)
#   because t2 = t1 is forced (see section 3.4)
```

### 3.4  SAS geometry: t2 = t1 is FORCED

In `set_fftgrid` (fastnc.py lines 219-222):

```python
# Allocate the same bin to second side
self.t2 = self.t1
self.t2_fft = self.t1_fft
self.ell2_fft = self.ell1_fft
```

So `FastNaturalComponents` always evaluates on isoceles triangles with `t2 = t1`.
The third side is `t3 = sqrt(t1^2 + t2^2 - 2*t1*t2*cos(phi)) = t1 * sqrt(2*(1-cos(phi)))`.
For `t1 = t2 = gamma` this gives:

```
t3 = 2 * gamma * sin(phi/2)       (law of cosines, t1=t2=gamma)
```

Source: `trigutils.x1x2phi_to_x1x2x3` (trigutils.py line 82):
```python
def x1x2phi_to_x1x2x3(x1, x2, phi):
    x3 = np.sqrt(x1**2 + x2**2 - 2*x1*x2*np.cos(phi))
    return x1, x2, x3
```

For the collapsed config (phi -> 0): t3 -> 0, which is the degenerate limit.
fastnc should NOT be evaluated exactly at phi=0 (t3=0 is a degenerate triangle);
instead, extrapolate the finite-phi curve to phi -> 0.

### 3.5  F2_eff location and why you must NOT reuse it

`BispectrumGilMarin.F2_eff` (bispectrum.py lines 1038-1043):

```python
def F2_eff(self, z, k1, k2, k3, knl):
    dot = (-k3 ** 2 + k1 ** 2 + k2 ** 2) / 2
    # f2 = 5 / 7 * self.agm(k1, z, knl) * self.agm(k2, z, knl) + ...
    f2 = 5 / 7 + 2 * dot ** 2 / (7 * k1 ** 2 * k2 ** 2) - dot * (1 / k1 ** 2 + 1 / k2 ** 2) / 2
    return f2
```

This is the Gil-Marin EFFECTIVE F2 which includes nonlinear correction factors `agm`, `bgm`,
`cgm` (commented out above in the production code, simplified back to SPT form). Importantly,
it is a method of `BispectrumGilMarin` which uses `halofit.get_pkhalofit` for the NONLINEAR
power spectrum. Do NOT inherit from or call `F2_eff`; implement the TRUE SPT F2 directly.

The TRUE SPT F2 kernel to implement in the subclass:

```python
def F2(self, k1, k2, k3):
    """
    Standard SPT second-order kernel.
    k3^2 = k1^2 + k2^2 + 2*(k1.k2) so (k1.k2) = (k3^2-k1^2-k2^2)/2
    cos_theta = (k1.k2)/(k1*k2) = (k3^2-k1^2-k2^2)/(2*k1*k2)
    """
    cos_theta = (k3**2 - k1**2 - k2**2) / (2*k1*k2)
    return (5/7
            + 0.5 * cos_theta * (k1/k2 + k2/k1)
            + (2/7) * cos_theta**2)

def matter_bispectrum_no_baryon(self, k1, k2, k3, z):
    """All k in h/Mpc; z is redshift array (same shape as k1/k2/k3)."""
    D = self.z2lgr(z)   # linear growth D(z)/D(0)
    pk1 = self._pklin_at(k1) * D**2
    pk2 = self._pklin_at(k2) * D**2
    pk3 = self._pklin_at(k3) * D**2
    b  =  2 * self.F2(k1, k2, k3) * pk1 * pk2
    b += 2 * self.F2(k2, k3, k1) * pk2 * pk3
    b += 2 * self.F2(k3, k1, k2) * pk3 * pk1
    return b
```

fastnc ships NO tree-level bispectrum class (only BiHalofit, BispectrumGilMarin, BispectrumNFW1Halo).

### 3.6  Projection options

Configured at `FastNaturalComponents` construction via `projection=`:

| Key     | Description                                | Notes |
|---------|--------------------------------------------|-------|
| `'x'`   | FFT-native x-projection (default)          | Computed first; all Gamma^i are derived in this basis |
| `'cent'`| Centroid projection (arXiv:2309.08601)     | Applied via `x2cent(mu, T1, T2, PHI)` phase factor |
| `'ortho'`| Orthocenter projection                    | Applied via `x2ortho(mu, T1, T2, PHI)` phase factor |

Conversion source (fastnc.py line 728):

```python
def ortho2cent(mu, t1, t2, phi):
    NotImplementedError('This function is not validated yet')
    ...
```

WARNING: `ortho2cent` is NOT validated (the `NotImplementedError` is raised as a bare expression,
not as a `raise` statement, so it does not crash — but the function body may produce incorrect
results). Avoid the `ortho -> cent` round-trip. Use `projection='cent'` directly (applied via
`x2cent` from x-native, which IS validated by arXiv:2309.08601 equations).

The `_change_shear_projection` path we use: `x -> cent` (fastnc.py lines 429-432):
```python
elif dept == 'x' and dest == 'cent':
    self.Gamma0 *= x2cent(0, T1, T2, PHI)
    self.Gamma1 *= x2cent(1, T1, T2, PHI)
    self.Gamma2 *= x2cent(2, T1, T2, PHI)
    self.Gamma3 *= x2cent(3, T1, T2, PHI)
```

where `x2cent` (fastnc.py lines 752-783) implements equations between Eq. (15) and (16) of
arXiv:2309.08601.

### 3.7  Mmax, Lmax_diag, epmu

These are the key resolution parameters:

| Parameter    | Default | Set where                                  | Meaning |
|--------------|---------|--------------------------------------------|---------|
| `Mmax`       | None    | `FastNaturalComponents(Mmax=...)`          | Max angular Fourier mode M for phi resummation; phi resolution = Mmax harmonics |
| `Lmax`       | None    | `BispectrumBase(Lmax=...)` AND `FastNaturalComponents(Lmax=...)` | Max Legendre multipole L of kappa bispectrum; must match between bs and fnc |
| `Lmax_diag`  | =Lmax   | `BispectrumBase(Lmax_diag=...)` AND `FastNaturalComponents(Lmax_diag=...)` | Higher Lmax for diagonal (t1=t2) elements where bispectrum is squeezed |
| `epmu`       | 1e-7    | `BispectrumBase(epmu=...)` (config_scale)  | Regularization: mumax = 1 - epmu avoids the squeezed (collinear) triangle limit |

**CRITICAL**: `Lmax` and `Lmax_diag` MUST be set at `BispectrumBase.__init__` time (passed as
kwargs or config dict). If omitted, `set_multipole_grid` returns early (line 228:
`if self.config_multipole['Lmax'] is None: return 0`) and `decompose()` will fail with
`AttributeError: 'BispectrumHalofit' object has no attribute 'ELL1_multipole'`.

---

## 4. Smoke Test Results (BiHalofit, STEP 1)

Cosmology: Omega_m=0.3160919980475834, h=0.6711, w0=-1.0, sigma8=0.809, n_s=0.97
Source: single delta plane z_s=5, `set_source_distribution([array([5.])], [array([1.])])`
Pipeline: `interpolate()` -> `decompose(Lmax=20)` -> `FastNaturalComponents(Lmax=20, Mmax=20, projection='cent')`

Test points: t1=t2 (isoceles, forced), phi in {60 deg, 90 deg}

```
t1=t2=10', phi=60deg:
  Gamma0 = 7.557e-07  (real; Im ~ 5e-23, numerical noise)
  Gamma1 = 5.018e-07  (real; Im ~ -3e-23, numerical noise)
  Gamma2 = (4.484e-07) + (-9.716e-09)j
  Gamma3 = (4.484e-07) + (+9.716e-09)j   [Gamma3 = Gamma2* for isoceles]

t1=t2=10', phi=90deg:
  Gamma0 = 5.102e-07
  Gamma1 = 3.705e-07
  Gamma2 = (2.526e-07) + (-4.445e-08)j
  Gamma3 = (2.526e-07) + (+4.445e-08)j

t1=t2=100', phi=60deg:
  Gamma0 = 1.441e-08
  Gamma1 = 7.636e-09
  Gamma2 = (5.283e-09) + (-1.221e-09)j
  Gamma3 = (5.283e-09) + (+1.221e-09)j

t1=t2=100', phi=90deg:
  Gamma0 = 8.839e-09
  Gamma1 = 5.082e-09
  Gamma2 = (3.490e-09) + (-3.658e-10)j
  Gamma3 = (3.490e-09) + (+3.658e-10)j
```

All 16 Gamma^i values are finite. `SMOKE TEST: PASS`.

The imaginary parts of Gamma0, Gamma1 are O(1e-23) -- pure numerical noise from the 2DFFTLog
(expected: Gamma0, Gamma1 are real for the isoceles t1=t2 configuration due to symmetry).
Gamma2, Gamma3 are complex conjugates of each other (also expected by symmetry).

The mode-coupling cache was built to `~/.fastnc/` on the first run (441 (L,M) pairs for Lmax=Mmax=20).

---

## 5. Install Gotchas

1. `Lmax` (and `Lmax_diag`) must be passed at `BispectrumBase.__init__` time. If you pass
   them only to `FastNaturalComponents` but not to the bispectrum object, `decompose()` fails.

2. pandas is required by fastnc (used in `coupling.py` for the mode-coupling cache dataframe).
   Not listed in requirements.txt; added manually.

3. mpi4py is imported unconditionally at the top of `coupling.py`. The brew open-mpi binary is
   at `/opt/homebrew/bin/mpicc`. A pre-built wheel was available for macosx_11_0_arm64, so
   no compilation was needed. If compilation is needed: set `MPICC=/opt/homebrew/bin/mpicc`.

4. There is a harmless `FutureWarning` from pandas about DataFrame concatenation in coupling.py
   (triggered when the cache is populated). Does not affect results.

5. fastnc has no v1.2.3 git tag; the version string lives only in `__init__.py`. The HEAD
   commit at install time is `d94fd2a`.
