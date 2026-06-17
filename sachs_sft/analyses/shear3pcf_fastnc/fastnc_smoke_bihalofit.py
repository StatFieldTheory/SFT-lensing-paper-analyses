#!/Users/zzhang/projects/fastnc_venv/bin/python
"""
fastnc_smoke_bihalofit.py
=========================
Interpreter: /Users/zzhang/projects/fastnc_venv/bin/python

BiHalofit smoke test for the fastnc installation.
PURPOSE: confirm that mpi4py loads, the ~/.fastnc mode-coupling cache builds,
the pipeline runs end-to-end, and Gamma^0..Gamma^3 are finite.
This is STEP 1 (plumbing check only); it uses the built-in nonlinear
BiHalofit, not the SPT tree-level subclass.

Cosmology: astropy FlatwCDM / wCDM with
    Omega_m = 0.3160919980475834, h = 0.6711, n_s = 0.97
Single source plane: z_s = 5

z_s=5 spline-edge handling: BispectrumBase.set_cosmology (bispectrum.py:292)
builds the chi<->z spline on linspace(0, 5, 100).  With z_s=5 the source is
at the grid endpoint; scipy IUS evaluates correctly at the boundary (no
extrapolation needed).  The lensing kernel integral runs from zmin=1e-4 to
zs=5 so z2chi is queried only at in-range points.  We verified no
extrapolation artifacts by inspecting the kernel near z=5 (g -> 0 as
expected).  NO patch needed; we run unmodified with z_s=5 and document.

Outputs:
    outputs/fastnc_bihalofit_smoke.npz
"""

import os
import sys
import numpy as np
from astropy.cosmology import wCDM

# --- Cosmology ---------------------------------------------------------------
Omega_m = 0.3160919980475834
h       = 0.6711
n_s     = 0.97
sigma8  = 0.809  # derived from PCAMBz0.txt (k in h/Mpc, P in (Mpc/h)^3)
w0      = -1.0
# Ode0 for flat cosmology: Ode0 = 1 - Omega_m (Omega_k=0)
Ode0    = 1.0 - Omega_m

cosmo = wCDM(H0=100*h, Om0=Omega_m, Ode0=Ode0, w0=w0,
             meta={'n': n_s, 'sigma8': sigma8})
print(f"Cosmology: Om0={cosmo.Om0}, h={cosmo.h}, w0={cosmo.w0}")
print(f"  sigma8={sigma8}, n_s={n_s}")

# --- Source distribution -----------------------------------------------------
z_s = 5.0   # single delta plane
print(f"Source plane: z_s = {z_s}")

# --- Load linear P(k) from PCAMBz0.txt (h-units) ----------------------------
pkfile = '/Users/zzhang/projects/canoes/examples/data/PCAMBz0.txt'
data = np.loadtxt(pkfile)
k_h  = data[:, 0]   # k [h/Mpc]
pk_h = data[:, 1]   # P [(Mpc/h)^3]
print(f"Loaded PCAMBz0.txt: {k_h.size} points, k in [{k_h.min():.3e}, {k_h.max():.3e}] h/Mpc")

# Linear growth rate: D(z)/D(0) for flat wCDM
# Simple numerical integration; fastnc only needs it via set_lgr
from scipy.integrate import odeint
def E(z):
    a = 1.0 / (1.0 + z)
    return np.sqrt(Omega_m / a**3 + Ode0 * a**(-3*(1+w0)))

def growth_integrand(y, z):
    a = 1.0 / (1.0 + z)
    return [y[1], -(2.5 - 1.5*Omega_m / (E(z)**2 * a**3) * 0.5) * y[1] / (1+z)
            + 1.5 * Omega_m / (E(z)**2 * a**3) / (1+z)**2 * y[0]]

# Actually use a simpler ODE form in a variable
from scipy.integrate import quad
def growth_factor_unnorm(z):
    """D(z) proportional to integral (unnormalized)"""
    def integrand(zp):
        return (1.0 + zp) / E(zp)**3
    val, _ = quad(integrand, z, 100.0)
    return E(z) * val

z_arr = np.linspace(0, 5, 100)
D_arr = np.array([growth_factor_unnorm(z) for z in z_arr])
D_arr /= D_arr[0]  # normalize D(z=0) = 1
print(f"Linear growth rate D(z=0)={D_arr[0]:.4f}, D(z=5)={D_arr[-1]:.6f}")

# --- Import fastnc -----------------------------------------------------------
print("\nImporting fastnc...")
import fastnc
from fastnc.bispectrum import BispectrumHalofit
from fastnc.fastnc import FastNaturalComponents
print(f"fastnc version: {fastnc.__version__}")

# --- Instantiate BiHalofit ---------------------------------------------------
# Lmax and Lmax_diag must be passed at construction time so that
# set_multipole_grid creates the ELL1_multipole grid (bispectrum.py:228).
# Without Lmax, decompose() fails with AttributeError.
print("\nInstantiating BispectrumHalofit(Lmax=20, Lmax_diag=20)...")
bs = BispectrumHalofit(Lmax=20, Lmax_diag=20)

# Set cosmology (also initializes halofit internals)
bs.set_cosmology(cosmo, ns=n_s, sigma8=sigma8)

# Feed the linear P(k) (h-units, as-is)
bs.set_pklin(k_h, pk_h.copy())

# Feed the linear growth rate
bs.set_lgr(z_arr, D_arr)

# Set single source plane at z_s = 5
bs.set_source_distribution([np.array([z_s])], [np.array([1.0])])

# Compute lensing kernel
print("Computing lensing kernel...")
bs.compute_kernel()

# Check kernel near z=5 (spline-edge check)
import warnings
zl_check = np.array([4.5, 4.8, 4.9, 4.95, 5.0])
g_check  = bs.chi2g_dict['0'](bs.z2chi(zl_check))
print("  Kernel g(chi(z)) near z=5 [spline-edge check]:")
for z_c, g_c in zip(zl_check, g_check):
    print(f"    z={z_c:.2f}  g={g_c:.6e}")

# --- Interpolate and decompose -----------------------------------------------
print("\nInterpolating bispectrum (builds on interpolation grid)...")
bs.interpolate()

print("Decomposing into multipoles (Lmax=20)...")
bs.decompose()

# --- FastNaturalComponents ---------------------------------------------------
print("\nSetting up FastNaturalComponents (Lmax=20, Mmax=20, projection='cent')...")
# Use two test angular separations in radians
t1_arcmin = np.array([10.0, 100.0])   # test points: 10 arcmin, 100 arcmin
t1_rad    = t1_arcmin * np.pi / 10800.0
phi_arr   = np.array([np.pi/3, np.pi/2])  # 60 deg, 90 deg opening angles

fnc = FastNaturalComponents(
    Lmax=20,
    Mmax=20,
    projection='cent',
    t1=t1_rad,
    phi=phi_arr,
    verbose=True,
)
fnc.set_bispectrum(bs)

print("\nRunning FastNaturalComponents.compute() [first run builds ~/.fastnc cache; may take minutes]...")
fnc.compute()

# --- Extract Gamma^i ---------------------------------------------------------
print("\nSmoke test results:")
print("  Gamma^0..3 shape:", fnc.Gamma0.shape)  # (nt1, nt1, nphi)
print("  t1 = t2 (arcmin):", t1_arcmin)
print("  phi (rad):", phi_arr, " = (", np.degrees(phi_arr), "deg)")

# Print all values (diagonal t1=t2, opening angles phi)
print("\n  --- Gamma values at (t1=t2, phi) ---")
labels = ['Gamma0', 'Gamma1', 'Gamma2', 'Gamma3']
gammas = [fnc.Gamma0, fnc.Gamma1, fnc.Gamma2, fnc.Gamma3]
results = {}

for i_t, t_am in enumerate(t1_arcmin):
    for i_phi, phi_d in enumerate(np.degrees(phi_arr)):
        print(f"\n  t1=t2={t_am:.0f}', phi={phi_d:.0f}deg:")
        for lbl, G in zip(labels, gammas):
            val = G[i_t, i_t, i_phi]
            results[f'{lbl}_t{t_am:.0f}_phi{phi_d:.0f}'] = val
            print(f"    {lbl} = {val:.6e}  (Re={val.real:.6e}, Im={val.imag:.6e})")

# --- Save outputs ------------------------------------------------------------
outdir = os.path.join(os.path.dirname(__file__), 'outputs')
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, 'fastnc_bihalofit_smoke.npz')

np.savez(outfile,
    t1_arcmin=t1_arcmin,
    t2_arcmin=t1_arcmin,   # t2 = t1 forced
    phi_rad=phi_arr,
    phi_deg=np.degrees(phi_arr),
    Gamma0=fnc.Gamma0,
    Gamma1=fnc.Gamma1,
    Gamma2=fnc.Gamma2,
    Gamma3=fnc.Gamma3,
    **{k: np.array([v.real, v.imag]) for k, v in results.items()},
)
print(f"\nSaved smoke results to: {outfile}")

# --- Final status -----------------------------------------------------------
all_finite = all(
    np.all(np.isfinite(G)) for G in gammas
)
print(f"\nAll Gamma^i finite: {all_finite}")
print("SMOKE TEST: PASS" if all_finite else "SMOKE TEST: FAIL (non-finite values)")
