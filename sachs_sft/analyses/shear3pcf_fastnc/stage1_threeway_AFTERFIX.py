#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""stage1_threeway_AFTERFIX.py
================================
Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

VERIFY THE SPIN-2 SUBSET CLOSES AFTER THE SCALAR FIX.

This is the AFTER-FIX three-way comparison.  The canoes side
(stage1_ours_v2_FIXED.npz) is a FRESH rebuild of the shear 3PCF natural
components Gamma^0..3 with the CURRENT canoes source, which carries TWO
corrections relative to the pre-fix stage1_ours_v2.npz baseline:

  (1) GROUP-ZEROING fix in _spin_aware_three_pt.py: the (l2,l3)-keyed
      group-batched m-sum no longer reads count_g / the W3j reshape from the
      group-FIRST triple (which has l1 < |s1| for the s1=+/-2 channels and so
      m_count=0 zeroed the whole group).  This restores the dropped (l1>=|s1|)
      cells of zeta_D.
  (2) UNITS fix in kappa3.py: the "physical" zeta scaling changed h^6 -> h^4
      across compute_kappa3_zeta_table / _sigma3_high / _mod_zeta_equal_time /
      _mod_sigma3_high / _bare_bl_table.  The equal-shell zeta density is
      natively (h/Mpc)^4 = A(a)^3 (P.P)/chi^4, so physical units multiply by
      h^4 (NOT h^6).  This is a CHANNEL-INDEPENDENT factor h^2 = 0.6711^2 =
      0.4504 multiplying every zeta channel identically.

The prompt's prediction: the (1+z)^4 x h^2 correction is channel-independent,
so the spin-2 frame-invariant canoes/ref MUST drop from [4.1, 6.1, 17.2, 9.5]
toward ~1.  This script reports the AFTER values and EXPLICITLY isolates any
residual phi-dependent J_2 angular factor that survives the scalar fix.

WHAT 'can/ref' MEASURES, AND WHAT THE J_2 RESIDUAL IS
-----------------------------------------------------
ref (stage1_bkappa_phase_ref.npz) builds Gamma^mu from the GENUINE 2D Fourier
integral of B_kappa with the centroid spin-2 phase exp[2i s_j (phi_{l_j} -
alpha_j)].  The phi_{l_j} integration against B_kappa produces the FULL spin-2
angular (J_2-Bessel) weighting of the triangle.

canoes (our side) builds Gamma^mu = -fold(zeta_D) * derived_phase(mu, phi),
where fold(zeta_D) is the REAL radial fold of the equal-shell driving-field
modulus cumulant (its phi-dependence enters only through the cosine-triple
t3 = 2 gamma sin(phi/2)) and derived_phase is the analytic great-circle ->
centroid spin-2 rotation.  The explicit J_2 angular weighting that the
reference carries is NOT reproduced cell-by-cell; it is approximated by the
zeta equal-shell collapse.

Therefore:
  * A CHANNEL- and PHI-INDEPENDENT mis-normalization (the h^6 vs h^4 slip)
    shows up as a CONSTANT offset in can/ref across all configs.
  * Any LEFTOVER phi-dependence in can/ref AFTER removing that constant is the
    residual spin-2 J_2 angular factor -- the genuine remaining spin-2 question
    (great-circle equal-shell collapse vs the exact flat-sky J_2 projection).

We quantify the residual as (can/ref) / median(can/ref): a flat ~1 vector means
the spin-2 subset CLOSES up to overall normalization and there is NO surviving
J_2 residual; a vector that still trends with phi/gamma is the J_2 residual,
reported explicitly.

Writes: outputs/stage1_threeway_AFTERFIX.npz
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"

ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
fnc = np.load(OUT / "stage1_fastnc.npz", allow_pickle=True)
can = np.load(OUT / "stage1_ours_v2_FIXED.npz", allow_pickle=True)        # AFTER (both fixes)
can_old = np.load(OUT / "stage1_ours_v2.npz", allow_pickle=True)          # BEFORE (pre-fix)

ref_gamma = ref["gamma_arcmin"]
ref_phi = ref["phi_deg"]
CONFIGS = list(zip(ref_gamma.tolist(), ref_phi.tolist()))
COMPS = ["Gamma0", "Gamma1", "Gamma2", "Gamma3"]


def grid_index(npz, g, p):
    ig = int(np.argmin(np.abs(npz["gamma_arcmin"] - g)))
    ip = int(np.argmin(np.abs(npz["phi_deg"] - p)))
    return ig, ip


def get_vec(npz, g, p, is_grid):
    if is_grid:
        ig, ip = grid_index(npz, g, p)
        return np.array([npz[c][ig, ip] for c in COMPS])
    idx = CONFIGS.index((g, p))
    return np.array([npz[c][idx] for c in COMPS])


def frame_inv(v):
    return np.sqrt(np.sum(np.abs(v) ** 2))


store = {k: [] for k in (
    "gamma_arcmin", "phi_deg",
    "ref", "fnc", "can", "can_old",
    "fi_ref", "fi_fnc", "fi_can", "fi_can_old",
    "r_ref_fnc",                 # self-check (ref vs fastnc)
    "can_over_ref", "can_old_over_ref",
    "rc_can_over_ref",           # per-component can/ref (AFTER)
)}

H = float(ref["h"])
print("=" * 102)
print("THREE-WAY AFTER-FIX: shear 3PCF natural components |Gamma^mu|  (spin-2 subset closure)")
print("  ref = B_kappa-phase reference | fnc = fastnc tree | can = canoes zeta_D (FIXED: group-zero + h^4)")
print(f"  channel-independent units factor applied by the fix: h^2 = {H**2:.4f}")
print("=" * 102)

for (g, p) in CONFIGS:
    vref = get_vec(ref, g, p, is_grid=False)
    vfnc = get_vec(fnc, g, p, is_grid=True)
    vcan = get_vec(can, g, p, is_grid=True)
    vcan_old = get_vec(can_old, g, p, is_grid=True)
    fi_ref, fi_fnc, fi_can, fi_can_old = (
        frame_inv(vref), frame_inv(vfnc), frame_inv(vcan), frame_inv(vcan_old))
    aref, acan = np.abs(vref), np.abs(vcan)
    rc = acan / np.maximum(aref, 1e-300)        # per-component can/ref (AFTER)

    tag = "  <-- PRIMARY" if (g, p) == (10.0, 60.0) else ""
    print(f"\ngamma={g:.0f}'  phi={p:.0f}deg{tag}")
    print(f"  frame-inv:  ref={fi_ref:.4e}  fnc={fi_fnc:.4e}  "
          f"can_FIXED={fi_can:.4e}  can_old={fi_can_old:.4e}")
    print(f"  can/ref  BEFORE={fi_can_old/fi_ref:7.3f}   AFTER(FIXED)={fi_can/fi_ref:7.3f}"
          f"   (drop x{fi_can_old/fi_can:.3f})")
    print(f"  self-check ref/fnc={fi_ref/fi_fnc:.3f}")
    for i, c in enumerate(COMPS):
        print(f"    {c}: |ref|={aref[i]:.3e} |can|={acan[i]:.3e}  can/ref={rc[i]:.3f}")

    store["gamma_arcmin"].append(g); store["phi_deg"].append(p)
    store["ref"].append(vref); store["fnc"].append(vfnc)
    store["can"].append(vcan); store["can_old"].append(vcan_old)
    store["fi_ref"].append(fi_ref); store["fi_fnc"].append(fi_fnc)
    store["fi_can"].append(fi_can); store["fi_can_old"].append(fi_can_old)
    store["r_ref_fnc"].append(fi_ref / fi_fnc)
    store["can_over_ref"].append(fi_can / fi_ref)
    store["can_old_over_ref"].append(fi_can_old / fi_ref)
    store["rc_can_over_ref"].append(rc)

co_before = np.array(store["can_old_over_ref"])
co_after = np.array(store["can_over_ref"])
phi_arr = np.array(store["phi_deg"])
gam_arr = np.array(store["gamma_arcmin"])

# ---- residual phi-dependent J_2 isolation --------------------------------
# Remove the channel-/phi-independent overall normalization (median can/ref);
# what remains is the residual angular (J_2) factor that survives the fix.
med_after = float(np.median(co_after))
resid = co_after / med_after          # flat ~1 => closes; trend with phi => J_2 residual

print("\n" + "=" * 102)
print("VERDICT (frame-invariant can/ref -> ~1 if the spin-2 subset CLOSES):")
print(f"  BEFORE fix:        {np.array2string(co_before, precision=3)}"
      f"  median={np.median(co_before):.3f}")
print(f"  AFTER  fix:        {np.array2string(co_after, precision=3)}"
      f"  median={med_after:.3f}  range=[{co_after.min():.3f},{co_after.max():.3f}]")
print(f"  self-check ref/fnc:{np.array2string(np.array(store['r_ref_fnc']), precision=3)}"
      f"  median={np.median(store['r_ref_fnc']):.3f}")
print("\nRESIDUAL phi-dependent J_2 angular factor (can/ref / median, AFTER fix):")
for (g, p, r) in zip(gam_arr, phi_arr, resid):
    print(f"    gamma={g:.0f}' phi={p:.0f}deg : residual = {r:.3f}")
print(f"  residual spread: min={resid.min():.3f} max={resid.max():.3f}  "
      f"(max/min = {resid.max()/resid.min():.2f}x)")
print("  -> a residual that trends with phi at FIXED gamma (the (10',30),(10',60),(10',120)")
print("     triple) is the surviving spin-2 J_2 angular factor; a flat residual means the")
print("     equal-shell great-circle collapse reproduces the exact flat-sky J_2 projection.")

# phi-trend at fixed gamma=10' (the three configs sharing gamma)
mask10 = np.isclose(gam_arr, 10.0)
if mask10.sum() >= 2:
    phis10 = phi_arr[mask10]
    res10 = resid[mask10]
    order = np.argsort(phis10)
    print(f"\n  phi-trend at gamma=10' (residual vs phi): "
          f"phi={np.array2string(phis10[order], precision=0)} -> "
          f"resid={np.array2string(res10[order], precision=3)}")

outpath = OUT / "stage1_threeway_AFTERFIX.npz"
np.savez(
    outpath,
    gamma_arcmin=gam_arr, phi_deg=phi_arr,
    ref=np.array(store["ref"]), fnc=np.array(store["fnc"]),
    can=np.array(store["can"]), can_old=np.array(store["can_old"]),
    fi_ref=np.array(store["fi_ref"]), fi_fnc=np.array(store["fi_fnc"]),
    fi_can=np.array(store["fi_can"]), fi_can_old=np.array(store["fi_can_old"]),
    r_ref_fnc=np.array(store["r_ref_fnc"]),
    can_over_ref=co_after, can_old_over_ref=co_before,
    rc_can_over_ref=np.array(store["rc_can_over_ref"]),
    residual_j2=resid, residual_j2_median_norm=med_after,
    comps=np.array(COMPS), h=H,
    note=("AFTER-FIX three-way shear 3PCF closure. canoes side rebuilt with "
          "BOTH the group-zeroing fix (_spin_aware_three_pt.py) and the "
          "h^6->h^4 units fix (kappa3.py, channel-independent h^2=0.4504). "
          "can_over_ref BEFORE=[4.1,6.1,17.2,9.5] must drop toward ~1. "
          "residual_j2 = (can/ref)/median(can/ref) isolates the surviving "
          "phi-dependent spin-2 J_2 angular factor (great-circle equal-shell "
          "collapse vs exact flat-sky J_2 projection)."),
)
print(f"\nsaved -> {outpath}")
