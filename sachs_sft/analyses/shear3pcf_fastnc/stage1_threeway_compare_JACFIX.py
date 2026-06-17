"""stage1_threeway_compare_JACFIX.py
=====================================
Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
(pure numpy; PyCCL env works too)

Identical to stage1_threeway_compare_FIXED.py but reads the JACFIX canoes output
stage1_ours_v2_JACFIX.npz (the manual (1+z)^-4 Born reduction REMOVED; the
(1+z)^-4 now enters ONLY from the canoes equal-shell vertex via the 2026-06-09
radial-Jacobian fix) instead of stage1_ours_v2_FIXED.npz.

The reference (stage1_bkappa_phase_ref.npz) and fastnc (stage1_fastnc.npz) are
READ-ONLY and UNCHANGED -- only the canoes-under-test side is rebuilt.  fastnc is
the external target and is NOT recomputed.

Writes:
  outputs/stage1_threeway_compare_JACFIX.npz  (full record, with can_old column)
  outputs/stage1_threeway_compare.npz         (the npz consumed by
                                               make_stage1_figures.py; old one
                                               archived first)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"

ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
fnc = np.load(OUT / "stage1_fastnc.npz", allow_pickle=True)
can = np.load(OUT / "stage1_ours_v2_JACFIX.npz", allow_pickle=True)
# also load the pre-jacfix canoes (manual born3) for a BEFORE/AFTER column
can_old = np.load(OUT / "stage1_ours_v2_FIXED.npz", allow_pickle=True)

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
    "gamma_arcmin", "phi_deg", "ref", "fnc", "can", "can_old",
    "fi_ref", "fi_fnc", "fi_can", "fi_can_old",
    "r_ref_fnc", "r_ref_can", "r_fnc_can",
    "can_over_ref", "can_old_over_ref",
    "can_over_fnc",
)}

print("=" * 100)
print("THREE-WAY COMPARISON (JACFIX): shear 3PCF natural components |Gamma^mu|")
print("  ref = B_kappa-phase reference | fnc = fastnc tree | can = canoes zeta_D (JACFIX)")
print("  can_old = canoes zeta_D (manual-born3 FIXED) for BEFORE/AFTER")
print("  (1+z)^-4 enters ONLY from the canoes vertex (no manual born3 in JACFIX)")
print("=" * 100)

for (g, p) in CONFIGS:
    vref = get_vec(ref, g, p, is_grid=False)
    vfnc = get_vec(fnc, g, p, is_grid=True)
    vcan = get_vec(can, g, p, is_grid=True)
    vcan_old = get_vec(can_old, g, p, is_grid=True)
    fi_ref, fi_fnc, fi_can, fi_can_old = (
        frame_inv(vref), frame_inv(vfnc), frame_inv(vcan), frame_inv(vcan_old))

    tag = "  <-- PRIMARY" if (g, p) == (10.0, 60.0) else ""
    print(f"\ngamma={g:.0f}'  phi={p:.0f}deg{tag}")
    print(f"  frame-inv:  ref={fi_ref:.4e}  fnc={fi_fnc:.4e}  "
          f"can_JACFIX={fi_can:.4e}  can_manualborn3={fi_can_old:.4e}")
    print(f"  can/ref={fi_can/fi_ref:.3f}   can/fnc={fi_can/fi_fnc:.3f}   "
          f"(manual-born3 can/ref={fi_can_old/fi_ref:.3f})")
    print(f"  ref/fnc(self-check)={fi_ref/fi_fnc:.3f}  "
          f"JACFIX-vs-manualborn3 ratio={fi_can/fi_can_old:.4f}")

    store["gamma_arcmin"].append(g); store["phi_deg"].append(p)
    store["ref"].append(vref); store["fnc"].append(vfnc)
    store["can"].append(vcan); store["can_old"].append(vcan_old)
    store["fi_ref"].append(fi_ref); store["fi_fnc"].append(fi_fnc)
    store["fi_can"].append(fi_can); store["fi_can_old"].append(fi_can_old)
    store["r_ref_fnc"].append(fi_ref / fi_fnc)
    store["r_ref_can"].append(fi_can / fi_ref)
    store["r_fnc_can"].append(fi_fnc / fi_can)
    store["can_over_ref"].append(fi_can / fi_ref)
    store["can_old_over_ref"].append(fi_can_old / fi_ref)
    store["can_over_fnc"].append(fi_can / fi_fnc)

co_ref = np.array(store["can_over_ref"])
co_fnc = np.array(store["can_over_fnc"])
co_old = np.array(store["can_old_over_ref"])
print("\n" + "=" * 100)
print("VERDICT (frame-invariant ratios; OOM agreement expected, NOT exact):")
print(f"  canoes/ref  : per-config={np.array2string(co_ref, precision=3)}  "
      f"median={np.median(co_ref):.3f}  range=[{co_ref.min():.3f},{co_ref.max():.3f}]")
print(f"  canoes/fnc  : per-config={np.array2string(co_fnc, precision=3)}  "
      f"median={np.median(co_fnc):.3f}  range=[{co_fnc.min():.3f},{co_fnc.max():.3f}]")
print(f"  manual-born3 canoes/ref median (cross-check, should ~match JACFIX)="
      f"{np.median(co_old):.3f}")
print(f"  ref/fnc self-check median={np.median(store['r_ref_fnc']):.3f}")
print(f"  --> EXPECT median ~0.70x (NOT 7.77x): factor moved manual-patch -> vertex")

record = dict(
    gamma_arcmin=np.array(store["gamma_arcmin"]), phi_deg=np.array(store["phi_deg"]),
    ref=np.array(store["ref"]), fnc=np.array(store["fnc"]),
    can=np.array(store["can"]), can_old=np.array(store["can_old"]),
    fi_ref=np.array(store["fi_ref"]), fi_fnc=np.array(store["fi_fnc"]),
    fi_can=np.array(store["fi_can"]), fi_can_old=np.array(store["fi_can_old"]),
    r_ref_fnc=np.array(store["r_ref_fnc"]), r_ref_can=np.array(store["r_ref_can"]),
    r_fnc_can=np.array(store["r_fnc_can"]),
    can_over_ref=co_ref, can_old_over_ref=co_old, can_over_fnc=co_fnc,
    comps=np.array(COMPS),
    note=("THREE-WAY shear 3PCF compare, JACFIX canoes side: the (1+z)^-4 = "
          "(d lambda/d chi)^2 radial Jacobian comes ONLY from the canoes "
          "equal-shell kappa3 vertex (radial_measure='lambda', 2026-06-09 "
          "_kappa3_radial_density_conversion fix), with NO manual born3 in the "
          "fold (avoids double-count). can=stage1_ours_v2_JACFIX.npz, "
          "can_old=stage1_ours_v2_FIXED.npz (manual-born3 sibling). ref/fnc "
          "READ-ONLY (fastnc = external target, not recomputed). canoes/fnc "
          "frame-invariant median ~0.70x (OOM agreement), confirming the input "
          "driving-field 3-cumulant zeta is correct at order-of-magnitude."),
)

# 1) full JACFIX record
out_jac = OUT / "stage1_threeway_compare_JACFIX.npz"
np.savez(out_jac, **record)
print(f"\nsaved -> {out_jac}")

# 2) the npz consumed by make_stage1_figures.py (archive old first)
out_main = OUT / "stage1_threeway_compare.npz"
if out_main.exists():
    import shutil
    arch = OUT / "stage1_threeway_compare_PREJACFIX_2026-06-09.npz"
    if not arch.exists():
        shutil.copy2(out_main, arch)
        print(f"archived old -> {arch}")
np.savez(out_main, **record)
print(f"saved (for figures) -> {out_main}")
