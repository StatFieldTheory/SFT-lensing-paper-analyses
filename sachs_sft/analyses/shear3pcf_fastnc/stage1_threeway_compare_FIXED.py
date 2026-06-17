#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""stage1_threeway_compare_FIXED.py
===================================
Identical to stage1_threeway_compare.py but reads the FRESH-REBUILT canoes
output stage1_ours_v2_FIXED.npz (group-zeroing fix applied) instead of
stage1_ours_v2.npz. Writes outputs/stage1_threeway_compare_FIXED.npz.

The reference (stage1_bkappa_phase_ref.npz) and fastnc (stage1_fastnc.npz) are
READ-ONLY and unchanged -- only the canoes-under-test side is rebuilt.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"

ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
fnc = np.load(OUT / "stage1_fastnc.npz", allow_pickle=True)
can = np.load(OUT / "stage1_ours_v2_FIXED.npz", allow_pickle=True)
# also load the pre-fix canoes for a BEFORE/AFTER column
can_old = np.load(OUT / "stage1_ours_v2.npz", allow_pickle=True)

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
)}

print("=" * 100)
print("THREE-WAY COMPARISON (FIXED): shear 3PCF natural components |Gamma^mu|")
print("  ref = B_kappa-phase reference | fnc = fastnc tree | can = canoes zeta_D (FIXED)")
print("  can_old = canoes zeta_D (pre-fix) for BEFORE/AFTER")
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
          f"can_FIXED={fi_can:.4e}  can_old={fi_can_old:.4e}")
    print(f"  can/ref:  BEFORE={fi_can_old/fi_ref:.3f}   AFTER(FIXED)={fi_can/fi_ref:.3f}")
    print(f"  ref/can_FIXED={fi_ref/fi_can:.3f}  ref/fnc(self-check)={fi_ref/fi_fnc:.3f}")

    store["gamma_arcmin"].append(g); store["phi_deg"].append(p)
    store["ref"].append(vref); store["fnc"].append(vfnc)
    store["can"].append(vcan); store["can_old"].append(vcan_old)
    store["fi_ref"].append(fi_ref); store["fi_fnc"].append(fi_fnc)
    store["fi_can"].append(fi_can); store["fi_can_old"].append(fi_can_old)
    store["r_ref_fnc"].append(fi_ref / fi_fnc)
    store["r_ref_can"].append(fi_ref / fi_can)
    store["r_fnc_can"].append(fi_fnc / fi_can)
    store["can_over_ref"].append(fi_can / fi_ref)
    store["can_old_over_ref"].append(fi_can_old / fi_ref)

co_before = np.array(store["can_old_over_ref"])
co_after = np.array(store["can_over_ref"])
print("\n" + "=" * 100)
print("VERDICT (frame-invariant can/ref, should -> ~1 if the shear CLOSES):")
print(f"  BEFORE fix: median={np.median(co_before):.3f}  range=[{co_before.min():.3f},{co_before.max():.3f}]")
print(f"  AFTER fix:  median={np.median(co_after):.3f}  range=[{co_after.min():.3f},{co_after.max():.3f}]")
print(f"  ref/fnc self-check: median={np.median(store['r_ref_fnc']):.3f}")

outpath = OUT / "stage1_threeway_compare_FIXED.npz"
np.savez(
    outpath,
    gamma_arcmin=np.array(store["gamma_arcmin"]), phi_deg=np.array(store["phi_deg"]),
    ref=np.array(store["ref"]), fnc=np.array(store["fnc"]),
    can=np.array(store["can"]), can_old=np.array(store["can_old"]),
    fi_ref=np.array(store["fi_ref"]), fi_fnc=np.array(store["fi_fnc"]),
    fi_can=np.array(store["fi_can"]), fi_can_old=np.array(store["fi_can_old"]),
    r_ref_fnc=np.array(store["r_ref_fnc"]), r_ref_can=np.array(store["r_ref_can"]),
    r_fnc_can=np.array(store["r_fnc_can"]),
    can_over_ref=co_after, can_old_over_ref=co_before,
    comps=np.array(COMPS),
    note=("THREE-WAY shear 3PCF compare with the FIXED canoes kernel "
          "(group-zeroing fix). can=stage1_ours_v2_FIXED.npz, "
          "can_old=stage1_ours_v2.npz (pre-fix). ref/fnc READ-ONLY. "
          "can_over_ref BEFORE vs AFTER quantifies whether the LOW-path fix "
          "closes the shear; HIGH (Limber) path dominates zeta_D so the fix "
          "is expected sub-dominant."),
)
print(f"\nsaved -> {outpath}")
