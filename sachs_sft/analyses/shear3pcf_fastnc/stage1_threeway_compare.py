#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""stage1_threeway_compare.py
==============================
Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

THE VERDICT: three-way comparison of the cosmic-shear 3PCF natural components
Gamma^0..3 from THREE independent objects, all on the SAME pinned configs:

  1. ref     = stage1_bkappa_phase_ref.npz   (DECISIVE definition-level reference:
               B_kappa -> Gamma via the standard flat-sky spin-2 phase projection;
               independent of BOTH canoes zeta_D and fastnc 2DFFTLog)
  2. canoes  = stage1_ours_v2.npz            (canoes driving-field zeta_D route -- under test)
  3. fastnc  = stage1_fastnc.npz             (fastnc tree-level 2DFFTLog)

Computes ref/fastnc, ref/canoes, fastnc/canoes per-component (magnitude) AND
frame-invariant sqrt(sum|Gamma|^2).  The ref-vs-fastnc ratio is the SELF-CHECK
(must be ~O(1), ideally <~30%, to trust the reference).
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"

ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
fnc = np.load(OUT / "stage1_fastnc.npz", allow_pickle=True)
can = np.load(OUT / "stage1_ours_v2.npz", allow_pickle=True)

# reference configs (the ones we actually built)
ref_gamma = ref["gamma_arcmin"]
ref_phi = ref["phi_deg"]
CONFIGS = list(zip(ref_gamma.tolist(), ref_phi.tolist()))
COMPS = ["Gamma0", "Gamma1", "Gamma2", "Gamma3"]


def grid_index(npz, g, p):
    """Locate (g, p) in a (n_gamma, n_phi) grid npz."""
    ig = int(np.argmin(np.abs(npz["gamma_arcmin"] - g)))
    ip = int(np.argmin(np.abs(npz["phi_deg"] - p)))
    return ig, ip


def get_vec(npz, g, p, is_grid):
    if is_grid:
        ig, ip = grid_index(npz, g, p)
        return np.array([npz[c][ig, ip] for c in COMPS])
    else:
        # ref npz is a flat list over configs
        idx = CONFIGS.index((g, p))
        return np.array([npz[c][idx] for c in COMPS])


def frame_inv(v):
    return np.sqrt(np.sum(np.abs(v) ** 2))


rows = []
print("=" * 92)
print("THREE-WAY COMPARISON: shear 3PCF natural components |Gamma^mu|")
print("  ref = B_kappa-phase reference (definition) | fnc = fastnc tree | can = canoes zeta_D")
print("=" * 92)

store = {
    "gamma_arcmin": [], "phi_deg": [],
    "ref": [], "fnc": [], "can": [],
    "fi_ref": [], "fi_fnc": [], "fi_can": [],
    "r_ref_fnc": [], "r_ref_can": [], "r_fnc_can": [],
    "rc_ref_fnc": [], "rc_ref_can": [], "rc_fnc_can": [],
}

for (g, p) in CONFIGS:
    vref = get_vec(ref, g, p, is_grid=False)
    vfnc = get_vec(fnc, g, p, is_grid=True)
    vcan = get_vec(can, g, p, is_grid=True)
    fi_ref, fi_fnc, fi_can = frame_inv(vref), frame_inv(vfnc), frame_inv(vcan)

    aref, afnc, acan = np.abs(vref), np.abs(vfnc), np.abs(vcan)
    # per-component magnitude ratios
    rc_ref_fnc = aref / np.maximum(afnc, 1e-300)
    rc_ref_can = aref / np.maximum(acan, 1e-300)
    rc_fnc_can = afnc / np.maximum(acan, 1e-300)

    tag = "  <-- PRIMARY" if (g, p) == (10.0, 60.0) else ""
    print(f"\ngamma={g:.0f}'  phi={p:.0f}deg{tag}")
    print(f"  frame-inv:  ref={fi_ref:.4e}  fnc={fi_fnc:.4e}  can={fi_can:.4e}")
    print(f"     ratios:  ref/fnc={fi_ref/fi_fnc:.3f}   ref/can={fi_ref/fi_can:.3f}"
          f"   fnc/can={fi_fnc/fi_can:.3f}")
    for i, c in enumerate(COMPS):
        print(f"    {c}: |ref|={aref[i]:.3e} |fnc|={afnc[i]:.3e} |can|={acan[i]:.3e}"
              f"  | ref/fnc={rc_ref_fnc[i]:.2f} ref/can={rc_ref_can[i]:.2f} "
              f"fnc/can={rc_fnc_can[i]:.2f}")

    store["gamma_arcmin"].append(g)
    store["phi_deg"].append(p)
    store["ref"].append(vref); store["fnc"].append(vfnc); store["can"].append(vcan)
    store["fi_ref"].append(fi_ref); store["fi_fnc"].append(fi_fnc); store["fi_can"].append(fi_can)
    store["r_ref_fnc"].append(fi_ref / fi_fnc)
    store["r_ref_can"].append(fi_ref / fi_can)
    store["r_fnc_can"].append(fi_fnc / fi_can)
    store["rc_ref_fnc"].append(rc_ref_fnc)
    store["rc_ref_can"].append(rc_ref_can)
    store["rc_fnc_can"].append(rc_fnc_can)

# ---- summary / verdict ----------------------------------------------------
r_ref_fnc = np.array(store["r_ref_fnc"])
r_ref_can = np.array(store["r_ref_can"])
r_fnc_can = np.array(store["r_fnc_can"])
print("\n" + "=" * 92)
print("SELF-CHECK (ref vs fastnc tree -- must be O(1) to trust the reference):")
print(f"  frame-inv ref/fnc over configs: {np.array2string(r_ref_fnc, precision=3)}")
print(f"    -> median={np.median(r_ref_fnc):.3f}  range=[{r_ref_fnc.min():.3f},{r_ref_fnc.max():.3f}]")
print("\nVERDICT inputs (frame-invariant):")
print(f"  ref/can (definition vs canoes zeta_D): median={np.median(r_ref_can):.3f}"
      f"  range=[{r_ref_can.min():.3f},{r_ref_can.max():.3f}]")
print(f"  fnc/can (fastnc vs canoes zeta_D):     median={np.median(r_fnc_can):.3f}"
      f"  range=[{r_fnc_can.min():.3f},{r_fnc_can.max():.3f}]")
print(f"  can/ref (canoes / definition):         median={np.median(1/r_ref_can):.3f}"
      f"  range=[{(1/r_ref_can).min():.3f},{(1/r_ref_can).max():.3f}]")

outpath = OUT / "stage1_threeway_compare.npz"
np.savez(
    outpath,
    gamma_arcmin=np.array(store["gamma_arcmin"]),
    phi_deg=np.array(store["phi_deg"]),
    ref=np.array(store["ref"]), fnc=np.array(store["fnc"]), can=np.array(store["can"]),
    fi_ref=np.array(store["fi_ref"]), fi_fnc=np.array(store["fi_fnc"]),
    fi_can=np.array(store["fi_can"]),
    r_ref_fnc=r_ref_fnc, r_ref_can=r_ref_can, r_fnc_can=r_fnc_can,
    rc_ref_fnc=np.array(store["rc_ref_fnc"]),
    rc_ref_can=np.array(store["rc_ref_can"]),
    rc_fnc_can=np.array(store["rc_fnc_can"]),
    comps=np.array(COMPS),
    note=("Three-way shear 3PCF natural-component comparison. ref=B_kappa-phase "
          "definition reference (stage1_bkappa_phase_ref.npz), fnc=fastnc tree "
          "(stage1_fastnc.npz), can=canoes zeta_D (stage1_ours_v2.npz). Ratios are "
          "magnitude ratios; frame-inv = sqrt(sum|Gamma^mu|^2). ref/fnc is the "
          "self-check; ref/can and fnc/can are the verdict."),
)
print(f"\nsaved -> {outpath}")
