"""Digit-by-digit comparison: clean-room vs fk_expect_exact.py."""
import sys
from pathlib import Path
_MC = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete")
sys.path.insert(0, str(_MC))
sys.path.insert(0, str(_MC / "review"))
import numpy as np

import _bootstrap
_DS, _BG, _CORE = _bootstrap.wire()
import fk_expect_exact as EX
import cleanroom_exact as CR

GAMMA, N = 1.0, 1000
cosg = float(np.cos(np.deg2rad(GAMMA / 60.0)))

# ---------------- reference side ------------------------------------------
def ref_side(apply_anchor):
    grid = EX.build_grid(n_lambda=N, apply_anchor=apply_anchor)
    Zt = np.array([_DS.zeta6(cosg, float(l)) for l in grid.lam])
    loc = EX.fk_reference(grid, Zt)
    out = {"local": float(loc[0, 3] + loc[3, 0]) / 2.0 * 2.0}   # M+M.T -> [0,3]
    out["local_03"] = float(loc[0, 3])
    for cal in ("nominal", "smeared"):
        V, A, Q, rho = EX.node_stats(grid, cosg, SIG, calibrate=cal)
        T1, T2 = EX.expectation(grid, cosg, SIG, V=V, A=A, Q=Q, rho=rho)
        FK1 = T1 + T1.T
        FKt = (T1 + T2) + (T1 + T2).T
        out[cal] = {"T1": float(FK1[0, 3]), "T1pT2": float(FKt[0, 3]),
                    "share": float(FK1[0, 3] / FKt[0, 3])}
    return out

# ---------------- clean-room side -----------------------------------------
def cr_side():
    o = CR.run(GAMMA, N, SIG, verbose=False)
    return {"local_03": o["local"]["total"],
            "nominal": {"T1": o["nominal"]["T1"], "T1pT2": o["nominal"]["total"],
                        "share": o["nominal"]["T1"] / o["nominal"]["total"]},
            "smeared": {"T1": o["smeared"]["T1"], "T1pT2": o["smeared"]["total"],
                        "share": o["smeared"]["T1"] / o["smeared"]["total"]}}

def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b))

for SIG in (8.0, 4.0):
    print("=" * 78)
    print(f"gamma = {GAMMA}'   N = {N}   sigma_lambda = {SIG}")
    R = ref_side(apply_anchor=False)
    Ra = ref_side(apply_anchor=True)
    C = cr_side()
    print(f"  ref local ref (apply_c0=False) = {R['local_03']:.10e}")
    print(f"  ref local ref (apply_c0=True ) = {Ra['local_03']:.10e}")
    print(f"  CR  local ref                  = {C['local_03']:.10e}"
          f"   rel diff vs ref(False) = {rel(C['local_03'], R['local_03']):.3e}")
    for cal in ("nominal", "smeared"):
        for key in ("T1", "T1pT2", "share"):
            a, b, c = R[cal][key], Ra[cal][key], C[cal][key]
            print(f"  {cal:8s} {key:6s} ref(c0=F)={a:+.10e}  ref(c0=T)={b:+.10e}  "
                  f"CR={c:+.10e}   rel(CR,refF)={rel(c, a):.3e}")
    print()
