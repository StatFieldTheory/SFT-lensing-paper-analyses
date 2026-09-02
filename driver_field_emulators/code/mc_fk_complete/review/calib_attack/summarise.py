import json, glob, os
import numpy as np
D = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete/review/calib_attack"
FFREE = {   # <kappa1^3>/local, channel (0,0,0), gamma=1.0152', sigma=4, N=1999
 "smeared": 0.9995, "kernel": 1.0000, "kernel_fk": 1.0910, "smeared_V": 1.0084,
 "causal": 3.2643, "causal2": 0.8161, "anticausal2": 1.1927, "trunc3.0": 1.0905,
 "rtrunc2.0": 1.0192, "rtrunc3.0": 1.0058, "geo": 1.0894, "nominal": 2.0487,
 "smeared_Vk": None}
for tag in ("g1", "g173", "g53"):
    p = os.path.join(D, f"_{tag}.json")
    if not os.path.exists(p):
        continue
    d = json.load(open(p)); fv = d["fold"]
    print(f"\n=== gamma = {d['gamma']}'   fold = {fv:.6e} ===")
    print(f"{'calib':>12} {'Ffree':>7} | {'s=16':>8} {'s=8':>8} {'s=4':>8} {'s=2':>8}"
          f" | {'[16,8]':>8} {'[8,4]':>8} {'[4,2]':>8}")
    for nm, sd in d["res"].items():
        sig = sorted((float(k) for k in sd), reverse=True)
        r = {s: sd[str(s)] / fv for s in sig}
        ff = FFREE.get(nm)
        ffs = f"{ff:7.4f}" if ff else "      -"
        vals = " ".join(f"{r.get(s, float('nan')):8.4f}" for s in (16.0, 8.0, 4.0, 2.0))
        ext = []
        for a, b in ((16., 8.), (8., 4.), (4., 2.)):
            if a in r and b in r:
                ext.append(f"{2*r[b]-r[a]:8.4f}")
            else:
                ext.append("       -")
        print(f"{nm:>12} {ffs} | {vals} | {' '.join(ext)}")
