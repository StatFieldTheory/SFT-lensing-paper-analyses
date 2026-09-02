"""Fixed sigma/dlam scan of the exact expectation under MANY calibrations."""
import sys, time, json
sys.path.insert(0, "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete/review/calib_attack")
import numpy as np
import altcal, fold

gam = float(sys.argv[1])
pairs = sys.argv[2]                       # "500:16,1000:8,1999:4"
names = sys.argv[3].split(",")
tag = sys.argv[4] if len(sys.argv) > 4 else "scan"

fv, exact_node = fold.fold_at(gam)
print(f"gamma = {gam}'   fold = {fv:.6e}  (exact node: {exact_node})", flush=True)
res = {}
for item in pairs.split(","):
    Ns, sigs = item.split(":")
    N, sig = int(Ns), float(sigs)
    t0 = time.time()
    grid, out = altcal.run(gam, N, sig, names)
    print(f"  N={N} dlam={grid.dlam:.4f} sigma={sig} s/dl={sig/grid.dlam:.3f}"
          f"   [{time.time()-t0:.0f} s]", flush=True)
    for nm in names:
        r = out[nm]["tot"] / fv
        res.setdefault(nm, {})[sig] = out[nm]["tot"]
        print(f"      {nm:>12s}  tot={out[nm]['tot']:+.6e}  /fold={r:+.6f}"
              f"  T1share={out[nm]['t1']/out[nm]['tot']:+.4f}", flush=True)
json.dump({"gamma": gam, "fold": fv, "res": {k: {str(s): v for s, v in d.items()}
           for k, d in res.items()}},
          open(f"/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/"
               f"code/mc_fk_complete/review/calib_attack/_{tag}.json", "w"), indent=1)
print("\nlinear sigma->0 extrapolations (2*r(s/2) - r(s)):", flush=True)
for nm in names:
    sigs = sorted(res[nm].keys(), reverse=True)
    line = f"  {nm:>12s}: "
    for a, b in zip(sigs[:-1], sigs[1:]):
        line += f" [{a:g},{b:g}]->{(2*res[nm][b]-res[nm][a])/fv:+.5f}"
    print(line, flush=True)
