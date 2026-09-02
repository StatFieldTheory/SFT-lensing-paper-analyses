import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, json
import cleanroom_exact as CR

res = {}
for s in (8.0, 4.0):
    o = CR.run(1.0, 1000, s)
    res[f"sigma{s:g}"] = {
        "local_total": o["local"]["total"],
        "nominal": {k: o["nominal"][k] for k in ("T1_A","T1_B","T2_A","T2_B","T1","T2","total")},
        "smeared": {k: o["smeared"][k] for k in ("T1_A","T1_B","T2_A","T2_B","T1","T2","total")},
        "o0_exact": o["o0_exact"], "o0_local": o["o0_local"], "rho": o["rho"],
    }
    print()
Path(__file__).with_name("cleanroom_results_g1_N1000.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
