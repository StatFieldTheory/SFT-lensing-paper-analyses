import sys, os
from pathlib import Path
_MC2PT = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_sachs_2pt")
_CALLFIX = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware")
for p in (str(_MC2PT), str(_CALLFIX)):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np
import background as bg_mod
import driver_stats as ds
import perm_aware_kappa3_callable as pa

# --- SPEC section 4: bind the corrected vertex callable ------------------
pa.TABLE_PATH = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/table_permclosed.npz")
pa._CACHE = None
ds._k3 = pa
