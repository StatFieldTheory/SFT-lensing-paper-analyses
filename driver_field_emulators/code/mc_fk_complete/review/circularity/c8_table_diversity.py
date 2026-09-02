"""How different ARE the tables used in c3?  (guard against a fake replication)"""
import _env, importlib.util, sys, numpy as np
from pathlib import Path
import fk_expect_exact as ex
from _fold import fold_nodes

PROD = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products")
CALS = {
 "permclosed_cut15360": "_variant_table_permclosed_cut15360_permfix/callable__variant_table_permclosed_cut15360_permfix.py",
 "permclosed_cut1000": "_variant_table_permclosed_cut1000_permfix/callable__variant_table_permclosed_cut1000_permfix.py",
 "permclosed_cut1000_r1": "_variant_table_permclosed_cut1000_r1_permfix/callable__variant_table_permclosed_cut1000_r1_permfix.py",
 "permclosed_cut1000_LEGACY": "_variant_table_permclosed_cut1000_legacy/callable__variant_table_permclosed_cut1000_legacy.py",
 "tree_cut1000_r4": "_variant_table_tree_cut1000_r4/callable__variant_table_tree_cut1000_r4.py",
 "tree_cut15360_r4": "_variant_table_tree_cut15360_r4/callable__variant_table_tree_cut15360_r4.py",
 "bihalofit_cut960_r4": "_variant_table_bihalofit_cut960_r4/callable__variant_table_bihalofit_cut960_r4.py",
 "bihalofit_cut15360_r4": "_variant_table_bihalofit_cut15360_r4/callable__variant_table_bihalofit_cut15360_r4.py",
}
lam = np.linspace(406.0, 2313.029, 200)
cosg = float(np.cos(np.deg2rad(1.0155 / 60.0)))
prof = {}
for i, (nm, cal) in enumerate(CALS.items()):
    spec = importlib.util.spec_from_file_location(f"_tk3_{i}", PROD / cal)
    m = importlib.util.module_from_spec(spec); sys.modules[f"_tk3_{i}"] = m
    spec.loader.exec_module(m)
    ex._DS._k3 = m
    z = np.array([ex._DS.zeta6(cosg, float(l)) for l in lam])
    prof[nm] = -np.einsum("Abc,mbcB->m", ex.F6, z)[...] if False else \
        np.einsum("Abc,mbcB->mAB", ex.F6, z, optimize=True)[:, 0, 3]
ref = prof["permclosed_cut15360"]
print(f"{'table':<28}{'|F:zeta| at 1200 Mpc':>22}{'shape corr w/ prod':>20}"
      f"{'peak lam':>10}")
for nm, p in prof.items():
    s = p / np.abs(p).max()
    r = ref / np.abs(ref).max()
    print(f"{nm:<28}{p[np.argmin(abs(lam-1200))]:>22.4e}"
          f"{float(np.corrcoef(s, r)[0,1]):>20.6f}{lam[np.argmax(abs(p))]:>10.0f}")
print("\npairwise max |shape ratio - 1| against the production table:")
for nm, p in prof.items():
    q = (p / p.sum()) / (ref / ref.sum())
    good = np.abs(ref) > 0.02 * np.abs(ref).max()
    print(f"  {nm:<28}{np.abs(q[good]-1).max():>10.4f}")
