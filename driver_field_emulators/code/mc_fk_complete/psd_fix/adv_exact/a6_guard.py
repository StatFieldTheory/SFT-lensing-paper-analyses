"""(i) does the stock recipe converge to R5Exact as its node count grows?
   (ii) does R5Exact's structural guard fire on a cubic-interpolated table?"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bootstrap
ds, bgm, core = _bootstrap.wire()
import corr_op_C_callable as co
import candidates as K

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
L = K.TABLE_LAM
r5 = K.R5Exact(background=bg, apply_c0=False)
cg = float(np.cos(np.deg2rad(1.0/60.0)))

print("=== (i) stock builder vs R5Exact as n_nodes grows ===")
print("    'far' = >5 Mpc from any table node; 'near' = within 5 Mpc")
xs = np.sort(np.random.default_rng(3).uniform(406.2, 2327.8, 400))
d2n = np.min(np.abs(xs[:,None]-L[None,:]), axis=1)
far, near = d2n > 5.0, d2n <= 5.0
ref = np.array([r5.matrix(cg, float(x))[0,0] for x in xs])
print(f"{'n_nodes':>8} {'far max %':>11} {'far med %':>11} {'near max %':>12} {'near med %':>11}")
for n in (160, 320, 640, 1280, 2560):
    b = ds.Sigma2Builder(background=bg, apply_c0=False, n_nodes=n)
    got = np.array([b.matrix(cg, float(x))[0,0] for x in xs])
    e = np.abs(got/ref - 1)
    print(f"{n:>8} {100*e[far].max():>11.4f} {100*np.median(e[far]):>11.2e}"
          f" {100*e[near].max():>12.3f} {100*np.median(e[near]):>11.2e}")

print("\n=== (ii) structural guard ===")
w = r5.verify_structure(cg)
print("  bilinear (production) table : guard passes, max rel dev %.3e" % w)

from canoes.sachs.sft_input.corr_op.table import make_C_fn_lambda_project
cubic = make_C_fn_lambda_project(co._TABLE, interp_method="cubic",
                                 output_h_power=co.OUTPUT_H_POWER)
orig = co.C_fn
co.C_fn = lambda n1, t1, n2, t2: np.asarray(cubic(n1, float(t1), n2, float(t2)), float)
K._corr_op.C_fn = co.C_fn
try:
    r5b = K.R5Exact(background=bg, apply_c0=False)
    r5b.verify_structure(cg)
    print("  CUBIC table                : guard did NOT fire  <-- guard is useless")
except RuntimeError as e:
    print("  CUBIC table                : guard FIRES ->", str(e).split("\n")[0])
finally:
    co.C_fn = orig; K._corr_op.C_fn = orig

print("\n=== (iii) corr_op evaluations, whole production precompute (2 cos) ===")
for name, cls, kw in [("stock", ds.Sigma2Builder, {}),
                      ("R5 EXACT (transpose sym)", K.R5Exact, {}),
                      ("R5 EXACT (no sym shortcut)", K.R5Exact,
                       {"use_transpose_symmetry": False})]:
    cnt = {"n": 0}
    o = co.C_fn
    def counted(n1, t1, n2, t2, _o=o, _c=cnt):
        _c["n"] += 1; return _o(n1, t1, n2, t2)
    co.C_fn = counted; ds._corr_op.C_fn = counted; K._corr_op.C_fn = counted
    b = cls(background=bg, apply_c0=False, **kw)
    b.matrix(1.0, 1000.0); b.matrix(cg, 1000.0)
    co.C_fn = o; ds._corr_op.C_fn = o; K._corr_op.C_fn = o
    print(f"  {name:<28} {cnt['n']:>5d} corr_op calls")
