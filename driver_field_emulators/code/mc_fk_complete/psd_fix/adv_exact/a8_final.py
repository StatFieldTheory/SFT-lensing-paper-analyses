"""Guard test on an injected non-bilinear C; PSD at extreme gamma; and the
one-line summary numbers."""
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
import sigma2_repaired as R
import candidates as K

bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
cg = float(np.cos(np.deg2rad(1.0/60.0)))

print("=== guard vs an injected NON-bilinear C (control that it can fire) ===")
r5 = K.R5Exact(background=bg, apply_c0=False)
print("  production C          : max rel dev %.3e  -> passes" % r5.verify_structure(cg))
orig = co.C_fn
for amp in (1e-12, 1e-8, 1e-4, 5e-2):
    co.C_fn = lambda n1, t1, n2, t2, _o=orig, _a=amp: (
        np.asarray(_o(n1, t1, n2, t2), float) * (1.0 + _a*np.sin(float(t1)/100.0)))
    K._corr_op.C_fn = co.C_fn
    try:
        w = K.R5Exact(background=bg, apply_c0=False).verify_structure(cg)
        print("  perturbed by %.0e     : max rel dev %.3e  -> PASSES (guard blind)" % (amp, w))
    except RuntimeError as e:
        print("  perturbed by %.0e     : -> FIRES (%s)" % (amp, str(e).split("=")[1].split(">")[0].strip()))
co.C_fn = orig; K._corr_op.C_fn = orig

print("\n=== PSD on the dense 0.48 Mpc grid at extreme gamma ===")
gB = np.linspace(406.0, 2328.0, 4000)
print(f"{'gamma':>8} " + " ".join(f"{n:>22}" for n in ("stock", "R1 knots", "R5 EXACT")))
print(f"{'':>8} " + " ".join(f"{'nviol':>9}{'min/max eig':>13}" for _ in range(3)))
for gam in (0.5, 5.0, 114.3, 596.9):
    c = float(np.cos(np.deg2rad(gam/60.0)))
    cells = []
    for cls in (ds.Sigma2Builder, R.R1Knots, K.R5Exact):
        b = cls(background=bg, apply_c0=False)
        nv = 0; worst = np.inf
        for x in gB:
            M = core._assemble_6x6(b, c, float(x))
            w = np.linalg.eigvalsh(0.5*(M+M.T))
            if w[0] < 0: nv += 1
            worst = min(worst, w[0]/max(abs(w[-1]), 1e-300))
        cells.append(f"{nv:>9d}{worst:>13.3e}")
    print(f"{gam:>8.1f} " + " ".join(cells))

print("\n=== R1 vs R5: how much does the choice between them matter? ===")
xs = np.sort(np.random.default_rng(5).uniform(406.1, 2327.9, 800))
b1 = R.R1Knots(background=bg, apply_c0=False)
for gam in (1.0, 17.3):
    c = float(np.cos(np.deg2rad(gam/60.0)))
    a = np.array([b1.matrix(c, float(x)) for x in xs])
    b = np.array([r5.matrix(c, float(x)) for x in xs])
    sc = np.maximum(np.abs(b).max(axis=(1,2)), 1e-300)
    e = np.abs(a-b).max(axis=(1,2))/sc
    print("  gamma=%5.1f'  max %.4f%%  p99 %.4f%%  median %.2e%%"
          % (gam, 100*e.max(), 100*np.percentile(e,99), 100*np.median(e)))
