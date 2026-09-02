"""Sanity probe: node covariance, deformation tensor, and placement share."""
from __future__ import annotations
import numpy as np
import _bootstrap

ds, bgm, core = _bootstrap.wire()

gam = 1.0
cosg = float(np.cos(np.deg2rad(gam / 60.0)))
sig = 8.0
lam = 1200.0
bg = bgm.Background(lam_min=120.0, lam_max=bgm.LAM_SOURCE_BASELINE + 5.0)
builder = ds.Sigma2Builder(background=bg, apply_c0=False)

S2 = core._assemble_6x6(builder, cosg, lam)
V = S2 / (2.0 * sig)
Z6 = ds.zeta6(cosg, lam)
Zt = Z6 / (2.0 * sig) ** 2
Q = ds.solve_Q(V, Zt)
print("Sigma2 diag :", np.diag(S2))
print("V diag      :", np.diag(V))
print("Z6[0,0,3]   :", Z6[0, 0, 3], " Z6[0,0,0]:", Z6[0, 0, 0])
print("Q max/min   :", Q.max(), Q.min())
print("node skew   :", Zt[0, 0, 0] / V[0, 0] ** 1.5)
rt = ds.cum3_from_Q(Q, V)
print("roundtrip rel err:", np.abs(rt - Zt).max() / np.abs(Zt).max())

# placement share in the V-eigenbasis: (Q:VV)_pqr / zeta_pqr = d_q d_r / sum
d, U = np.linalg.eigh(V)
print("eigs d      :", d)
P = np.einsum("amn,mb,nc->abc", Q, V, V, optimize=True)     # Q on leg a only
print("share(0,0,3) P/Z:", P[0, 0, 3] / Zt[0, 0, 3], " P[3,0,0]/Z:", P[3, 0, 0] / Zt[0, 0, 3])
print("sum of 3 perms / Z:", (P[0,0,3] + P[0,0,3].__class__(P.transpose(1,0,2)[0,0,3]) + P.transpose(2,1,0)[0,0,3]) / Zt[0,0,3])
