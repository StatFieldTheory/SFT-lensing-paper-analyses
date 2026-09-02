"""Clean-room discrete FK reference, derived from SPEC.md section 1 alone.

DERIVATION (nothing below was taken from fk_expect_exact.py).

Grid  lam[k] = lam_min + k dlam,  k = 0..N-1,  lam[N-1] = lam_f.
State s[k] in R^6, s[-1] = 0,

    s[k] = resp[k] s[k-1] + F(s[k-1]) dlam + f[k] dlam,
    resp[0] = 1,  resp[k] = (D[k-1]/D[k])^2.

The homogeneous propagator telescopes exactly:

    R[k,j] = prod_{i=j+1}^{k} resp[i] = (D[j]/D[k])^2      (k >= j >= 0).

Order 0 in F (linear response to the drive):

    s^(1)[k] = dlam sum_{j<=k} (D[j]/D[k])^2 f[j],
    kappa^(1)_a = sum_k w_k s^(1)_a[k] = dlam sum_j Wd[j] f_a[j],

with the DISCRETE OBSERVABLE WINDOW

    Wd[j] = sum_{k>=j} w_k (D[j]/D[k])^2,   w_k = dlam, w_0 = w_{N-1} = dlam/2.

Order 1 in F.  F is evaluated at the OLD state, so the vertex at step m eats
s^(1)[m-1]:

    s^(2)[k] = dlam sum_{m<=k} (D[m]/D[k])^2 F(s^(1)[m-1]),
    kappa^(2)_a = dlam sum_m Wd[m] F_{aBC} s^(1)_B[m-1] s^(1)_C[m-1]
                = dlam^3 sum_m Wd[m] F_{aBC}
                  sum_{p,q<=m-1} (D[p]/D[m-1])^2 (D[q]/D[m-1])^2 f_B[p] f_C[q].

Local (white-noise) limit of the drive.  f is a DENSITY: delta(lam-lam') ->
delta_{kl}/dlam, so

    <f_B[p] f_C[q] f_b[j]>_c = zeta6_{BCb}[p] delta_{pq} delta_{pj} / dlam^2

(check: sum_{q,j} dlam^2 <f f f> = zeta6[p], the tabulated density).

Hence, keeping first order in F and first order in zeta,

    <kappa^(2)_a kappa^(1)_b>_c
      = dlam^4 sum_m Wd[m] F_{aBC} sum_{p<=m-1} (D[p]/D[m-1])^4 Wd[p]
        zeta_{BCb}[p] / dlam^2
      = sum_p dlam Wd[p] Hd[p] F_{aBC} zeta_{BCb}[p],

with the DISCRETE F-LEG RESPONSE INTEGRAL

    Hd[p] = dlam sum_{m=p+1}^{N-1} Wd[m] (D[p]/D[m-1])^4.

Adding the mirror term <kappa^(1) kappa^(2)> (the two ray orderings) gives

    FK_ref = 2 sum_p dlam Wd[p] Hd[p] F_{Abc} zeta_{bcB}[p].

For the kappa-kappa entry: A = 0 (Phi00 of ray 1) so F_{0bc} = -delta_{bc}
over the ray-1 components b,c in {0,1,2}; B = 3 (Phi00 of ray 2).  So the
contraction is  -sum_{c=0}^{2} zeta6_{c,c,3},  which is a single ray-triple
(n1,n1,n2) evaluation of the (3,3,3) vertex, entry [c,c,0].

Continuum limits: Wd -> W(la) = int_la^{lam_f} [D(la)/D(t)]^2 dt,
Hd -> H(v) = D(v)^4 int_v^{lam_f} W(u)/D(u)^4 du, sum_p dlam -> int dv.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _bootstrap import wire, FOLD  # noqa: E402

LAM_MIN = 406.0


def discrete_kernels(bg, N: int, lam_f: float):
    """Wd, Hd, dlam, lam on the SPEC grid.  Pure algebra, no quadrature rule."""
    lam = np.linspace(LAM_MIN, lam_f, int(N))
    dlam = float(lam[1] - lam[0])
    D = np.asarray(bg.D(lam), dtype=float)

    w = np.full(N, dlam)
    w[0] = w[-1] = 0.5 * dlam

    # Wd[j] = D[j]^2 * sum_{k>=j} w[k]/D[k]^2
    tail = np.cumsum((w / D**2)[::-1])[::-1]
    Wd = D**2 * tail

    # Hd[p] = dlam * D[p]^4 * sum_{m=p+1}^{N-1} Wd[m]/D[m-1]^4
    G = np.zeros(N)                      # G[m] for m = 1..N-1
    G[1:] = Wd[1:] / D[:-1] ** 4
    tail4 = np.zeros(N + 1)
    tail4[:-1] = np.cumsum(G[::-1])[::-1]   # tail4[p] = sum_{m>=p} G[m]
    Hd = dlam * D**4 * tail4[1:N + 1][:N]   # sum_{m>=p+1}
    return lam, dlam, Wd, Hd


def zeta_contraction(pa, cos_gamma: float, lam: np.ndarray) -> np.ndarray:
    """-sum_c zeta6_{c,c,3}(cos_gamma; lam), vectorised over lam."""
    n1 = np.array([0.0, 0.0, 1.0])
    c = float(np.clip(cos_gamma, -1.0, 1.0))
    n2 = np.array([math.sqrt(max(0.0, 1.0 - c * c)), 0.0, c])
    n = len(lam)
    n_arr = np.stack([np.repeat(n1[None, :], n, 0),
                      np.repeat(n1[None, :], n, 0),
                      np.repeat(n2[None, :], n, 0)], axis=0)
    t_arr = np.stack([lam, lam, lam], axis=0)
    z = pa.coupling_fn_batch(n_arr, t_arr)          # (n,3,3,3)
    return -(z[:, 0, 0, 0] + z[:, 1, 1, 0] + z[:, 2, 2, 0])


def fk_ref(bg, pa, gamma_arcmin: float, N: int, lam_f: float) -> float:
    lam, dlam, Wd, Hd = discrete_kernels(bg, N, lam_f)
    tru = zeta_contraction(pa, math.cos(math.radians(gamma_arcmin / 60.0)), lam)
    return 2.0 * dlam * float(np.sum(Wd * Hd * tru))


def load_fold_kk(path: Path):
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return g[o], v[o]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 5.0, 17.3])
    ap.add_argument("--Ns", type=int, nargs="+",
                    default=[250, 500, 1000, 2000, 4000])
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    ds, bgm, _core = wire()
    import perm_aware_kappa3_callable as pa
    lam_f = bgm.LAM_SOURCE_BASELINE
    bg = bgm.Background(lam_min=120.0, lam_max=lam_f + 5.0)

    g_fold, v_fold = load_fold_kk(FOLD)
    print(f"[ref] lam in [{LAM_MIN}, {lam_f}]")
    hdr = "  ".join(f"{'N=%d' % N:>13}" for N in a.Ns)
    print(f"{'gamma':>7} {hdr}   {'fold':>13}")
    rows = []
    for g in a.gammas:
        vals = [fk_ref(bg, pa, g, N, lam_f) for N in a.Ns]
        fold = float(np.interp(g, g_fold, v_fold))
        rows.append([g] + vals + [fold])
        print(f"{g:7.2f} " + "  ".join(f"{v:13.6e}" for v in vals)
              + f"   {fold:13.6e}")
    print()
    print("ratio to fold")
    print(f"{'gamma':>7} " + "  ".join(f"{'N=%d' % N:>13}" for N in a.Ns))
    for r in rows:
        print(f"{r[0]:7.2f} " + "  ".join(f"{v / r[-1]:13.6f}" for v in r[1:-1]))
    if a.out:
        np.savez(a.out, gamma=np.array([r[0] for r in rows]),
                 Ns=np.array(a.Ns),
                 vals=np.array([r[1:-1] for r in rows]),
                 fold=np.array([r[-1] for r in rows]))
        print(f"[ref] -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
