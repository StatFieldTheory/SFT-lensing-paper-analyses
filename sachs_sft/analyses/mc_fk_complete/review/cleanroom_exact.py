"""Exact (Wick) evaluation of the same FK estimator the Monte-Carlo samples.

Same model object as ``cleanroom_mc``.  This is not an independent physics
derivation -- it is the noise-free value of

    <E> = <X_0 Kz_3> + <Y_0 Kg_3> + <X_3 Kz_0> + <Y_3 Kg_0>,

so it (a) verifies that the MC code implements the estimator it claims to,
and (b) gives a zero-variance handle on the sigma_lambda dependence.

Ingredients (all with the MC's own discrete kernels):

    Wd[j]      = sum_{k>=j} w_k (D[j]/D[k])^2
    R[k,i]     = (D[i]/D[k])^2  Theta(k-i)
    C(i,l)     = <z[i] z[l]^T> = rho^|i-l| P[min(i,l)]
    A_cm(k,l)  = <sz_c[k] z_m[l]> = sum_{i<=k} R[k,i] dlam C_cm(i,l)
                 (recursion  A(k,.) = resp[k] A(k-1,.) + dlam C(k,.))
    Chatf_n(i) = sum_l Wd[l] dlam C_{n,f}(i,l)

With C[k] = P[k] the subtraction is exact, so <g> = 0 and the pair (z_m z_n)
never contracts with itself.  Wick then gives, for vertex ray r (components
c in {3r, 3r+1, 3r+2}) and free leg f = 3(1-r):

    <Y_{3r} Kg_f> = - sum_{k,l} Wd[k+1] Wd[l] dlam^2
                      Q_{f,mn}(l) sum_c A_cm(k,l) A_cn(k,l)

    <X_{3r} Kz_f> = -2 sum_{k,i} Wd[k+1] dlam R[k,i] dlam
                      sum_c Q_{c,mn}(i) A_cm(k,i) Chatf_n(i)

(k = v-1 is the state index the Euler step evaluates F at.)

NOTE: the Wick reduction above assumes C[k] = P[k], i.e. that the subtraction
constant is the chain's own covariance so the (z_m z_n) self-pairing cancels.
It is therefore only valid for ``--c-sub P``; the ``V`` case must be measured
with the MC.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import cleanroom_mc as cm  # noqa: E402


def fk_exact(m: cm.Model, verbose: bool = True) -> dict:
    N, dl = m.N, m.dlam
    D, w, rho = m.D, m.w, m.rho
    tail = np.cumsum((w / D ** 2)[::-1])[::-1]
    Wd = D ** 2 * tail

    P = m.P
    idx = np.arange(N)

    # Chatf_n(i) = sum_l Wd[l] dlam C_{n,f}(i,l), for f = 0 and f = 3
    Chat_f = np.zeros((2, N, 6))
    for i in range(N):
        rr = rho ** np.abs(idx - i)
        mn = np.minimum(idx, i)
        Ci = rr[:, None, None] * P[mn]              # (N,6,6) = C(i,l) over l
        wl = Wd * dl
        for fi, f in enumerate((0, 3)):
            Chat_f[fi, i] = np.einsum("l,ln->n", wl, Ci[:, :, f])

    Q = m.Q
    # accumulators, index 0 -> vertex ray 0 (free leg 3), 1 -> vertex ray 1 (leg 0)
    Yterm = np.zeros(2)
    Xterm = np.zeros(2)
    T = np.zeros((2, N, 6, 6))                      # T_mn(l) per vertex ray

    A = np.zeros((N, 6, 6))                         # A[l, c, m] at current k
    D2 = D ** 2
    t0 = time.time()
    for k in range(N):
        rr = rho ** np.abs(idx - k)
        mn = np.minimum(idx, k)
        Ck = rr[:, None, None] * P[mn]              # C(k,l) over l, (N,6,6)
        A = m.resp[k] * A + dl * Ck                 # A_cm(k,l)
        if k + 1 >= N:
            break
        wk = Wd[k + 1] * dl
        Rki = np.zeros(N)
        Rki[:k + 1] = D2[:k + 1] / D2[k]            # R[k,i], i<=k

        for r in (0, 1):
            cs = slice(3 * r, 3 * r + 3)
            # Y: T_mn(l) += wk * sum_c A_cm A_cn
            T[r] += wk * np.einsum("lcm,lcn->lmn", A[:, cs, :], A[:, cs, :])
            # X: sum_i R[k,i] dlam sum_c Q_cmn(i) A_cm(k,i) Chatf_n(i)
            f_other = 1 - r
            pc = np.einsum("icmn,icm,in->i", Q[:, cs, :, :], A[:, cs, :],
                           Chat_f[f_other], optimize=True)
            Xterm[r] += -2.0 * wk * dl * float(np.dot(Rki, pc))

    for r in (0, 1):
        f = 3 * (1 - r)
        Yterm[r] = -float(np.einsum("l,lmn,lmn->", Wd * dl, Q[:, f, :, :], T[r]))

    if verbose:
        print(f"[exact] {time.time()-t0:.1f} s   X(ray0)={Xterm[0]:.5e} "
              f"Y(ray0)={Yterm[0]:.5e}  X(ray1)={Xterm[1]:.5e} "
              f"Y(ray1)={Yterm[1]:.5e}", flush=True)
    total = float(Xterm.sum() + Yterm.sum())
    return {"fk_exact": total, "X": Xterm.tolist(), "Y": Yterm.tolist()}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument("--N", type=int, default=1000)
    p.add_argument("--sigma", type=float, nargs="+", default=[8.0, 4.0])
    p.add_argument("--q-basis", default="chat", choices=["chat", "sv", "sigma2"])
    p.add_argument("--c-sub", default="P", choices=["P", "V"])
    p.add_argument("--smooth-sigma", type=float, default=0.0)
    p.add_argument("--out", default=None)
    a = p.parse_args()

    tgt = cm.ANALYTIC_TARGET.get(round(a.gamma, 3))
    res = {"gamma": a.gamma, "N": a.N, "q_basis": a.q_basis, "target": tgt,
           "runs": {}}
    ref = False
    for sig in a.sigma:
        m = cm.Model(a.gamma, a.N, sig, q_basis=a.q_basis, c_sub=a.c_sub,
                     smooth_sigma=a.smooth_sigma)
        if not ref:
            res["fk_discrete"] = m.fk_discrete()
            res["fk_continuum"] = m.fk_continuum()
            print(f"[ref] target {tgt:.6e} | continuum {res['fk_continuum']:.6e} "
                  f"| discrete {res['fk_discrete']:.6e}", flush=True)
            ref = True
        else:
            m.fk_discrete()
        r = fk_exact(m)
        r["ratio_target"] = r["fk_exact"] / tgt
        r["ratio_discrete"] = r["fk_exact"] / res["fk_discrete"]
        res["runs"][str(sig)] = r
        print(f"[exact] sigma={sig}: FK = {r['fk_exact']:.6e}  "
              f"ratio/target = {r['ratio_target']:.4f}  "
              f"ratio/discrete = {r['ratio_discrete']:.4f}", flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1, default=float))
        print("wrote", a.out)


if __name__ == "__main__":
    main()
