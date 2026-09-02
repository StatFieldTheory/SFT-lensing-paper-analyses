"""Clean-room Monte-Carlo extraction of the FK channel.

Written from ``SPEC.md`` alone (plus the module APIs of ``driver_stats`` /
``background`` / ``_bootstrap`` that SPEC.md names).  ``fk_complete_core.py``,
``fk_expect_exact.py`` and ``NOTES.md`` were not read until every number in
``summary.json`` had been produced.

===========================================================================
1.  What has to be measured
===========================================================================
With ``s[k] = resp[k] s[k-1] + F(s[k-1]) dlam + f[k] dlam`` and
``kappa_a = sum_k w_k s_a[k]``, expand in the vertex F:

    s0[k] = resp[k] s0[k-1] + f[k] dlam
    s1[k] = resp[k] s1[k-1] + F(s0[k-1]) dlam

so that, with the discrete window ``Wd[j] = sum_{k>=j} w_k (D[j]/D[k])^2``,

    kappa0_a = sum_j Wd[j] dlam f_a[j]
    kappa1_a = sum_v Wd[v] dlam F_a(s0[v-1]).

The F-linear part of the cross-ray correlator is
``<kappa1_0 kappa0_3> + <kappa0_0 kappa1_3>``.  ``F`` is quadratic, so each
term is a *cubic* functional of the drive; its Gaussian part vanishes with
``<f> = 0``, so the whole thing is already first order in the three-point
cumulant of ``f``.  That is FK.

===========================================================================
2.  Why a naive simulation fails, and the construction used here
===========================================================================
``f = z + 0.5 Q:(zz - C)`` with Q sized so that the field carries the
tabulated zeta.  Q is huge (per-node skewness of order 10), so simulating the
deformed field directly buries the O(Q^1) piece we want under the O(Q^2) and
O(Q^3) pieces we do not.  Make the zeta-order explicit instead of statistical:
introduce a bookkeeping amplitude eps,

    f_eps[k] = z[k] + eps g[k],   g[k] = 0.5 Q[k]:(z[k] z[k]^T - C[k]),

and take the *pathwise* derivative at eps = 0 (one Gaussian path, no finite
differencing).  FK is exactly the eps^1 coefficient.  Because the correlator
is a cubic functional of the drive, the product rule distributes the
deformation over all three legs by itself:

    kappa0_a = Kz_a + eps Kg_a
    kappa1_a = Y_a + eps X_a + O(eps^2),
        Y_a = sum_v Wd[v] dlam F_a(sz[v-1])          (both vertex legs z)
        X_a = sum_v Wd[v] dlam 2 B_a(sz[v-1], sg[v-1])   (one vertex leg g)

with ``B`` the symmetric bilinear form of ``F``.  The estimator

    E = X_0 Kz_3 + Y_0 Kg_3 + X_3 Kz_0 + Y_3 Kg_0

therefore carries the deformation on the free leg (Y*Kg) *and* on either of
the two vertex legs (X*Kz; the factor 2 in ``2 B(sz, sg)`` is exactly those
two placements).  That is all three permutations of
``<fff> = Q:CC + 2 perms``.  Deforming only one leg would return a third.

At eps = 0 the path is purely Gaussian, so E is a quartic polynomial in
Gaussians: bounded variance, and Q enters strictly linearly, so the size of Q
no longer drives the noise.

Variance reduction: ``<Kz> = <Kg> = 0`` exactly (the second because C is the
chain's own covariance, see section 4), so X and Y may be shifted by any
constant.  Leave-one-out sample means are used, which is exactly unbiased.

===========================================================================
3.  The white-noise sum rules of the discrete AR(1) chain
===========================================================================
The chain is ``z[k] = rho z[k-1] + sqrt(1-rho^2) L[k] xi_k``, ``rho =
exp(-dlam/sigma_lambda)``, ``L L^T = V = Sigma2/(2 sigma_lambda)``.  Its
exact node variance and node-to-node covariance are

    P[0] = V[0],  P[k] = rho^2 P[k-1] + (1-rho^2) V[k]
    <z[k] z[j]^T> = rho^{|k-j|} P[min(k,j)]

so the drive's *integrated* covariance density at node k is

    Chat[k] = sum_j dlam <z[k] z[j]^T>
            = dlam ( sum_{j<=k} rho^{k-j} P[j]
                     + P[k] (rho - rho^{N-k}) / (1 - rho) ).

For a slowly varying V this is ``S V`` with ``S = dlam (1+rho)/(1-rho)``,
which tends to ``2 sigma_lambda`` as ``dlam -> 0``: the chain then carries
exactly ``Sigma2`` as a density.  For the *cumulant*, summing the three
Wick placements against smooth kernels gives (kernels vary on ~1000 Mpc,
the chain correlates over sigma_lambda ~ 8 Mpc)

    sum_{legs} kernels x <f f f>_c  ~  sum_k kernels(k) cum3_from_Q(Q[k], Chat[k])

so ``Q[k] = solve_Q(Chat[k], zeta6[k])`` is what makes the field carry the
tabulated three-point cumulant density.

===========================================================================
4.  Two places where SPEC.md is not self-contained, and what is done here
===========================================================================
(a) ``solve_Q(M, Z)``: SPEC.md does not say which M.  In the smooth,
    dlam -> 0 limit ``Chat -> Sigma2`` and every candidate coincides; the
    driver_stats Sigma2 is however NOT smooth (see (b)), so the choice is
    numerically decisive.  ``--q-basis`` selects:
      chat    : M = Chat[k]  -- the chain's own integrated covariance.  DEFAULT.
      sv      : M = S V[k]   -- the exact discrete sum rule for smooth V.
      sigma2  : M = Sigma2[k] = 2 sigma_lambda V[k]  -- the continuum reading.
    ``sv`` and ``sigma2`` differ by the exactly-known factor
    ``(S / 2 sigma_lambda)^2`` on FK: 1.0095 at sigma_lambda = 8 and 1.0380
    at sigma_lambda = 4.

(b) The Sigma2 density is glitchy.  ``Sigma2Builder`` differentiates a
    160-knot cubic spline of ``D^4 C(t,t)``, and the corr_op diagonal has
    small kinks at its own table nodes; the derivative turns them into deep
    V-shaped dips roughly every 130 Mpc, an order of magnitude down over
    ~10 Mpc, and at gamma = 1' two nodes near lambda = 1310 come out
    negative-definite outright.  FK is in exact arithmetic *independent* of
    Sigma2 (Q ~ zeta / M^2 while the estimator ~ Q V^2), but only if M is the
    covariance the chain actually has.  Using ``M = Sigma2`` instead makes
    ``Q ~ 1/Sigma2^2`` blow up in the dips while the chain's own covariance
    stays smoothed by the AR memory, and FK comes out a factor ~5 too large
    at sigma_lambda = 8.  That is why ``chat`` is the default.  ``L[k]`` is
    built by eigen-clipping V to PSD, and ``C[k] = P[k]`` (SPEC.md: "any
    choice shifts the drive by a constant") so that ``<g> = 0`` exactly --
    with ``C = V`` the residual mean drive would leak a Gaussian two-point
    term into the eps^1 coefficient.

===========================================================================
5.  References computed here
===========================================================================
``fk_continuum``  : SPEC.md section 3 quadrature, fine grid.
``fk_discrete``   : the same diagram with the MC's own discrete kernels
                    ``Wd`` and ``Hd[u] = sum_v dlam Wd[v] (D[u]/D[v-1])^4``,
                    i.e. the white-noise (sigma_lambda -> 0) limit of the MC.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

_HERE = Path(__file__).resolve().parent
if str(_HERE.parent) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))

import _bootstrap  # noqa: E402

LAM_MIN = 406.0
LAM_F = 2313.029
ARCMIN = np.pi / (180.0 * 60.0)

# Paper target, SPEC.md section 3 (order 2, a = b = 0).
ANALYTIC_TARGET = {0.5: 1.948e-05, 1.0: 1.598e-05, 2.0: 1.223e-05,
                   5.0: 7.038e-06, 17.3: 1.966e-06}

_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy 1.x / 2.x


# ---------------------------------------------------------------------------
def _F6(s: NDArray[np.float64]) -> NDArray[np.float64]:
    """F(s)_0 = -(s0^2+s1^2+s2^2), F_1 = -2 s0 s1, F_2 = -2 s0 s2, per ray."""
    out = np.empty_like(s)
    for r in (0, 3):
        s0, s1, s2 = s[:, r], s[:, r + 1], s[:, r + 2]
        out[:, r] = -(s0 * s0 + s1 * s1 + s2 * s2)
        out[:, r + 1] = -2.0 * s0 * s1
        out[:, r + 2] = -2.0 * s0 * s2
    return out


def _B6(x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    """Symmetric bilinear form with B(x, x) = F(x)."""
    out = np.empty_like(x)
    for r in (0, 3):
        x0, x1, x2 = x[:, r], x[:, r + 1], x[:, r + 2]
        y0, y1, y2 = y[:, r], y[:, r + 1], y[:, r + 2]
        out[:, r] = -(x0 * y0 + x1 * y1 + x2 * y2)
        out[:, r + 1] = -(x0 * y1 + x1 * y0)
        out[:, r + 2] = -(x0 * y2 + x2 * y0)
    return out


# ---------------------------------------------------------------------------
class Model:
    """Everything that does not depend on the random path."""

    def __init__(self, gamma_arcmin: float, n_nodes: int, sigma_lambda: float,
                 q_basis: str = "chat", c_sub: str = "P", smooth_sigma: float = 0.0,
                 verbose: bool = True):
        self.ds, self.bgmod, _ = _bootstrap.wire()
        self.gamma_arcmin = float(gamma_arcmin)
        self.cos_gamma = float(np.cos(gamma_arcmin * ARCMIN))
        self.N = int(n_nodes)
        self.sigma_lambda = float(sigma_lambda)
        self.q_basis = q_basis
        self.c_sub = c_sub
        self.smooth_sigma = float(smooth_sigma)

        self.bg = self.bgmod.Background(lam_min=120.0, lam_max=LAM_F + 5.0)
        self.dlam = (LAM_F - LAM_MIN) / (self.N - 1)
        self.lam = LAM_MIN + self.dlam * np.arange(self.N)
        self.D = np.asarray(self.bg.D(self.lam), dtype=float)

        resp = np.empty(self.N)
        resp[0] = 1.0
        resp[1:] = (self.D[:-1] / self.D[1:]) ** 2
        self.resp = resp

        w = np.full(self.N, self.dlam)
        w[0] = w[-1] = 0.5 * self.dlam
        self.w = w

        self.rho = float(np.exp(-self.dlam / self.sigma_lambda))
        self.S_chain = self.dlam * (1.0 + self.rho) / (1.0 - self.rho)

        t0 = time.time()
        self._build(verbose)
        if verbose:
            print(f"[model] gamma={gamma_arcmin}' N={self.N} dlam={self.dlam:.6f} "
                  f"sigma={self.sigma_lambda} rho={self.rho:.6f} "
                  f"S={self.S_chain:.5f} (S/2sig)^2="
                  f"{(self.S_chain/(2*self.sigma_lambda))**2:.6f} "
                  f"q_basis={q_basis} C={c_sub}  [{time.time()-t0:.1f}s]",
                  flush=True)

    # -- tables ------------------------------------------------------------
    def _build(self, verbose: bool) -> None:
        ds, N = self.ds, self.N
        Sig = np.empty((N, 6, 6))
        Zet = np.empty((N, 6, 6, 6))
        for k, l in enumerate(self.lam):
            Sig[k] = ds.sigma2_6x6(self.cos_gamma, float(l))
            Zet[k] = ds.zeta6(self.cos_gamma, float(l))
        self.Sigma2, self.zeta6 = Sig, Zet

        if self.smooth_sigma > 0.0:
            # optional robustness knob: Gaussian-smooth the glitchy Sigma2 in
            # lambda.  FK is Sigma2-independent in exact arithmetic, so this
            # must not move the answer if the construction is sound.
            from scipy.ndimage import gaussian_filter1d
            npix = self.smooth_sigma / self.dlam
            Sig = gaussian_filter1d(Sig, npix, axis=0, mode="nearest")
            self.Sigma2_smooth = Sig

        V = Sig / (2.0 * self.sigma_lambda)
        L = np.empty_like(V)
        Vc = np.empty_like(V)          # actual covariance injected per step
        n_clip = 0
        for k in range(N):
            d, U = np.linalg.eigh(0.5 * (V[k] + V[k].T))
            if d.min() < 0:
                n_clip += 1
            dc = np.clip(d, 0.0, None)
            L[k] = U * np.sqrt(dc)
            Vc[k] = (U * dc) @ U.T
        self.V, self.Vc, self.L = V, Vc, L
        self.n_clip = n_clip

        # exact AR node variance P and integrated covariance Chat
        rho = self.rho
        P = np.empty_like(Vc)
        P[0] = Vc[0]
        for k in range(1, N):
            P[k] = rho * rho * P[k - 1] + (1.0 - rho * rho) * Vc[k]
        self.P = P
        acc = np.zeros((6, 6))
        Chat = np.empty_like(Vc)
        for k in range(N):
            acc = rho * acc + P[k]
            fut = (rho - rho ** (N - k)) / (1.0 - rho)
            Chat[k] = self.dlam * (acc + P[k] * fut)
        self.Chat = Chat

        if self.q_basis == "chat":
            M = Chat
        elif self.q_basis == "sv":
            M = self.S_chain * Vc
        elif self.q_basis == "sigma2":
            M = 2.0 * self.sigma_lambda * Vc
        else:
            raise ValueError(f"bad q_basis {self.q_basis!r}")

        Q = np.empty((N, 6, 6, 6))
        dead = 0
        for k in range(N):
            d = np.linalg.eigvalsh(0.5 * (M[k] + M[k].T))
            if d.max() <= 0 or d.min() <= 1e-10 * d.max():
                Q[k] = 0.0          # node cannot carry a cumulant
                dead += 1
                continue
            Q[k] = ds.solve_Q(M[k], Zet[k])
        self.Q = Q
        self.Qflat = Q.reshape(N, 6, 36)
        self.n_dead = dead

        # subtraction constant
        self.Csub = P if self.c_sub == "P" else Vc

        # per-node skewness diagnostic (component 0, ray 0)
        c3 = np.einsum("kamn,kmb,knc->kabc", Q, P, P, optimize=True)
        c3 = (c3 + c3.transpose(0, 2, 1, 3) + c3.transpose(0, 3, 2, 1))[:, 0, 0, 0]
        var = np.clip(P[:, 0, 0], 1e-300, None)
        self.node_skew = c3 / var ** 1.5
        if verbose:
            print(f"[model] V eigen-clipped at {n_clip}/{N} nodes; "
                  f"Q zeroed at {dead}/{N} nodes; per-node skewness comp0: "
                  f"median {np.median(self.node_skew):.3g}, "
                  f"|max| {np.max(np.abs(self.node_skew)):.3g}", flush=True)

    # -- analytic references ----------------------------------------------
    def fk_continuum(self, n: int = 40001) -> float:
        """SPEC.md section 3: FK = 2 int dv W H [-sum_c zeta6_{cc3}]."""
        u = np.linspace(LAM_MIN, LAM_F, n)
        Du = np.asarray(self.bg.D(u), dtype=float)
        du = np.diff(u)
        i2 = 1.0 / Du ** 2
        c2 = np.concatenate(([0.0], np.cumsum(0.5 * (i2[1:] + i2[:-1]) * du)))
        W = Du ** 2 * (c2[-1] - c2)
        i4 = W / Du ** 4
        c4 = np.concatenate(([0.0], np.cumsum(0.5 * (i4[1:] + i4[:-1]) * du)))
        H = Du ** 4 * (c4[-1] - c4)
        # zeta on a coarser grid (smooth in lambda), splined up
        nz = 2001
        lz = np.linspace(LAM_MIN, LAM_F, nz)
        zz = np.empty(nz)
        for i, lv in enumerate(lz):
            Z = self.ds.zeta6(self.cos_gamma, float(lv))
            zz[i] = -(Z[0, 0, 3] + Z[1, 1, 3] + Z[2, 2, 3])
        from scipy.interpolate import CubicSpline
        zc = CubicSpline(lz, zz)(u)
        return float(2.0 * _trapz(W * H * zc, u))

    def fk_discrete(self) -> float:
        """Same diagram with the MC's own discrete kernels (white-noise limit)."""
        D, w, dl = self.D, self.w, self.dlam
        tail = np.cumsum((w / D ** 2)[::-1])[::-1]
        Wd = D ** 2 * tail                       # Wd[j] = sum_{k>=j} w_k (D_j/D_k)^2
        self.Wd = Wd
        # Hd[u] = sum_{v : v-1 >= u} dlam Wd[v] (D[u]/D[v-1])^4
        # vertex node v contributes at leg index i = v - 1
        inv = np.zeros(self.N)
        inv[:-1] = dl * Wd[1:] / D[:-1] ** 4
        tail4 = np.cumsum(inv[::-1])[::-1]
        Hd = D ** 4 * tail4
        self.Hd = Hd
        Z = self.zeta6
        zc = -(Z[:, 0, 0, 3] + Z[:, 1, 1, 3] + Z[:, 2, 2, 3])
        self.fk_integrand = Wd * Hd * zc
        return float(2.0 * np.sum(dl * Wd * Hd * zc))


# ---------------------------------------------------------------------------
def run_batch(m: Model, n_real: int, rng: np.random.Generator) -> dict:
    N, dl, B = m.N, m.dlam, n_real
    a = m.rho
    b = float(np.sqrt(1.0 - a * a))

    sz = np.zeros((B, 6)); sg = np.zeros((B, 6))
    uY = np.zeros((B, 6)); uX = np.zeros((B, 6))
    Kz = np.zeros((B, 6)); Kg = np.zeros((B, 6))
    Yk = np.zeros((B, 6)); Xk = np.zeros((B, 6))
    z = np.zeros((B, 6))

    for k in range(N):
        r = m.resp[k]
        # F is evaluated at the OLD state s[k-1]
        uY *= r
        uY += _F6(sz) * dl
        uX *= r
        uX += (2.0 * _B6(sz, sg)) * dl

        xi = rng.standard_normal((B, 6))
        n = xi @ m.L[k].T
        z = n if k == 0 else a * z + b * n
        g = 0.5 * (((z[:, :, None] * z[:, None, :]) - m.Csub[k]).reshape(B, 36)
                   @ m.Qflat[k].T)

        sz = r * sz + z * dl
        sg = r * sg + g * dl

        wk = m.w[k]
        Kz += wk * sz; Kg += wk * sg
        Yk += wk * uY; Xk += wk * uX

    return {"X0": Xk[:, 0].copy(), "X3": Xk[:, 3].copy(),
            "Y0": Yk[:, 0].copy(), "Y3": Yk[:, 3].copy(),
            "Kz0": Kz[:, 0].copy(), "Kz3": Kz[:, 3].copy(),
            "Kg0": Kg[:, 0].copy(), "Kg3": Kg[:, 3].copy()}


def estimator(d: dict, centre: bool = True) -> NDArray[np.float64]:
    n = d["X0"].size

    def loo(x):
        return (x.sum() - x) / (n - 1) if centre else 0.0

    return ((d["X0"] - loo(d["X0"])) * d["Kz3"]
            + (d["Y0"] - loo(d["Y0"])) * d["Kg3"]
            + (d["X3"] - loo(d["X3"])) * d["Kz0"]
            + (d["Y3"] - loo(d["Y3"])) * d["Kg0"])


def run_seed(m: Model, n_real: int, seed: int, batch_size: int = 2000) -> dict:
    rng = np.random.default_rng(seed)
    parts, done = [], 0
    while done < n_real:
        nb = min(batch_size, n_real - done)
        parts.append(run_batch(m, nb, rng))
        done += nb
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    est = estimator(d, True)
    raw = estimator(d, False)
    n = est.size
    null = d["Kz0"] * d["Kg3"] + d["Kg0"] * d["Kz3"]
    return {"seed": seed, "n": n,
            "fk": float(est.mean()),
            "fk_sem": float(est.std(ddof=1) / np.sqrt(n)),
            "fk_raw": float(raw.mean()),
            "fk_raw_sem": float(raw.std(ddof=1) / np.sqrt(n)),
            "null_F0": float(null.mean()),
            "null_F0_sem": float(null.std(ddof=1) / np.sqrt(n)),
            "O0_cross": float(np.mean(d["Kz0"] * d["Kz3"])),
            "FF_like": float(np.mean(d["Y0"] * d["Kz3"] + d["Y3"] * d["Kz0"]))}


# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument("--N", type=int, default=1000)
    p.add_argument("--sigma", type=float, nargs="+", default=[8.0, 4.0])
    p.add_argument("--n-real", type=int, default=24000)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--seed0", type=int, default=1000)
    p.add_argument("--batch", type=int, default=2000)
    p.add_argument("--q-basis", default="chat", choices=["chat", "sv", "sigma2"])
    p.add_argument("--c-sub", default="P", choices=["P", "V"])
    p.add_argument("--smooth-sigma", type=float, default=0.0)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    out = {"gamma_arcmin": args.gamma, "N": args.N, "n_real": args.n_real,
           "seeds": args.seeds, "q_basis": args.q_basis, "c_sub": args.c_sub,
           "runs": {}}
    tgt = ANALYTIC_TARGET.get(round(args.gamma, 3))
    out["analytic_target"] = tgt

    ref_done = False
    for isig, sig in enumerate(args.sigma):
        m = Model(args.gamma, args.N, sig, q_basis=args.q_basis, c_sub=args.c_sub,
                  smooth_sigma=args.smooth_sigma)
        if not ref_done:
            fkd = m.fk_discrete()
            fkc = m.fk_continuum()
            out["fk_discrete"] = fkd
            out["fk_continuum"] = fkc
            print(f"[ref] SPEC target       = {tgt:.6e}")
            print(f"[ref] continuum quad    = {fkc:.6e}  (/target {fkc/tgt:.4f})")
            print(f"[ref] discrete kernels  = {fkd:.6e}  (/target {fkd/tgt:.4f})",
                  flush=True)
            ref_done = True

        per = []
        for j in range(args.seeds):
            t0 = time.time()
            r = run_seed(m, args.n_real, args.seed0 + 100 * isig + j, args.batch)
            r["t"] = time.time() - t0
            per.append(r)
            print(f"  sig={sig} seed={r['seed']} fk={r['fk']:.5e} "
                  f"(+/-{r['fk_sem']:.1e}) ratio={r['fk']/tgt:.4f} "
                  f"raw={r['fk_raw']:.4e}(+/-{r['fk_raw_sem']:.1e}) "
                  f"null={r['null_F0']:.2e}(+/-{r['null_F0_sem']:.1e}) "
                  f"[{r['t']:.0f}s]", flush=True)

        v = np.array([r["fk"] for r in per])
        mean, sem = float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))
        out["runs"][str(sig)] = {
            "per_seed": per, "mean": mean, "sem": sem,
            "ratio_target": mean / tgt, "ratio_target_sem": sem / tgt,
            "ratio_discrete": mean / out["fk_discrete"],
            "S_over_2sig_sq": (m.S_chain / (2 * sig)) ** 2,
            "O0_cross": float(np.mean([r["O0_cross"] for r in per])),
            "FF_like": float(np.mean([r["FF_like"] for r in per])),
            "n_clip": m.n_clip, "n_dead": m.n_dead}
        print(f"[sigma={sig}]  FK = {mean:.6e} +/- {sem:.2e} (seed SEM over "
              f"{args.seeds} seeds x {args.n_real})")
        print(f"              ratio to SPEC target {tgt:.4e} = "
              f"{mean/tgt:.4f} +/- {sem/tgt:.4f}")
        print(f"              ratio to discrete kernels        = "
              f"{mean/out['fk_discrete']:.4f} +/- {sem/out['fk_discrete']:.4f}",
              flush=True)

    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1, default=float))
        print("wrote", args.out)


if __name__ == "__main__":
    main()
