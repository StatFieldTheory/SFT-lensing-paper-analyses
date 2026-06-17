"""NumPy Monte-Carlo engine for the two-ray Sachs fluctuation.

For a fixed separation ``gamma`` we evolve the Sachs fluctuation 3-vector
``s = (s_0, s_1, s_2) = (kappa-rate, gamma_+-rate, gamma_x-rate)`` along two rays
``n1, n2`` (``n1.n2 = cos gamma``) under a shared, correlated driving field, and
read the convergence/shear as the line-of-sight integral ``kappa_a = int s_a
dlambda`` (the ``integrate_over='all'`` observable; see DESIGN.md sec. 2-3).

Dynamics (exponential integrator; the linear part is the closed-form response):

    s_{k+1} = [D(lam_k)/D(lam_{k+1})]^2 s_k + (F_abc s_b s_c) dlam + dW_k
    kappa_a = sum_k s_a(lam_k) dlam

with dW_k a correlated 6-vector increment (3 components x 2 rays) of covariance
``Sigma2_6x6(cos gamma; lam_k) * dlam`` and, in the skewed run, a third cumulant
``zeta`` injected by a local-quadratic deformation (demo2 generalization).

The linear part telescopes analytically to ``<kappa kappa> = int Sigma2 G^2 dlam``
(driver_stats.order0_mc), so the Gaussian run reproduces the analytic Order-0 by
construction.  The F-vertex adds the FF channel; the skewness adds FK.

Engine mirrors sft-wick/examples/demo2/run_simulation.py (batched realizations,
common-random-number alpha-toggle, blow-up guard).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

import numpy as np
from numpy.typing import NDArray

import background as _bg
import driver_stats as _ds

# --- paper F-vertex (sections/sachs_dynamics.tex, Table 1; 0-indexed) --------
# F_111=-1, F_122=-1, F_133=-1, F_212=-2, F_313=-2 (others 0).
_F: Final[NDArray[np.float64]] = np.zeros((3, 3, 3))
_F[0, 0, 0] = -1.0
_F[0, 1, 1] = -1.0
_F[0, 2, 2] = -1.0
_F[1, 0, 1] = -2.0
_F[2, 0, 2] = -2.0


def _f_vertex(s: NDArray[np.float64]) -> NDArray[np.float64]:
    """F_abc s_b s_c for s shaped (3, ...). Returns (3, ...)."""
    s0, s1, s2 = s[0], s[1], s[2]
    out = np.empty_like(s)
    out[0] = -(s0 * s0 + s1 * s1 + s2 * s2)
    out[1] = -2.0 * s0 * s1
    out[2] = -2.0 * s0 * s2
    return out


@dataclass(frozen=True)
class MCConfig:
    """Monte-Carlo run configuration."""
    n_real: int = 50_000
    batch_size: int = 2_000
    n_lambda: int = 600              # lambda grid nodes (observer -> source)
    lam_min: float = 406.0           # Mpc; = Sigma2Builder.lam_lo (corr_op support)
    lam_source: float = _bg.LAM_SOURCE_BASELINE  # 2313.029, matches analysis-3
    seed: int = 12345
    use_f_vertex: bool = True        # include the F-vertex (FF channel)
    apply_anchor: bool = False       # honest equal-time Sigma2 (no 1.135 fudge)
    skew_scale: float = 0.0          # 0 = Gaussian; 1 = full zeta injection
    sigma_lambda: float = 0.0        # 0 = white noise; >0 = AR(1) corr length (Mpc)
    blowup: float = 1e6


@dataclass
class MCResult:
    gamma_arcmin: float
    cos_gamma: float
    xi: NDArray[np.float64]          # (3,3) <kappa_a(n1) kappa_b(n2)>
    xi_err: NDArray[np.float64]      # (3,3) standard error of the mean
    n_good: int
    o0_analytic: float              # driver_stats.order0_mc (kk), same Sigma2 mode


def _cholesky_psd(M: NDArray[np.float64]) -> NDArray[np.float64]:
    """Cholesky with a tiny eigen-floor fallback for marginal PSD matrices."""
    try:
        return np.linalg.cholesky(M)
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(M)
        w = np.maximum(w, 0.0)
        return V * np.sqrt(w)[None, :]


def _assemble_6x6(builder: "_ds.Sigma2Builder", cos_gamma: float,
                  lam: float) -> NDArray[np.float64]:
    """Two-ray (6,6) covariance from a specific (anchor-aware) builder."""
    within = builder.matrix(1.0, lam)
    cross = builder.matrix(float(cos_gamma), lam)
    out = np.zeros((6, 6))
    out[0:3, 0:3] = within
    out[3:6, 3:6] = within
    out[0:3, 3:6] = cross
    out[3:6, 0:3] = cross.T
    return out


def _precompute(cfg: MCConfig, cos_gamma: float):
    """Lambda grid, response ratios, and per-node noise factors (6-dim)."""
    bg = _bg.Background(lam_min=120.0, lam_max=cfg.lam_source + 5.0)
    lam = np.linspace(cfg.lam_min, cfg.lam_source, cfg.n_lambda)
    dlam = float(lam[1] - lam[0])
    D = bg.D(lam)
    # Exact discrete propagator: response INTO node k from node k-1,
    # resp[k] = [D(lam_{k-1})/D(lam_k)]^2, resp[0]=1 (s_{-1}=0).  Injecting the
    # noise at node k with unit response then yields, telescoping,
    #   s_k = sum_{j<=k} [D(lam_j)/D(lam_k)]^2 dW_j   (exact, no O(dlam) bias).
    resp = np.empty(cfg.n_lambda)
    resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    # per-node 6x6 white-noise increment covariance and its Cholesky
    builder = _ds.Sigma2Builder(background=bg, apply_c0=cfg.apply_anchor)
    L_white = np.empty((cfg.n_lambda, 6, 6))
    Q_node = np.zeros((cfg.n_lambda, 6, 6, 6))
    sig = max(cfg.sigma_lambda, 0.0)
    for k, lk in enumerate(lam):
        S2 = _assemble_6x6(builder, float(cos_gamma), float(lk))  # density (per dlam)
        if sig <= 0.0:
            M2 = S2 * dlam                                 # white-noise increment cov
        else:
            # AR(1) stationary node covariance V = Sigma2/(2 sigma_lambda)
            M2 = S2 / (2.0 * sig)
        L_white[k] = _cholesky_psd(M2)
        if cfg.skew_scale != 0.0:
            Z6 = _ds.zeta6(float(cos_gamma), float(lk))
            # node 3rd-cumulant target consistent with M2 (see DESIGN sec.3)
            if sig <= 0.0:
                Z_target = Z6 * dlam
            else:
                Z_target = Z6 / (2.0 * sig)
            Q_node[k] = cfg.skew_scale * _ds.solve_Q(M2, Z_target)
    return bg, lam, dlam, resp, L_white, Q_node, sig, builder


def _draw_increment(rng, Lk, Qk, m: int, skew: bool):
    """Draw m correlated 6-vectors with cov Lk Lk^T and optional skewness."""
    z = rng.standard_normal((6, m))
    dW = Lk @ z                                  # (6, m), Gaussian
    if skew:
        # f_tilde = f + 0.5 Q (f x f - <f x f>); <f x f> = Lk Lk^T
        M2 = Lk @ Lk.T
        ff = np.einsum("im,jm->ijm", dW, dW) - M2[:, :, None]
        dW = dW + 0.5 * np.einsum("aij,ijm->am", Qk, ff)
    return dW                                    # (6, m)


def simulate(cfg: MCConfig, gamma_arcmin: float) -> MCResult:
    """Run the two-ray MC at one separation; return <kappa_a(n1) kappa_b(n2)>."""
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, L_white, Q_node, sig, builder = _precompute(cfg, cos_gamma)
    rng = np.random.default_rng(cfg.seed)
    skew = cfg.skew_scale != 0.0
    colored = sig > 0.0
    rho = np.exp(-dlam / sig) if colored else 0.0

    n_batches = max(1, cfg.n_real // cfg.batch_size)
    # accumulators for <kappa_a(n1) kappa_b(n2)> over the 6-dim kappa vector
    acc = np.zeros((6, 6))
    acc_sq = np.zeros((6, 6))
    n_good = 0

    for _ in range(n_batches):
        m = cfg.batch_size
        s = np.zeros((6, m))               # fluctuation (ray1[3], ray2[3])
        kap = np.zeros((6, m))             # accumulated kappa = sum s dlam
        f_prev = np.zeros((6, m))          # AR(1) field state (colored only)
        blown = np.zeros(m, dtype=bool)
        for k in range(cfg.n_lambda):
            if colored:
                # AR(1) update of the driving field, then deform for skewness
                innov = _draw_increment(rng, L_white[k] * np.sqrt(1.0 - rho * rho),
                                        Q_node[k], m, skew)
                f = rho * f_prev + innov if k > 0 else _draw_increment(
                    rng, L_white[k], Q_node[k], m, skew)
                f_prev = f
                noise = f * dlam            # contribution to s over the step
            else:
                noise = _draw_increment(rng, L_white[k], Q_node[k], m, skew)
            # exact-propagator update: response INTO node k, then inject noise
            # (and the F-vertex drift) at node k with unit response.
            # s laid out (ray0[3], ray1[3]); F-vertex acts per ray on its 3-vector.
            if cfg.use_f_vertex:
                s_cr = s.reshape(2, 3, m).transpose(1, 0, 2)   # (comp, ray, m)
                fv = _f_vertex(s_cr).transpose(1, 0, 2).reshape(6, m)
                s = resp[k] * s + fv * dlam + noise
            else:
                s = resp[k] * s + noise
            # trapezoidal line-of-sight integral kappa = int s dlam
            w = 0.5 * dlam if (k == 0 or k == cfg.n_lambda - 1) else dlam
            kap += w * s
            bad = np.any(np.abs(s) > cfg.blowup, axis=0)
            blown |= bad
        good = ~blown
        kg = kap[:, good]
        acc += kg @ kg.T
        acc_sq += (kg * kg) @ (kg * kg).T  # for a crude SE estimate
        n_good += int(good.sum())

    xi6 = acc / max(n_good, 1)
    # standard error of the mean per (a,b): sqrt(var(prod)/n)
    var6 = acc_sq / max(n_good, 1) - xi6 * xi6
    se6 = np.sqrt(np.maximum(var6, 0.0) / max(n_good, 1))
    # cross-ray block: kappa_a(n1) = kap[a], kappa_b(n2) = kap[3+b]
    xi = xi6[0:3, 3:6]
    xi_err = se6[0:3, 3:6]
    o0 = float(_ds.order0_mc(cos_gamma, lam_f=cfg.lam_source, builder=builder,
                             n_gauss=256))
    return MCResult(gamma_arcmin, cos_gamma, xi, xi_err, n_good, o0)


@dataclass
class CRNResult:
    gamma_arcmin: float
    cos_gamma: float
    sigma_lambda: float
    xi_gauss: NDArray[np.float64]    # (3,3) <kappa_a(n1) kappa_b(n2)> Gaussian arm
    xi_skew: NDArray[np.float64]     # (3,3) skewed arm (O0+FF+FK+...)
    xi_fk: NDArray[np.float64]       # (3,3) CRN difference skew - gauss (= FK + ...)
    xi_fk_err: NDArray[np.float64]   # (3,3) SE of the CRN difference (low variance)
    n_good: int


def _field_precompute(cfg: MCConfig, cos_gamma: float):
    """Colored-noise node stats: V=Sigma2/(2 sigma), its Cholesky, and Q from zeta."""
    if cfg.sigma_lambda <= 0.0:
        raise ValueError("CRN/skewed run requires sigma_lambda > 0 (colored noise)")
    bg = _bg.Background(lam_min=120.0, lam_max=cfg.lam_source + 5.0)
    lam = np.linspace(cfg.lam_min, cfg.lam_source, cfg.n_lambda)
    dlam = float(lam[1] - lam[0])
    D = bg.D(lam)
    resp = np.empty(cfg.n_lambda)
    resp[0] = 1.0
    resp[1:] = (D[:-1] / D[1:]) ** 2
    builder = _ds.Sigma2Builder(background=bg, apply_c0=cfg.apply_anchor)
    sig = cfg.sigma_lambda
    V = np.empty((cfg.n_lambda, 6, 6))         # node stationary covariance
    Lf = np.empty((cfg.n_lambda, 6, 6))        # chol(V)
    Q = np.zeros((cfg.n_lambda, 6, 6, 6))      # deformation tensor
    for k, lk in enumerate(lam):
        S2 = _assemble_6x6(builder, float(cos_gamma), float(lk))
        Vk = S2 / (2.0 * sig)                  # one-delta collapse: ∫k dΔ = 2 sigma
        V[k] = Vk
        Lf[k] = _cholesky_psd(Vk)
        if cfg.skew_scale != 0.0:
            Z6 = _ds.zeta6(float(cos_gamma), float(lk))
            # two-delta collapse: node mu3 = Zeta6 / (2 sigma)^2
            Q[k] = cfg.skew_scale * _ds.solve_Q(Vk, Z6 / (2.0 * sig) ** 2)
    return bg, lam, dlam, resp, V, Lf, Q, sig, builder


def simulate_crn(cfg: MCConfig, gamma_arcmin: float) -> CRNResult:
    """Common-random-number run: Gaussian and skewed arms share the base AR(1)
    field (skewed = demo2-style deformed Gaussian), so xi_skew - xi_gauss cancels
    the O0 variance floor per realization and isolates the FK (3-cumulant) channel.
    """
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, V, Lf, Q, sig, builder = _field_precompute(cfg, cos_gamma)
    rho = float(np.exp(-dlam / sig))
    sqrt1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    rng = np.random.default_rng(cfg.seed)
    skew = cfg.skew_scale != 0.0
    N = cfg.n_lambda
    n_batches = max(1, cfg.n_real // cfg.batch_size)

    acc_g = np.zeros((6, 6)); acc_s = np.zeros((6, 6))
    acc_d = np.zeros((6, 6)); acc_d2 = np.zeros((6, 6))   # diff mean + sq for SE
    n_good = 0
    for _ in range(n_batches):
        m = cfg.batch_size
        z = np.zeros((6, m))                  # base Gaussian AR(1) field
        s_g = np.zeros((6, m)); s_s = np.zeros((6, m))
        kap_g = np.zeros((6, m)); kap_s = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            innov = (Lf[k] @ rng.standard_normal((6, m)))
            z = innov if k == 0 else (rho * z + sqrt1mr2 * innov)
            f_g = z
            if skew:
                ff = np.einsum("im,jm->ijm", z, z) - V[k][:, :, None]
                f_s = z + 0.5 * np.einsum("aij,ijm->am", Q[k], ff)
            else:
                f_s = z
            # evolve both arms with the same response; field enters as f*dlam
            s_g = _step(s_g, resp[k], f_g, dlam, cfg.use_f_vertex, m)
            s_s = _step(s_s, resp[k], f_s, dlam, cfg.use_f_vertex, m)
            w = 0.5 * dlam if (k == 0 or k == N - 1) else dlam
            kap_g += w * s_g; kap_s += w * s_s
            blown |= np.any(np.abs(s_s) > cfg.blowup, axis=0)
        good = ~blown
        gg = kap_g[:, good]; ss = kap_s[:, good]
        acc_g += gg @ gg.T; acc_s += ss @ ss.T
        # per-realization difference of the products (CRN)
        dprod = np.einsum("im,jm->ijm", ss, ss) - np.einsum("im,jm->ijm", gg, gg)
        acc_d += dprod.sum(axis=2)
        acc_d2 += (dprod * dprod).sum(axis=2)
        n_good += int(good.sum())

    ng = max(n_good, 1)
    xi_g6 = acc_g / ng; xi_s6 = acc_s / ng
    fk6 = acc_d / ng
    fk_var = acc_d2 / ng - fk6 * fk6
    fk_se6 = np.sqrt(np.maximum(fk_var, 0.0) / ng)
    return CRNResult(gamma_arcmin, cos_gamma, sig,
                     xi_g6[0:3, 3:6], xi_s6[0:3, 3:6],
                     fk6[0:3, 3:6], fk_se6[0:3, 3:6], n_good)


def _step(s, resp_k, f, dlam, use_f, m):
    """One exact-propagator step: response into node, then inject F-vertex + field."""
    if use_f:
        s_cr = s.reshape(2, 3, m).transpose(1, 0, 2)
        fv = _f_vertex(s_cr).transpose(1, 0, 2).reshape(6, m)
        return resp_k * s + fv * dlam + f * dlam
    return resp_k * s + f * dlam


# ---------------------------------------------------------------------------
# Path A -- perturbative-in-zeta FK estimator (no Riccati blow-up)
# ---------------------------------------------------------------------------
def _f_vertex_lin(sg: NDArray[np.float64],
                  sd: NDArray[np.float64]) -> NDArray[np.float64]:
    """Linearised F-vertex DF(sg).sd = F_abc (sg_b sd_c + sd_b sg_c), per ray.

    sg, sd shaped (3, ...).  This is the directional derivative of ``_f_vertex``
    about the Gaussian trajectory sg, acting on the increment-driven field sd.
    """
    g0, g1, g2 = sg[0], sg[1], sg[2]
    d0, d1, d2 = sd[0], sd[1], sd[2]
    out = np.empty_like(sd)
    out[0] = -2.0 * (g0 * d0 + g1 * d1 + g2 * d2)
    out[1] = -2.0 * (g0 * d1 + d0 * g1)
    out[2] = -2.0 * (g0 * d2 + d0 * g2)
    return out


@dataclass
class PathAResult:
    gamma_arcmin: float
    cos_gamma: float
    sigma_lambda: float
    xi_gauss: NDArray[np.float64]     # (3,3) <kappa_g(n1) kappa_g(n2)> (O0+FF)
    fk: NDArray[np.float64]           # (3,3) FK = <kg x kd + kd x kg> cross-ray
    fk_err: NDArray[np.float64]       # (3,3) standard error of FK
    n_good: int


def simulate_fk_pathA(cfg: MCConfig, gamma_arcmin: float) -> PathAResult:
    """Perturbative-in-zeta FK estimator (Path A).

    Evolves two coupled fields per realization, sharing the AR(1) Gaussian base
    ``z`` (so the estimator is intrinsically CRN / low-variance):

      * ``s_g`` : the full nonlinear Gaussian arm (F-vertex on), kappa_g = int s_g.
      * ``s_d`` : a LINEAR field driven by the skew increment
        ``df = 0.5 Q (z x z - V)`` through the F-vertex LINEARISED about s_g:
        ``ds_d = resp s_d + DF(s_g).s_d dlam + df dlam``.

    To leading order in the injected 3-cumulant the skewed convergence is
    ``kappa_s ~ kappa_g + kappa_d`` (kappa_d linear in df), so

        FK_ab = <kappa_g_a(n1) kappa_d_b(n2)> + <kappa_d_a(n1) kappa_g_b(n2)>

    which this routine accumulates as the symmetric cross-ray block of
    ``<kappa_g x kappa_d + kappa_d x kappa_g>``.  Because s_d obeys a LINEAR ODE
    with bounded coefficient DF(s_g)~s_g, the Levy-large node skewness cannot
    blow up the Riccati (the failure mode of ``simulate_crn``).

    Requires ``sigma_lambda > 0`` and ``skew_scale != 0``.
    """
    if cfg.skew_scale == 0.0:
        raise ValueError("Path A requires skew_scale != 0 (inject the 3-cumulant)")
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, V, Lf, Q, sig, builder = _field_precompute(cfg, cos_gamma)
    rho = float(np.exp(-dlam / sig))
    sqrt1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    rng = np.random.default_rng(cfg.seed)
    N = cfg.n_lambda
    n_batches = max(1, cfg.n_real // cfg.batch_size)

    acc_g = np.zeros((6, 6))           # <kappa_g x kappa_g>  (O0 + FF)
    acc_p = np.zeros((6, 6))           # <kg x kd + kd x kg>   (symmetric FK)
    acc_p2 = np.zeros((6, 6))          # squared, for the SE
    n_good = 0
    for _ in range(n_batches):
        m = cfg.batch_size
        z = np.zeros((6, m))
        s_g = np.zeros((6, m)); s_d = np.zeros((6, m))
        kap_g = np.zeros((6, m)); kap_d = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            innov = Lf[k] @ rng.standard_normal((6, m))
            z = innov if k == 0 else (rho * z + sqrt1mr2 * innov)
            # skew increment df = 0.5 Q (z x z - V)
            ff = np.einsum("im,jm->ijm", z, z) - V[k][:, :, None]
            df = 0.5 * np.einsum("aij,ijm->am", Q[k], ff)
            # Gaussian arm: full F-vertex, driven by z
            s_g = _step(s_g, resp[k], z, dlam, cfg.use_f_vertex, m)
            # increment arm: LINEARISED F-vertex about s_g, driven by df
            if cfg.use_f_vertex:
                sg_cr = s_g.reshape(2, 3, m).transpose(1, 0, 2)
                sd_cr = s_d.reshape(2, 3, m).transpose(1, 0, 2)
                dfv = _f_vertex_lin(sg_cr, sd_cr).transpose(1, 0, 2).reshape(6, m)
                s_d = resp[k] * s_d + dfv * dlam + df * dlam
            else:
                s_d = resp[k] * s_d + df * dlam
            w = 0.5 * dlam if (k == 0 or k == N - 1) else dlam
            kap_g += w * s_g; kap_d += w * s_d
            blown |= np.any(np.abs(s_g) > cfg.blowup, axis=0)
        good = ~blown
        gg = kap_g[:, good]; dd = kap_d[:, good]
        acc_g += gg @ gg.T
        prod = np.einsum("im,jm->ijm", gg, dd)
        prod = prod + prod.transpose(1, 0, 2)         # kg x kd + kd x kg
        acc_p += prod.sum(axis=2)
        acc_p2 += (prod * prod).sum(axis=2)
        n_good += int(good.sum())

    ng = max(n_good, 1)
    xi_g6 = acc_g / ng
    fk6 = acc_p / ng
    fk_var = acc_p2 / ng - fk6 * fk6
    fk_se6 = np.sqrt(np.maximum(fk_var, 0.0) / ng)
    return PathAResult(gamma_arcmin, cos_gamma, sig,
                       xi_g6[0:3, 3:6], fk6[0:3, 3:6], fk_se6[0:3, 3:6], n_good)


@dataclass
class FKVRResult:
    gamma_arcmin: float
    cos_gamma: float
    sigma_lambda: float
    fk: NDArray[np.float64]           # (3,3) cross-ray FK
    fk_err: NDArray[np.float64]       # (3,3) batch-jackknife SE
    n_good: int


def simulate_fk_vr(cfg: MCConfig, gamma_arcmin: float, n_win: int = 64) -> FKVRResult:
    """Variance-reduced FK estimator (Q pulled OUTSIDE the MC average).

    The FK diagram connects obs_b to one zeta leg F-FREE, so the no-F windowed
    response ``kappa_d^(0)_B = int W delta_f_B`` suffices (no s_d ODE). Writing
    ``delta_f_B(k) = 0.5 sum_ij Q[k]_Bij w_ij(k)``, ``w_ij = z_i z_j - V_ij``,

        <kappa_g_A kappa_d^(0)_B> = 0.5 dlam sum_k W(k) sum_ij Q[k]_Bij C_{A,ij}(k),
        C_{A,ij}(k) = <kappa_g_A . w_ij(k)>   (all O(1) -> controlled variance).

    Q (huge for the Levy zeta) multiplies the AVERAGED C, never the per-sample w.
    FK = T + T^T (the x<->y symmetrisation; T_AB = the contraction above).
    """
    if cfg.skew_scale == 0.0:
        raise ValueError("FK VR estimator requires skew_scale != 0")
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, V, Lf, Q, sig, builder = _field_precompute(cfg, cos_gamma)
    rho = float(np.exp(-dlam / sig))
    sqrt1mr2 = float(np.sqrt(max(1.0 - rho * rho, 0.0)))
    W = _ds.order0_window(bg, lam, cfg.lam_source, n_gauss=n_win)   # (N,)
    V_knn = np.transpose(V, (1, 2, 0))                              # (6,6,N)
    rng = np.random.default_rng(cfg.seed)
    N = cfg.n_lambda
    n_batches = max(1, cfg.n_real // cfg.batch_size)

    T_batches = []
    n_good = 0
    for _ in range(n_batches):
        m = cfg.batch_size
        z = np.zeros((6, m))
        s_g = np.zeros((6, m)); s_g0 = np.zeros((6, m))
        kap_g = np.zeros((6, m)); kap_g0 = np.zeros((6, m))
        zstore = np.empty((N, 6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            innov = Lf[k] @ rng.standard_normal((6, m))
            z = innov if k == 0 else (rho * z + sqrt1mr2 * innov)
            zstore[k] = z
            s_g = _step(s_g, resp[k], z, dlam, cfg.use_f_vertex, m)
            s_g0 = _step(s_g0, resp[k], z, dlam, False, m)         # F-OFF arm (shared z)
            w_ = 0.5 * dlam if (k == 0 or k == N - 1) else dlam
            kap_g += w_ * s_g; kap_g0 += w_ * s_g0
            blown |= np.any(np.abs(s_g) > cfg.blowup, axis=0)
        good = ~blown
        # CRN variance reduction: the FK couples only to the F-INDUCED part of
        # kappa_g (Cov(L1, z z) = 0 in mean but adds noise); kap_g - kap_g0 keeps
        # the L2+ signal and drops the dominant L1 (Order-0) variance (~1000x).
        kg = (kap_g - kap_g0)[:, good]                             # (6, mb)
        mb = int(good.sum())
        meankap = kg.mean(axis=1)                                  # (6,)
        # T_AB = 0.5 dlam sum_k W(k) sum_ij Q[k]_Bij <kappa_g_A w_ij(k)>_connected,
        # with <kg_A w_ij(k)>_c = Cov(kappa_g_A, z_i z_j) = the SAMPLE covariance.
        # The disconnected <kg_A><z_i z_j> is ~1e7x the connected signal, so it
        # MUST be removed with the matched SAMPLE <z_i z_j> (not the theoretical
        # V[k]; the AR(1) sample covariance differs from V[k] by transient/
        # discretization terms that would otherwise swamp the signal).
        T = np.zeros((6, 6))
        for k in range(N):
            zk = zstore[k][:, good]                                # (6, mb)
            C_k = np.einsum("Am,im,jm->Aij", kg, zk, zk) / mb      # <kg z z>_batch
            zz_k = (zk @ zk.T) / mb                                # <z z>_batch (6,6)
            Cc_k = C_k - meankap[:, None, None] * zz_k[None, :, :]  # connected cov
            T += (0.5 * dlam * W[k]) * np.einsum("Bij,Aij->AB", Q[k], Cc_k)
        T_batches.append(T)
        n_good += mb

    T_arr = np.array(T_batches)                                    # (nb,6,6)
    T_mean = T_arr.mean(axis=0)
    FK6 = T_mean + T_mean.T
    # batch SE of the symmetrised cross-ray block
    FK6_batches = T_arr + np.transpose(T_arr, (0, 2, 1))
    se6 = FK6_batches.std(axis=0) / np.sqrt(max(len(T_batches), 1))
    return FKVRResult(gamma_arcmin, cos_gamma, sig,
                      FK6[0:3, 3:6], se6[0:3, 3:6], n_good)


@dataclass
class FFResult:
    gamma_arcmin: float
    cos_gamma: float
    ff_moment: NDArray[np.float64]    # (3,3) <kk>_Fon - <kk>_Foff (cross-ray MOMENT)
    ff_conn: NDArray[np.float64]      # (3,3) moment - disconnected (CONNECTED FF)
    ff_disc: NDArray[np.float64]      # (3,3) <k>_Fon x <k>_Fon - <k>_Foff x <k>_Foff
    ff_moment_err: NDArray[np.float64]
    mean_kappa_on: NDArray[np.float64]  # (3,) <kappa>_Fon ray-1 (mean 2nd-order conv)
    n_good: int


def simulate_ff_crn(cfg: MCConfig, gamma_arcmin: float) -> FFResult:
    """Isolate the MC FF channel via an F-on/F-off CRN difference (Gaussian, white).

    Both arms share the SAME white-noise draws; arm A applies the F-vertex, arm B
    does not.  The per-realization product difference gives the FF MOMENT
    ``<kk>_Fon - <kk>_Foff`` (all orders of F) at low variance (the O0 floor
    cancels).  We also track the arm MEANS so we can split off the DISCONNECTED
    piece ``<kappa>_Fon x <kappa>_Fon`` (the square of the mean 2nd-order
    convergence), leaving the CONNECTED FF correlation.  ``skew_scale`` and
    ``sigma_lambda`` are forced to the Gaussian / white-noise values.
    """
    from dataclasses import replace
    cfg = replace(cfg, skew_scale=0.0, sigma_lambda=0.0)
    cos_gamma = float(np.cos(np.deg2rad(gamma_arcmin / 60.0)))
    bg, lam, dlam, resp, L_white, _Q, _sig, builder = _precompute(cfg, cos_gamma)
    rng = np.random.default_rng(cfg.seed)
    N = cfg.n_lambda
    n_batches = max(1, cfg.n_real // cfg.batch_size)

    acc_dmom = np.zeros((6, 6)); acc_dmom2 = np.zeros((6, 6))
    sum_on = np.zeros(6); sum_off = np.zeros(6)
    n_good = 0
    for _ in range(n_batches):
        m = cfg.batch_size
        s_on = np.zeros((6, m)); s_off = np.zeros((6, m))
        kap_on = np.zeros((6, m)); kap_off = np.zeros((6, m))
        blown = np.zeros(m, dtype=bool)
        for k in range(N):
            noise = _draw_increment(rng, L_white[k], None, m, False)  # shared
            s_cr = s_on.reshape(2, 3, m).transpose(1, 0, 2)
            fv = _f_vertex(s_cr).transpose(1, 0, 2).reshape(6, m)
            s_on = resp[k] * s_on + fv * dlam + noise
            s_off = resp[k] * s_off + noise
            w = 0.5 * dlam if (k == 0 or k == N - 1) else dlam
            kap_on += w * s_on; kap_off += w * s_off
            blown |= np.any(np.abs(s_on) > cfg.blowup, axis=0)
        good = ~blown
        on = kap_on[:, good]; off = kap_off[:, good]
        dprod = np.einsum("im,jm->ijm", on, on) - np.einsum("im,jm->ijm", off, off)
        acc_dmom += dprod.sum(axis=2); acc_dmom2 += (dprod * dprod).sum(axis=2)
        sum_on += on.sum(axis=1); sum_off += off.sum(axis=1)
        n_good += int(good.sum())

    ng = max(n_good, 1)
    ff_mom6 = acc_dmom / ng
    mean_on6 = sum_on / ng; mean_off6 = sum_off / ng
    disc6 = np.outer(mean_on6, mean_on6) - np.outer(mean_off6, mean_off6)
    conn6 = ff_mom6 - disc6
    mom_var = acc_dmom2 / ng - ff_mom6 * ff_mom6
    mom_se6 = np.sqrt(np.maximum(mom_var, 0.0) / ng)
    return FFResult(gamma_arcmin, cos_gamma,
                    ff_mom6[0:3, 3:6], conn6[0:3, 3:6], disc6[0:3, 3:6],
                    mom_se6[0:3, 3:6], mean_on6[0:3], n_good)


if __name__ == "__main__":
    # Validation step 3 (Gaussian): MC kk (F off) must reproduce order0_mc.
    print("Gaussian MC vs analytic Order-0 (F off, white noise):")
    print(f"{'gamma[arcmin]':>14} {'MC kk':>14} {'+/-SE':>12} "
          f"{'order0_mc':>14} {'ratio':>8}")
    cfg = MCConfig(n_real=40_000, batch_size=2_000, n_lambda=500,
                   use_f_vertex=False, skew_scale=0.0, sigma_lambda=0.0)
    for g in (0.5, 1.0, 10.0, 100.0, 1000.0):
        res = simulate(cfg, g)
        kk, se = res.xi[0, 0], res.xi_err[0, 0]
        ratio = kk / res.o0_analytic if res.o0_analytic else float("nan")
        print(f"{g:>14.3f} {kk:>14.6e} {se:>12.2e} "
              f"{res.o0_analytic:>14.6e} {ratio:>8.4f}  (n_good={res.n_good})")
