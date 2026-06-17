#!/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python
"""DECISIVE: deployed canoes HIGH zeta_D fold vs a COMPLEX-kernel variant vs ref.

Interpreter: /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python

Faithfully RE-IMPLEMENTS the deployed compute_kappa3_mod_sigma3_high /
compute_kappa3_sigma3_high HIGH inner loop (kappa3.py:3945-4053 / 3945-3623),
byte-for-byte in the integrand assembly, with ONE switch: the orientation kernel
returns either np.real(I) (deployed) or the full COMPLEX I.  Everything else
(B_delta, response, growth^4, measure_weights, prefactor, radial_factor, h^6) is
identical.  HIGH only (HIGH dominates zeta_D by ~8 orders).

The deployed function casts `float(np.sum(weighted*alpha))` (kappa3.py:4039),
which silently drops Im; that float() IS the operative truncation.  Here we keep
complex through the (u,v,phi) sum and the radial fold, taking |.| only at the
very end (as the natural components are genuinely complex).

We fold like stage1_ours_v2: |Gamma^mu| = |fold(K^3 * zeta_mu)| (derived_phase
is unit-modulus, so it cannot change |Gamma|).  We validate the deployed_real
branch reproduces the deployed compute_kappa3_mod_sigma3_high to ~1e-12 first.

Loads only our own freshly-written reference npz (trusted) with allow_pickle.
"""
from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
OUT = ANALYSIS / "outputs"

sys.path.insert(0, str(Path("/Users/zzhang/projects/canoes")))
_ETL = (Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
             "callables/kappa3_vertex/equal_time_limber"))
sys.path.insert(0, str(_ETL))
from _local_cosmo_pk import FiducialCosmology, load_pk_delta  # noqa: E402
from canoes.cosmo.background import chi as chi_of_z  # noqa: E402
import canoes.sachs.kappa3 as k3  # noqa: E402
from canoes.sachs.kappa3 import (  # noqa: E402
    _validate_limber_high_quadrature, chi_of_lambda, build_lightcone,
    _kappa3_flat_triangle_vertices, _f2_kernel_numpy,
    _kappa3_limber_response_product_h, _lambda_jacobian,
    _kappa3_radial_density_conversion, _H0_H_PER_MPC,
)

_N_LAMBDA = 40
ARCMIN2RAD = math.pi / 10800.0
_DMOD_SPINS = {1: (-2, 2, 2), 2: (2, -2, 2), 3: (2, 2, -2)}
_PPP = (2, 2, 2)


def alpha_kernel(spins, u, v, phi, w, r2, r3, keep_complex):
    """Orientation integral; real (deployed) or complex, same algebra."""
    s1, s2, s3 = spins
    S = int(s1 + s2 + s3)
    Q = u * r2 + v * np.exp(-1j * phi) * r3
    qabs = np.abs(Q)
    if S == 0:
        out = 2.0 * math.pi * jv(0, qabs)            # deployed scalar branch
        return out.astype(np.float64) if not keep_complex else out.astype(np.complex128)
    beta = np.where(qabs > 0.0, np.angle(Q), 0.0)
    l1_rel = -(u + v * np.exp(1j * phi))
    psi = np.where(w > 0.0, np.angle(l1_rel), 0.0)
    I = (2.0 * math.pi * (1j ** S) * np.exp(1j * S * beta)
         * jv(S, qabs) * np.exp(1j * (s1 * psi + s3 * phi)))
    return np.real(I).astype(np.float64) if not keep_complex else I.astype(np.complex128)


def high_zeta(triples, lam_h, spins, response_channel, cosmo, pk, *, keep_complex):
    """Faithful re-impl of the deployed HIGH zeta for one spin triple/channel."""
    ell_nodes, ell_weights, phi_nodes, phi_weights = _validate_limber_high_quadrature(
        ell_cut=60, ell_high_max=1000, ell_min=1e-3, n_ell=96, n_phi=64)
    cos_arr = np.clip(np.atleast_2d(np.asarray(triples, float)), -1.0, 1.0)
    lam = np.atleast_1d(np.asarray(lam_h, float))
    chi_arr = np.maximum(chi_of_lambda(cosmo, lam, convention="project"), 1e-6)
    lightcone, one_plus_z = build_lightcone(cosmo, chi_arr)
    poisson_amp = (-1.5 * float(cosmo.Omega_m) * (_H0_H_PER_MPC ** 2) * one_plus_z)
    n_cos, n_lam = cos_arr.shape[0], chi_arr.size

    u = ell_nodes[:, None, None]; v = ell_nodes[None, :, None]; phi = phi_nodes[None, None, :]
    cph = np.cos(phi)
    w = np.sqrt(np.maximum(u * u + v * v + 2 * u * v * cph, 0.0))
    max_leg = np.maximum(np.maximum(u, v), w)
    win = ((max_leg > 60.0) & (max_leg <= 1000.0)).astype(np.float64)
    meas = (ell_weights[:, None, None] * ell_weights[None, :, None]
            * phi_weights[None, None, :] * u * v * win)
    sw = np.maximum(w, 1e-3); su = np.maximum(u, 1e-3); sv = np.maximum(v, 1e-3)
    c12 = np.where(sw > 0, -(u + v * cph) / sw, 0.0)
    c13 = np.where(sw > 0, -(v + u * cph) / sw, 0.0)
    c23 = cph
    pref = 1.0 / ((2 * math.pi) ** 4)
    t12 = np.arccos(cos_arr[:, 0]); t23 = np.arccos(cos_arr[:, 1]); t31 = np.arccos(cos_arr[:, 2])

    alphas = []
    for ic in range(n_cos):
        r2, r3 = _kappa3_flat_triangle_vertices(float(t12[ic]), float(t23[ic]), float(t31[ic]))
        alphas.append(alpha_kernel(spins, u, v, phi, w, r2, r3, keep_complex))

    dt = np.complex128 if keep_complex else np.float64
    out = np.zeros((n_cos, n_lam), dtype=dt)
    for il, chi_val in enumerate(chi_arr):
        cf = float(chi_val)
        k1, k2, k3v = sw / cf, su / cf, sv / cf
        p1 = np.asarray(pk(k1), float); p2 = np.asarray(pk(k2), float); p3 = np.asarray(pk(k3v), float)
        f12 = _f2_kernel_numpy(k1, k2, c12); f13 = _f2_kernel_numpy(k1, k3v, c13); f23 = _f2_kernel_numpy(k2, k3v, c23)
        bdelta = 2.0 * (f12 * p1 * p2 + f13 * p1 * p3 + f23 * p2 * p3)
        growth = float(lightcone.M[0, il] / one_plus_z[il])
        radial_base = meas * (growth ** 4) * bdelta / (cf ** 4) * pref
        resp = _kappa3_limber_response_product_h(response_channel, float(poisson_amp[il]), float(one_plus_z[il]))
        weighted = radial_base * resp
        for ic in range(n_cos):
            s = np.sum(weighted * alphas[ic])
            out[ic, il] = complex(s) if keep_complex else float(np.real(s))

    d_lambda_d_chi = _lambda_jacobian(one_plus_z, convention="project")
    radial_factor, _, _ = _kappa3_radial_density_conversion(d_lambda_d_chi, "lambda")
    out = out * radial_factor[None, :]
    out = out * (float(cosmo.h) ** 6)   # units=physical
    return out


def build_lambda_grid(cosmo, z_source, n_lambda):
    from scipy.integrate import cumulative_trapezoid
    Om, h = cosmo.Omega_m, cosmo.h
    zf = np.linspace(0.0, z_source, 6000)
    chif = np.array([chi_of_z(zi, Omega_m=Om, h=h) for zi in zf])
    af = 1.0 / (1.0 + zf)
    lam_f = cumulative_trapezoid(af ** 2, chif, initial=0.0)
    z = np.linspace(z_source / n_lambda, z_source, n_lambda)
    chi_phys = np.interp(z, zf, chif)
    lam_phys = np.interp(z, zf, lam_f)
    return z, chi_phys, lam_phys, lam_phys * h, float(chif[-1])


def conv_kernel(chi_phys, z, chi_s):
    a = 1.0 / (1.0 + z)
    return a ** 2 * chi_phys * (chi_s - chi_phys) / chi_s


def sas_triples(gamma_rad, phi_rad):
    t3 = 2.0 * gamma_rad * np.sin(phi_rad / 2.0)
    cg = np.cos(gamma_rad); ct3 = np.cos(t3)
    return np.stack([ct3, np.full_like(ct3, cg), np.full_like(ct3, cg)], axis=1)


def main():
    cosmo = FiducialCosmology()
    pk = load_pk_delta(n_s=cosmo.n_s)
    z, chi_phys, lam_phys, lam_h, chi_s = build_lambda_grid(cosmo, 5.0, _N_LAMBDA)
    K = conv_kernel(chi_phys, z, chi_s)

    def fold(zeta):
        return np.trapezoid((K ** 3)[None, :] * zeta, lam_phys, axis=1)

    ref = np.load(OUT / "stage1_bkappa_phase_ref.npz", allow_pickle=True)
    configs = list(zip(ref["gamma_arcmin"].tolist(), ref["phi_deg"].tolist()))

    # ---- self-validation: deployed_real re-impl == deployed function ----------
    print("VALIDATION: re-impl deployed_real == compute_kappa3_mod_sigma3_high?")
    gr = math.radians(10.0 / 60.0); tri = sas_triples(gr, np.array([math.radians(60.0)]))
    mine = high_zeta(tri, lam_h, (2, 2, -2), "PPP", cosmo, pk, keep_complex=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dep = k3.compute_kappa3_mod_sigma3_high(
            tri, lam_h, pk_matter_today=pk, cosmo=cosmo, ell_cut=60, ell_high_max=1000,
            ell_min=1e-3, n_ell=96, n_phi=64, units="physical",
            lambda_convention="project", radial_measure="lambda",
            _spin_weights_dmod=(2, 2, -2))[1]
    rel = np.max(np.abs(mine - dep) / (np.abs(dep) + 1e-300))
    print(f"  max rel-err (2,2,-2): {rel:.2e}  {'OK' if rel < 1e-9 else 'MISMATCH'}")

    print("\n" + "=" * 104)
    print("DEPLOYED (np.real) vs COMPLEX-kernel HIGH zeta_D fold vs B_kappa-phase ref")
    print("=" * 104)
    rows = []
    for variant, kc in (("deployed_real", False), ("complex", True)):
        print(f"\n############ VARIANT = {variant} ############")
        for (gam, phid) in configs:
            grr = math.radians(gam / 60.0)
            tri = sas_triples(grr, np.array([math.radians(phid)]))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                zP = high_zeta(tri, lam_h, _PPP, "PPP", cosmo, pk, keep_complex=kc)
                G = {0: fold(zP)}
                for mu, spw in _DMOD_SPINS.items():
                    zd = high_zeta(tri, lam_h, spw, "PPP", cosmo, pk, keep_complex=kc)
                    G[mu] = fold(zd)
            fi = math.sqrt(sum(abs(G[m][0]) ** 2 for m in range(4)))
            idx = configs.index((gam, phid))
            fr = math.sqrt(sum(abs(ref[f"Gamma{m}"][idx]) ** 2 for m in range(4)))
            tag = " <-- PRIMARY" if (gam, phid) == (10.0, 60.0) else ""
            print(f"  g={gam:.0f}' phi={phid:.0f}{tag}: "
                  + " ".join(f"|G{m}|={abs(G[m][0]):.3e}" for m in range(4))
                  + f" | fi={fi:.3e} ref={fr:.3e} can/ref={fi/fr:.3f}")
            rows.append((variant, gam, phid, fi, fr, fi / fr,
                         *[complex(G[m][0]) for m in range(4)]))

    print("\n" + "=" * 104)
    print("VERDICT (frame-invariant can/ref):")
    for variant in ("deployed_real", "complex"):
        rr = [r[5] for r in rows if r[0] == variant]
        print(f"  {variant:14s}: median={np.median(rr):.3f}  range=[{min(rr):.3f}, {max(rr):.3f}]")

    np.savez(
        OUT / "HIGH_deployed_vs_complex.npz",
        variant=np.array([r[0] for r in rows]), gam=np.array([r[1] for r in rows]),
        phi=np.array([r[2] for r in rows]), fi=np.array([r[3] for r in rows]),
        fr=np.array([r[4] for r in rows]), can_over_ref=np.array([r[5] for r in rows]),
        G0=np.array([r[6] for r in rows], dtype=np.complex128),
        G1=np.array([r[7] for r in rows], dtype=np.complex128),
        G2=np.array([r[8] for r in rows], dtype=np.complex128),
        G3=np.array([r[9] for r in rows], dtype=np.complex128),
        note="deployed np.real(I) vs complex-I HIGH zeta_D folded like stage1_ours_v2.")
    print(f"\nsaved -> {OUT / 'HIGH_deployed_vs_complex.npz'}")


if __name__ == "__main__":
    main()
