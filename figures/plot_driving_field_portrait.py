"""Driving-field map chain schematic for the STF_lensing paper.

The figure starts from a synthetic density contrast map drawn from a CAMB
linear matter power spectrum, solves the scalar Poisson potential, then
derives the Ricci-focusing trace and Weyl-shear traceless Hessian.

Output:
    figures/driving_field_portrait.{pdf,png}
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Circle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIGURES = ROOT / "figures"

sys.path.insert(0, str(HERE))
from _plot_style import apply_rcparams  # noqa: E402

apply_rcparams()
mpl.rcParams.update({"axes.grid": False})

BLUE = "#2f6f9f"
CORAL = "#b95f58"


@dataclass(frozen=True)
class PkMetadata:
    source: str
    redshift: float
    h: float
    omega_b_h2: float
    omega_c_h2: float
    n_s: float
    sigma8: float


def _normalise(field: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=float)
    field = field - np.mean(field)
    sigma = np.std(field)
    if sigma == 0:
        return field
    return field / sigma


def _robust_symmetric_norm(field: np.ndarray, percentile: float = 99.0) -> TwoSlopeNorm:
    vmax = float(np.percentile(np.abs(field), percentile))
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def _k_grid_3d(
    n_pix: int,
    box_mpc: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dx = box_mpc / n_pix
    kx = 2.0 * np.pi * np.fft.fftfreq(n_pix, d=dx)[:, None, None]
    ky = 2.0 * np.pi * np.fft.fftfreq(n_pix, d=dx)[None, :, None]
    kz = 2.0 * np.pi * np.fft.rfftfreq(n_pix, d=dx)[None, None, :]
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    return kx, ky, kz, k


def _camb_linear_pk(k: np.ndarray, redshift: float) -> tuple[np.ndarray, PkMetadata]:
    """Evaluate a Planck-like linear matter P(k) in Mpc units."""
    import camb
    from camb import model

    h = 0.6711
    omega_b_h2 = 0.02207
    omega_c_h2 = 0.12029
    n_s = 0.97
    sigma8_target = 0.81

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=100.0 * h,
        ombh2=omega_b_h2,
        omch2=omega_c_h2,
        mnu=0.06,
        omk=0.0,
        tau=0.06,
    )
    pars.InitPower.set_params(As=2.1e-9, ns=n_s)
    pars.set_matter_power(redshifts=[redshift], kmax=max(6.0, float(np.nanmax(k)) * 1.05))
    pars.NonLinear = model.NonLinear_none

    interpolator = camb.get_matter_power_interpolator(
        pars,
        nonlinear=False,
        hubble_units=False,
        k_hunit=False,
        kmax=max(6.0, float(np.nanmax(k)) * 1.05),
        zmax=max(2.1, redshift + 0.1),
    )
    k_eval = np.clip(k, 1e-4, None)
    pk = np.asarray(interpolator.P(redshift, k_eval), dtype=float)
    pk[k == 0] = 0.0
    meta = PkMetadata(
        source="CAMB linear Pm(k)",
        redshift=redshift,
        h=h,
        omega_b_h2=omega_b_h2,
        omega_c_h2=omega_c_h2,
        n_s=n_s,
        sigma8=sigma8_target,
    )
    return pk, meta


def _bbks_fallback_pk(k: np.ndarray, redshift: float) -> tuple[np.ndarray, PkMetadata]:
    h = 0.6711
    omega_b_h2 = 0.02207
    omega_c_h2 = 0.12029
    omega_m = (omega_b_h2 + omega_c_h2) / h**2
    n_s = 0.97
    gamma = omega_m * h
    q = np.clip(k / gamma, 1e-5, None)
    transfer = (
        np.log1p(2.34 * q)
        / (2.34 * q)
        / (1.0 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** 0.25
    )
    growth = 1.0 / (1.0 + redshift)
    pk = growth * growth * np.clip(k, 1e-5, None) ** n_s * transfer * transfer
    pk[k == 0] = 0.0
    meta = PkMetadata(
        source="BBKS fallback Pm(k)",
        redshift=redshift,
        h=h,
        omega_b_h2=omega_b_h2,
        omega_c_h2=omega_c_h2,
        n_s=n_s,
        sigma8=0.81,
    )
    return pk, meta


def linear_matter_pk(k: np.ndarray, redshift: float) -> tuple[np.ndarray, PkMetadata]:
    try:
        return _camb_linear_pk(k, redshift)
    except Exception as exc:
        print(f"[driving-field] CAMB unavailable, using BBKS fallback: {exc}")
        return _bbks_fallback_pk(k, redshift)


def make_fields(
    seed: int,
    *,
    n_pix: int,
    box_mpc: float,
    redshift: float,
) -> tuple[dict[str, np.ndarray], PkMetadata]:
    rng = np.random.default_rng(seed)
    dx = box_mpc / n_pix
    kx, ky, _kz, k = _k_grid_3d(n_pix, box_mpc)
    pk, meta = linear_matter_pk(k, redshift)

    omega_m = (meta.omega_b_h2 + meta.omega_c_h2) / (meta.h * meta.h)
    # Geometric units (c=1): the Hubble rate enters as an inverse comoving
    # length.  Numerically that is H0 [km/s/Mpc] / c [km/s] = H0 in Mpc^-1.
    h0_inv_mpc = 100.0 * meta.h / 299792.458
    a = 1.0 / (1.0 + redshift)
    poisson_a = 1.5 * omega_m * h0_inv_mpc * h0_inv_mpc / a  # A = (3/2) Om H0^2 / a  [Mpc^-2]

    white = rng.normal(size=(n_pix, n_pix, n_pix))
    white_k = np.fft.rfftn(white, axes=(0, 1, 2))
    k_ny = np.pi / dx
    pixel_window = np.exp(-0.5 * (k / (0.90 * k_ny)) ** 8)
    delta_k = white_k * np.sqrt(np.maximum(pk, 0.0) / dx**3) * pixel_window
    delta_k[0, 0, 0] = 0.0
    delta = np.fft.irfftn(delta_k, s=(n_pix, n_pix, n_pix), axes=(0, 1, 2))

    phi_k = np.zeros_like(delta_k, dtype=np.complex128)
    mask = k > 0
    phi_k[mask] = -poisson_a * delta_k[mask] / (k[mask] * k[mask])
    phi = np.fft.irfftn(phi_k, s=(n_pix, n_pix, n_pix), axes=(0, 1, 2))

    ricci = -poisson_a * delta
    gamma1_k = -0.5 * (kx * kx - ky * ky) * phi_k
    gamma2_k = -(kx * ky) * phi_k
    gamma1 = np.fft.irfftn(gamma1_k, s=(n_pix, n_pix, n_pix), axes=(0, 1, 2))
    gamma2 = np.fft.irfftn(gamma2_k, s=(n_pix, n_pix, n_pix), axes=(0, 1, 2))

    slab = n_pix // 2
    psi = gamma1[:, :, slab] + 1j * gamma2[:, :, slab]
    return {
        "delta": delta[:, :, slab],
        "phi": phi[:, :, slab],
        "ricci": ricci[:, :, slab],
        "psi": psi,
        "poisson_a": np.array(poisson_a),
    }, meta


def _field_position(ix: int, iy: int, n_pix: int, extent: tuple[float, float, float, float]) -> tuple[float, float]:
    x0, x1, y0, y1 = extent
    x = x0 + (ix + 0.5) / n_pix * (x1 - x0)
    y = y0 + (iy + 0.5) / n_pix * (y1 - y0)
    return x, y


def _glyph_idx(n_pix: int, n_ticks: int) -> np.ndarray:
    """Regular sub-sampling indices for the panel-D Psi_0 glyph grid.

    Shared by ``render`` (to fix the single length scale from the sampled
    maximum) and by ``_component_segments`` (to draw the glyphs), so that both
    operate on the identical set of sample points.
    """
    return np.linspace(10, n_pix - 11, n_ticks).astype(int)


def _peak_marker(delta: np.ndarray, extent: tuple[float, float, float, float]) -> tuple[float, float]:
    margin = max(16, delta.shape[0] // 8)
    inner = delta[margin:-margin, margin:-margin]
    iy, ix = np.unravel_index(np.argmax(inner), inner.shape)
    return _field_position(ix + margin, iy + margin, delta.shape[0], extent)


def _component_segments(
    psi: np.ndarray,
    *,
    extent: tuple[float, float, float, float],
    component: str,
    scale: float,
    max_len: float,
    n_ticks: int = 11,
    cutoff: float = 0.05,
) -> LineCollection:
    x0, x1, y0, y1 = extent
    n_pix = psi.shape[0]
    idx = _glyph_idx(n_pix, n_ticks)
    xs = x0 + (idx + 0.5) / n_pix * (x1 - x0)
    ys = y0 + (idx + 0.5) / n_pix * (y1 - y0)

    segments: list[list[tuple[float, float]]] = []
    widths: list[float] = []
    if component == "real":
        values = np.real(psi)
        positive_angle, negative_angle = 0.0, 0.5 * np.pi
        color = BLUE
    elif component == "imag":
        values = np.imag(psi)
        positive_angle, negative_angle = 0.25 * np.pi, -0.25 * np.pi
        color = CORAL
    else:
        raise ValueError(f"unknown component {component!r}")

    for iy, y in zip(idx, ys):
        for ix, x in zip(idx, xs):
            val = values[iy, ix]
            # Length is strictly proportional to the local component amplitude.
            # ``scale`` is the sampled maximum over both components (set in
            # ``render``), so ``strength`` already lies in [0, 1]: the strongest
            # sampled cell draws exactly ``max_len`` and there is no saturation
            # clamp or length floor to compress the dynamic range.
            strength = abs(val) / scale
            if strength < cutoff:
                continue
            theta = positive_angle if val >= 0.0 else negative_angle
            length = max_len * strength
            dx = 0.5 * length * np.cos(theta)
            dy = 0.5 * length * np.sin(theta)
            segments.append([(x - dx, y - dy), (x + dx, y + dy)])
            widths.append(0.52 + 1.18 * strength)

    return LineCollection(
        segments,
        colors=color,
        linewidths=widths,
        alpha=0.88,
        capstyle="round",
        zorder=4,
    )


def _panel_label(ax, text: str, formula: str) -> None:
    ax.text(
        0.035,
        0.965,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14.5,
        color="0.12",
        bbox=dict(facecolor="white", edgecolor="0.82", pad=4.0),
        zorder=8,
    )
    ax.text(
        0.035,
        0.050,
        formula,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12.5,
        color="0.20",
        bbox=dict(facecolor="white", edgecolor="0.86", pad=3.2, alpha=0.92),
        zorder=8,
    )


def render(seed: int, out_stem: Path, *, n_pix: int, box_mpc: float, redshift: float) -> None:
    fields, meta = make_fields(seed, n_pix=n_pix, box_mpc=box_mpc, redshift=redshift)
    half_box = 0.5 * box_mpc
    extent = (-half_box, half_box, -half_box, half_box)
    peak_x, peak_y = _peak_marker(fields["delta"], extent)

    fig = plt.figure(figsize=(14.4, 5.6))
    gs = fig.add_gridspec(
        2,
        4,
        height_ratios=(1.0, 0.06),
        left=0.060,
        right=0.985,
        bottom=0.150,
        top=0.890,
        wspace=0.110,
        hspace=0.205,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    caxes = [fig.add_subplot(gs[1, i]) for i in range(4)]

    panels = [
        ("A  overdensity δₘ", "CAMB linear Pₘ(k)", fields["delta"], "RdBu_r",
         _robust_symmetric_norm(fields["delta"]), "δₘ"),
        ("B  potential Φ", "Φ(k) = -A δₘ(k)/k²", fields["phi"] / 1e-5, "RdBu_r",
         _robust_symmetric_norm(fields["phi"] / 1e-5), "Φ [10⁻⁵]"),
        ("C  Ricci Φ₀₀", "Φ₀₀ ≃ -A δₘ", fields["ricci"] / 1e-8, "RdBu_r",
         _robust_symmetric_norm(fields["ricci"] / 1e-8), "Φ₀₀ [10⁻⁸ Mpc⁻²]"),
    ]

    for i, (label, formula, field, cmap, norm, cbar_label) in enumerate(panels):
        im = axes[i].imshow(
            field,
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=norm,
            interpolation="bilinear",
            rasterized=True,
        )
        _panel_label(axes[i], label, formula)
        axes[i].add_patch(Circle((peak_x, peak_y), 17.5, fill=False,
                                 edgecolor="0.12", linewidth=0.85, zorder=7))
        cb = fig.colorbar(im, cax=caxes[i], orientation="horizontal")
        cb.set_label(cbar_label, fontsize=13.5, labelpad=2.5)
        cb.ax.tick_params(labelsize=12.0, length=3.2, pad=2.0)

    psi = fields["psi"]
    psi_mag = np.abs(psi)
    psi_scale = 1e-8
    im_psi = axes[3].imshow(
        np.clip(psi_mag / psi_scale, 0.0, 8.0),
        origin="lower",
        extent=extent,
        cmap="Greys",
        norm=Normalize(vmin=0.0, vmax=8.0),
        alpha=0.42,
        interpolation="bilinear",
        rasterized=True,
    )
    # Coarse, regular glyph grid for Re/Im Psi_0.  Fix a single shared length
    # scale from the sampled maximum of both components so that stick length is
    # strictly proportional to the local amplitude (the strongest sampled cell
    # maps to max_len; everything else is linearly shorter).
    glyph_n_ticks = 11
    glyph_idx = _glyph_idx(psi.shape[0], glyph_n_ticks)
    re_sampled = np.abs(np.real(psi)[np.ix_(glyph_idx, glyph_idx)])
    im_sampled = np.abs(np.imag(psi)[np.ix_(glyph_idx, glyph_idx)])
    component_scale = max(float(re_sampled.max()), float(im_sampled.max()), 1e-12)
    # Longest stick spans 0.85 of the grid pitch, so even two adjacent maxima
    # stop short of each other and no glyphs overlap.
    glyph_pitch = (glyph_idx[1] - glyph_idx[0]) / psi.shape[0] * box_mpc
    glyph_max_len = 0.85 * glyph_pitch
    axes[3].add_collection(
        _component_segments(
            psi, extent=extent, component="real",
            scale=component_scale, max_len=glyph_max_len, n_ticks=glyph_n_ticks,
        )
    )
    axes[3].add_collection(
        _component_segments(
            psi, extent=extent, component="imag",
            scale=component_scale, max_len=glyph_max_len, n_ticks=glyph_n_ticks,
        )
    )
    axes[3].add_patch(Circle((peak_x, peak_y), 17.5, fill=False,
                             edgecolor="0.12", linewidth=0.85, zorder=7))
    _panel_label(axes[3], "D  Weyl Ψ₀", "Re/Im Ψ₀ component samples")
    cb = fig.colorbar(im_psi, cax=caxes[3], orientation="horizontal")
    cb.set_label("|Ψ₀| [10⁻⁸ Mpc⁻²]", fontsize=13.5, labelpad=2.5)
    cb.set_ticks([0, 2, 4, 6, 8])
    cb.ax.tick_params(labelsize=12.0, length=3.2, pad=2.0)

    for i, ax in enumerate(axes):
        ax.set_aspect("equal")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xticks([-100, 0, 100])
        ax.set_yticks([-100, 0, 100])
        ax.tick_params(labelsize=12.5)
        ax.set_xlabel("x [Mpc]", labelpad=2, fontsize=14.0)
        ax.grid(False)
        if i == 0:
            ax.set_ylabel("y [Mpc]", labelpad=2, fontsize=14.0)
        else:
            ax.set_yticklabels([])

    # Flow arrows between panels were removed: the formula caption at the
    # bottom of each panel (e.g.~"Phi(k) = -A delta_m(k)/k^2", "Phi_00 ~= -A delta_m",
    # "Re/Im Psi_0 component samples") already encodes the full transform chain.

    print(
        f"[driving-field] {meta.source}, z={meta.redshift:.2f}, "
        f"h={meta.h:.4f}, ombh2={meta.omega_b_h2:.5f}, "
        f"omch2={meta.omega_c_h2:.5f}, ns={meta.n_s:.2f}, "
        f"sigma8={meta.sigma8:.2f}, "
        f"A={float(fields['poisson_a']):.3e} Mpc^-2, seed={seed}"
    )

    pdf = out_stem.with_suffix(".pdf")
    png = out_stem.with_suffix(".png")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-pix", type=int, default=192)
    parser.add_argument("--box-mpc", type=float, default=300.0)
    parser.add_argument("--redshift", type=float, default=0.7)
    parser.add_argument("--out-stem", type=Path, default=FIGURES / "driving_field_portrait")
    args = parser.parse_args()
    render(
        args.seed,
        args.out_stem,
        n_pix=args.n_pix,
        box_mpc=args.box_mpc,
        redshift=args.redshift,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
