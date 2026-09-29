"""Plot signed, input-aligned lensing functions and save stable residual metrics."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _plot_style import (  # noqa: E402
    PALETTE, apply_rcparams, plot_signed_line, plot_signed_markers,
)

apply_rcparams()
import matplotlib.pyplot as plt  # noqa: E402

CHANNELS = (
    ("kk", r"$\xi_{\kappa}$"),
    ("xip", r"$\xi_{+}$"),
    ("xim", r"$\xi_{-}$"),
    ("kgt", r"$\xi_{\kappa\gamma_t}$"),
)


def load_functions(path: Path) -> dict[str, np.ndarray]:
    """Read numeric arrays only, rejecting ambiguous grids or nonfinite values."""
    with np.load(path, allow_pickle=False) as saved:
        result = {key: np.asarray(saved[key], dtype=float) for key in
                  ("theta_arcmin", *(item[0] for item in CHANNELS))}
    theta = result["theta_arcmin"]
    if theta.ndim != 1 or len(theta) < 2 or np.any(theta <= 0):
        raise ValueError(f"Invalid angular grid in {path}")
    if np.any(np.diff(theta) <= 0):
        raise ValueError(f"Angles must be unique and increasing in {path}")
    for key, value in result.items():
        if value.shape != theta.shape or not np.all(np.isfinite(value)):
            raise ValueError(f"Invalid {key} in {path}")
    return result


def zero_crossings(theta: np.ndarray, value: np.ndarray) -> list[float]:
    """Estimate crossings by linear interpolation of signed values in log theta."""
    zeros = theta[value == 0].tolist()
    crossing = np.flatnonzero(np.signbit(value[:-1]) != np.signbit(value[1:]))
    for index in crossing:
        if value[index] == 0 or value[index + 1] == 0:
            continue
        fraction = -value[index] / (value[index + 1] - value[index])
        zeros.append(float(np.exp(np.log(theta[index]) + fraction *
                                  np.log(theta[index + 1] / theta[index]))))
    return sorted(zeros)


def channel_metrics(theta: np.ndarray, actual: np.ndarray,
                    reference: np.ndarray, relative_floor: float) -> dict:
    """Quantify signed differences without ratios at nearly zero reference signal."""
    peak = float(np.max(np.abs(reference)))
    if peak == 0:
        raise ValueError("Cannot normalize a channel with identically zero reference")
    delta = actual - reference
    mask = np.abs(reference) >= relative_floor * peak
    relative = np.abs(delta[mask] / reference[mask])
    actual_zeros = zero_crossings(theta, actual)
    reference_zeros = zero_crossings(theta, reference)
    maximum = int(np.argmax(np.abs(delta)))
    return {
        "reference_peak_abs": peak,
        "max_abs_difference": float(np.max(np.abs(delta))),
        "max_abs_difference_theta_arcmin": float(theta[maximum]),
        "max_abs_difference_over_peak_percent": float(100 * np.max(np.abs(delta)) / peak),
        "rms_difference_over_peak_percent": float(100 * np.sqrt(np.mean(delta**2)) / peak),
        "relative_mask_floor_fraction_of_peak": relative_floor,
        "relative_mask_included_samples": int(mask.sum()),
        "relative_mask_excluded_fraction": float(1 - mask.mean()),
        "max_abs_masked_relative_percent": float(100 * relative.max()) if relative.size else None,
        "median_abs_masked_relative_percent": float(100 * np.median(relative)) if relative.size else None,
        "sft_zero_crossings_arcmin": actual_zeros,
        "reference_zero_crossings_arcmin": reference_zeros,
        "zero_crossing_shifts_arcmin": [left - right for left, right in
                                         zip(actual_zeros, reference_zeros)]
        if len(actual_zeros) == len(reference_zeros) else None,
        "zero_crossing_method": "Linear signed interpolation in log angle on the displayed samples",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sft", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--theta-max", type=float, default=2000.0)
    parser.add_argument("--relative-floor", type=float, default=0.01)
    parser.add_argument("--archive-existing", action="store_true",
                        help="Move existing plot artifacts into a dated archive before writing")
    args = parser.parse_args()
    if not 0 < args.relative_floor <= 1 or not np.isfinite(args.theta_max):
        parser.error("Require 0 < relative-floor <= 1 and finite theta-max")
    actual = load_functions(args.sft)
    reference = load_functions(args.reference)
    theta = actual["theta_arcmin"]
    if theta.shape != reference["theta_arcmin"].shape or not np.allclose(
        theta, reference["theta_arcmin"], rtol=1e-8, atol=1e-7,
    ):
        raise ValueError("Use identical angular samples, rather than interpolating one result")
    mask = (theta >= 0.5 - 1e-6) & (theta <= args.theta_max)
    if mask.sum() < 2:
        raise ValueError("Displayed interval contains fewer than two samples")
    theta = theta[mask]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    artifact_names = (
        "analysis1_O0_vs_pyccl_aligned.pdf", "analysis1_O0_vs_pyccl_aligned.png",
        "residuals.pdf", "residuals.png", "comparison.json",
    )
    previous = [args.output_dir / name for name in artifact_names
                if (args.output_dir / name).exists()]
    if previous:
        if not args.archive_existing:
            raise FileExistsError("Plot artifacts exist. Use a fresh directory or --archive-existing.")
        archive = args.output_dir / datetime.now(timezone.utc).strftime("archive_%Y-%m-%dT%H%M%S_%fZ")
        archive.mkdir()
        for path in previous:
            path.rename(archive / path.name)
    metadata = {
        "sft": str(args.sft.resolve()),
        "reference": str(args.reference.resolve()),
        "sft_sha256": hashlib.sha256(args.sft.read_bytes()).hexdigest(),
        "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        "theta_arcmin": theta.tolist(),
        "rms_weighting": "Equal weights for the displayed logarithmic angular samples",
        "channels": {},
    }
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.4), sharex=True, sharey=True)
    diagnostic, residual_axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for index, ((key, label), ax, residual_ax) in enumerate(
        zip(CHANNELS, axes.flat, residual_axes.flat)
    ):
        left, right = actual[key][mask], reference[key][mask]
        metrics = channel_metrics(theta, left, right, args.relative_floor)
        metadata["channels"][key] = metrics
        plot_signed_line(ax, theta, right, label="PyCCL", sign_marker_size=3.5)
        plot_signed_markers(ax, theta, left, color=PALETTE[index], marker="o", label="SFT Order-0")
        ax.set_title(label)
        ax.set_xlim(0.5, args.theta_max)
        ax.set_ylim(1e-9, 2e-3)
        ax.legend(loc="best")
        residual_ax.semilogx(theta, 100 * (left - right) / metrics["reference_peak_abs"],
                            color=PALETTE[index], marker="o", markersize=3)
        residual_ax.axhline(0, color="0.5", lw=0.7)
        residual_ax.set_title(label)
        residual_ax.set_ylabel(r"$100\,\Delta\xi / \max|\xi_{\rm PyCCL}|$")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\gamma$ [arcmin]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|\xi|$")
    for ax in residual_axes[-1]:
        ax.set_xlabel(r"$\gamma$ [arcmin]")
    fig.tight_layout()
    diagnostic.tight_layout()
    pdf_metadata = {"CreationDate": None, "ModDate": None}
    fig.savefig(args.output_dir / "analysis1_O0_vs_pyccl_aligned.pdf", metadata=pdf_metadata)
    fig.savefig(args.output_dir / "analysis1_O0_vs_pyccl_aligned.png")
    diagnostic.savefig(args.output_dir / "residuals.pdf", metadata=pdf_metadata)
    diagnostic.savefig(args.output_dir / "residuals.png")
    plt.close(fig)
    plt.close(diagnostic)
    (args.output_dir / "comparison.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps(metadata["channels"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
