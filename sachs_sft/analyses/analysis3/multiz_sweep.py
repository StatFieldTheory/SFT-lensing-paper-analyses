"""Multi-redshift xi_kappa(gamma) sweep helper: z -> lambda mapping and the per-plane
sft-wick sweep used by compute_multiz_kappa_2pcf.py.

Moved here 2026-09-02 from talk/scripts/run_multiz_components.py (written 2026-06-10 for the
conference-talk animation); the talk keeps a thin wrapper importing this module. Original docstring:

Compute multi-redshift xi_kappa(gamma) decomposition for the conference talk gif.

Produces ``talk/assets/figures/_data/multiz_components.npz`` holding the
convergence two-point correlation xi_kappa(gamma), decomposed into

    O0  -- order-0 (textbook) term            (vertex_type = none, orders = [0])
    FF  -- order-2 nonlinear-propagation term  (vertex_type = F)
    FK  -- order-2 kappa3 vertex term          (vertex_type = FK)

on ~20 source redshifts z_s in [0.5, 5], reusing the validated analysis-3
sft-wick configs under ``SFT-lensing-paper-analyses/sachs_sft/sftwick_outputs/2PCF/``. Only the
[0, 0] (convergence auto) component pair is extracted.

The three configs (C_corr_op_O0, C_corr_op_K_limber_FF, C_corr_op_K_limber_FK)
reference all z-independent callable tables BY PATH; the only per-z change is the
sweep ``t_final_grid`` (= list of source affine distances lambda(z_s)). Adding
more t_final values is cheap: the propagator tables are built up to t_max=2360
independent of t_final, so a single sweep over all 20 lambda values reuses one
propagator build per vertex_type.

z -> lambda mapping uses the project closed-form background
``SFT-lensing-paper-analyses/sachs_sft/scripts/D_callable.py`` (d lambda / d chi = a^2,
D = a*chi). lambda(z=1) = 1822.721 matches the existing fk_multiz cache.

NOTE on the z=5 convention: D_callable's self-consistent lambda(5) = 2317.96,
but the production analysis-3 baseline (and the boundary-check npz
mc_vs_analysis3_xi_kappa.npz) was computed at the hard-coded t_final = 2313.029.
For the boundary check (--boundary-check) we therefore append the baseline
lambda value as an extra plane and compare it against the archived analysis-3
slice. The 20-point gif grid itself uses the self-consistent mapping so the
redshift evolution is smooth and physical.

The output npz emitted by sft-wick stores the x/y direction cosines as object
arrays (tuples), so reading it back requires allow_pickle=True. The file is the
pipeline's own freshly written /tmp output, not untrusted external content.

Run (sft-wick env), SMOKE first then full:
    .../envs/sft-wick/bin/python run_multiz_components.py --smoke
    .../envs/sft-wick/bin/python run_multiz_components.py --boundary-check
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Fixed project paths.
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[4]
SACHS = REPO / "SFT-lensing-paper-analyses" / "sachs_sft"
TWO_PCF = SACHS / "sftwick_outputs" / "2PCF"

CONFIGS = {
    "o0": TWO_PCF / "C_corr_op_O0" / "config_L2.yaml",
    "ff": TWO_PCF / "C_corr_op_K_limber_FF" / "config_L2.yaml",
    "fk": TWO_PCF / "C_corr_op_K_limber_FK" / "config_L2.yaml",
}

OUT_NPZ = Path(__file__).resolve().parent / "outputs" / "multiz_components_20planes.npz"
BOUNDARY_NPZ = (
    SACHS
    / "analyses"
    / "mc_sachs_2pt"
    / "figures"
    / "_archive_2026-06-06"
    / "mc_vs_analysis3_xi_kappa.npz"
)

# Production analysis-3 baseline source distance (z_s = 5 plane in the archived
# boundary-check npz). Used only for the boundary comparison, not the gif grid.
LAM_BASELINE_Z5 = 2313.0288751857356

# Scratch dir for per-vertex output npz emitted by run_workflow.
SCRATCH = Path(__file__).resolve().parent / "outputs" / "_multiz_scratch"


# ---------------------------------------------------------------------------
# z -> lambda mapping (project closed-form background).
# ---------------------------------------------------------------------------
def build_lambda_of_z():
    """Return a CubicSpline lambda(z) consistent with D_callable's background."""
    sys.path.insert(0, str(SACHS / "scripts"))
    import D_callable as Dc  # noqa: E402
    from scipy.interpolate import CubicSpline  # noqa: E402

    z_tab = Dc._CACHE["z"]
    lam_tab = Dc._CACHE["lambda"]
    return CubicSpline(z_tab, lam_tab)


def make_z_grid(n_z: int) -> np.ndarray:
    return np.linspace(0.5, 5.0, n_z)


# ---------------------------------------------------------------------------
# Single-vertex sweep over a list of lambda (= t_final) values.
# ---------------------------------------------------------------------------
def run_vertex_sweep(
    config_path: Path,
    lam_values: np.ndarray,
    label: str,
    n_jobs: int = 4,
    vertex_types=None,
):
    """Run one config over all lambda values; return (gamma[40], xi[n_lam, 40]).

    The sft-wick driver must chdir to the YAML's parent (callable-module paths
    are YAML-relative; here they are absolute, but cache_path / output paths are
    CWD-relative in general -- we keep the chdir for robustness). We point the
    output at a private /tmp npz and reuse the existing propagator/expand caches
    referenced (by absolute path) inside each config.
    """
    from sft_wick.workflow.config import (  # noqa: E402
        load_workflow_config,
        run_workflow,
    )

    here = config_path.parent
    prev_cwd = Path.cwd()
    os.chdir(here)
    try:
        cfg = load_workflow_config(config_path)

        out_path = SCRATCH / f"xi_multiz_{label}.npz"
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Cap parallelism at n_jobs (default 4): the YAMLs ship with
        # sweep.n_jobs=-1 (all cores -> 32 loky workers), but each worker
        # re-loads every callable table, so per-worker RAM x n_workers can OOM
        # (project memory: the May-14 OOM-reboot). expand.n_jobs stays 1 so the
        # joblib mutual-exclusion invariant (exactly one of {expand, sweep}
        # n_jobs > 1) holds. propagators.n_jobs is also capped for safety.
        sweep_overrides = dict(
            t_final_grid=[float(v) for v in lam_values],
            n_jobs=int(n_jobs),
        )
        if vertex_types is not None:
            sweep_overrides["vertex_types"] = list(vertex_types)
        new_sweep = dataclasses.replace(cfg.sweep, **sweep_overrides)
        new_props = dataclasses.replace(cfg.propagators, n_jobs=int(n_jobs))
        new_output = [dataclasses.replace(cfg.output[0], path=str(out_path))]
        cfg = dataclasses.replace(
            cfg, sweep=new_sweep, propagators=new_props, output=new_output
        )

        print(
            f"[{label}] vertex_types={cfg.sweep.vertex_types} "
            f"orders={cfg.expand.orders} n_lambda={len(lam_values)} "
            f"sweep.n_jobs={cfg.sweep.n_jobs} "
            f"props.n_jobs={cfg.propagators.n_jobs}",
            flush=True,
        )
        t0 = time.perf_counter()
        run_workflow(cfg)
        dt = time.perf_counter() - t0
        print(f"[{label}] run_workflow done in {dt:.1f} s", flush=True)
    finally:
        os.chdir(prev_cwd)

    return extract_kk(out_path, lam_values)


def extract_kk(out_path: Path, lam_values: np.ndarray):
    """Extract the [0,0] xi_kappa(gamma) per t_final from an emitted npz.

    The emitted record arrays carry one row per (gamma, component_pair,
    t_final). We mask a==0 & b==0, then group by t_final preserving gamma order.
    gamma is recovered from the y direction-cosine (angle from x = [0,0,1]).
    allow_pickle is required because x/y are stored as object tuples; the file is
    our own freshly written /tmp output.
    """
    d = np.load(out_path, allow_pickle=True)
    a = d["a"]
    b = d["b"]
    t_final = d["t_final"]
    value = d["value"]
    y = d["y"]

    mask = (a == 0) & (b == 0)
    a_t = t_final[mask]
    a_v = value[mask]
    a_y = y[mask]

    lam_sorted = np.array([float(v) for v in lam_values])
    first_lam = lam_sorted[0]
    blk0 = np.isclose(a_t, first_lam, rtol=0, atol=1e-6)
    n_gamma = int(blk0.sum())
    if n_gamma == 0:
        raise RuntimeError(
            f"no [0,0] rows at first t_final={first_lam} in {out_path}"
        )

    def gamma_of_y(yvec) -> float:
        nx = np.array([0.0, 0.0, 1.0])
        ny = np.array(list(yvec), dtype=float)
        cos = float(np.clip(np.dot(nx, ny) / np.linalg.norm(ny), -1.0, 1.0))
        return np.degrees(np.arccos(cos)) * 60.0  # arcmin

    gamma = np.array([gamma_of_y(yv) for yv in a_y[blk0]])

    xi = np.zeros((len(lam_sorted), n_gamma), dtype=float)
    for i, lam in enumerate(lam_sorted):
        sel = np.isclose(a_t, lam, rtol=0, atol=1e-6)
        if sel.sum() != n_gamma:
            raise RuntimeError(
                f"t_final={lam}: got {int(sel.sum())} rows, expected {n_gamma}"
            )
        xi[i] = a_v[sel]
    return gamma, xi


# ---------------------------------------------------------------------------
# Boundary check against the archived analysis-3 z=5 slice.
# ---------------------------------------------------------------------------
def boundary_check(comp) -> None:
    if not BOUNDARY_NPZ.exists():
        print(f"[boundary] MISSING {BOUNDARY_NPZ}; skipping", flush=True)
        return
    ref = np.load(BOUNDARY_NPZ)
    print(f"[boundary] ref keys: {ref.files}", flush=True)
    for key in ("o0", "ff", "fk"):
        ours = comp[f"{key}_baseline"]
        theirs = ref[key]
        denom = np.maximum(np.abs(theirs), 1e-30)
        rel = np.abs(ours - theirs) / denom
        print(
            f"[boundary] {key}: max rel-diff = {rel.max():.3e} "
            f"(median {np.median(rel):.3e}); "
            f"ours[0]={ours[0]:.6e} theirs[0]={theirs[0]:.6e}",
            flush=True,
        )


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="run only 3 z values to validate end-to-end",
    )
    ap.add_argument(
        "--n-z", type=int, default=20, help="number of source redshifts"
    )
    ap.add_argument(
        "--boundary-check",
        action="store_true",
        help="append the baseline lambda(z=5)=2313.029 plane and compare to "
        "the archived analysis-3 slice",
    )
    ap.add_argument(
        "--no-save",
        action="store_true",
        help="do not write the output npz (smoke/diagnostic only)",
    )
    args = ap.parse_args()

    lam_of_z = build_lambda_of_z()

    if args.smoke:
        z_grid = np.array([0.5, 2.75, 5.0])
    else:
        z_grid = make_z_grid(args.n_z)
    lam = np.asarray(lam_of_z(z_grid), dtype=float)

    print("z grid :", np.round(z_grid, 4).tolist(), flush=True)
    print("lambda :", np.round(lam, 4).tolist(), flush=True)
    print(
        f"sanity: lambda(1)={float(lam_of_z(1.0)):.4f} (expect 1822.7), "
        f"lambda(5)={float(lam_of_z(5.0)):.4f} (D_callable self-consistent)",
        flush=True,
    )

    lam_run = lam.copy()
    if args.boundary_check:
        lam_run = np.concatenate([lam, [LAM_BASELINE_Z5]])

    comp = {}
    gamma_ref = None
    for label in ("o0", "ff", "fk"):
        gamma, xi = run_vertex_sweep(CONFIGS[label], lam_run, label)
        if gamma_ref is None:
            gamma_ref = gamma
        else:
            if not np.allclose(gamma, gamma_ref, rtol=1e-6, atol=1e-6):
                raise RuntimeError(f"gamma grid mismatch for {label}")
        comp[label] = xi[: len(lam)]
        if args.boundary_check:
            comp[f"{label}_baseline"] = xi[len(lam)]

    assert gamma_ref is not None
    gamma = gamma_ref

    o0, ff, fk = comp["o0"], comp["ff"], comp["fk"]
    iz = len(z_grid) - 1
    g0 = 0
    print("--- magnitude sanity (highest-z plane, smallest gamma) ---")
    print(f"  z={z_grid[iz]:.3f}  O0={o0[iz, g0]:.4e}  "
          f"FF={ff[iz, g0]:.4e}  FK={fk[iz, g0]:.4e}", flush=True)
    print(f"  |FK|/O0 = {abs(fk[iz, g0]) / o0[iz, g0]:.4f} "
          f"(expect ~0.06 near z=5)", flush=True)
    print(f"  FK/FF   = {fk[iz, g0] / ff[iz, g0]:.1f} (expect FK >> FF)",
          flush=True)
    print(f"  O0 falls with gamma? O0[small]={o0[iz, 0]:.3e} > "
          f"O0[large]={o0[iz, -1]:.3e} : {bool(o0[iz, 0] > o0[iz, -1])}",
          flush=True)

    if args.boundary_check:
        boundary_check(comp)

    if not args.no_save:
        OUT_NPZ.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            OUT_NPZ,
            z=z_grid,
            gamma=gamma,
            o0=o0,
            ff=ff,
            fk=fk,
            lam=lam,
        )
        print(
            f"saved {OUT_NPZ}\n  keys: z{z_grid.shape} gamma{gamma.shape} "
            f"o0{o0.shape} ff{ff.shape} fk{fk.shape} lam{lam.shape}",
            flush=True,
        )


if __name__ == "__main__":
    main()
