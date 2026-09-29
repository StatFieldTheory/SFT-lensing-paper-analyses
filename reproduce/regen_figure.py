#!/usr/bin/env python
"""Regenerate the manuscript's figures from the analysis package, one command each.

Every generator writes into its own ``outputs/`` folder; nothing here touches the
repository ``figures/``. Compare the results with ``check_figures.py`` and copy
them into place with ``deploy.py``.

    python reproduce/regen_figure.py --list
    python reproduce/regen_figure.py --figure 2 3 4
    python reproduce/regen_figure.py --all

Interpreters (override with environment variables): ``STF_PYCCL`` (the PyCCL
conda env), ``STF_SFTWICK`` (the sft-wick conda env), ``CANOES_ROOT`` (the
canoes checkout, needed by figures 6, 8 and 9 through the corr_op callable).
The figures that reproduce byte-for-byte do so only in the interpreter listed
for them (the two envs differ at the 1 pt bounding-box level), see REPRODUCE.md.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from product_paths import active_manifest, resolve_product

PKG = Path(__file__).resolve().parents[1]
PAPER = PKG.parent
PYCCL = os.environ.get("STF_PYCCL", "/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python")
SFTW = os.environ.get("STF_SFTWICK", "/opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python")
CANOES = os.environ.get("CANOES_ROOT", "/Users/zzhang/projects/angular_statistics/canoes")

A3 = "sachs_sft/analyses/analysis3"
MC = "sachs_sft/analyses/mc_sachs_2pt"
ETL = "sachs_sft/callables/kappa3_vertex/equal_time_limber"
RB = "sachs_sft/callables/kappa3_vertex/rebuild"

#: figure number -> (generator, interpreter key, extra args, needs canoes, output(s), paper file(s))
FIGURES = {
    1: (f"sachs_sft/analyses/analysis1/plot_analysis1_O0_vs_pyccl_fkem.py", "sftw", [], False,
        ["sachs_sft/analyses/analysis1/outputs/analysis1_O0_vs_pyccl_fkem.pdf"], ["analysis1_O0_vs_pyccl.pdf"]),
    2: (f"{A3}/plot_analysis3_nlo_decomposition.py", "pyccl", [], False,
        [f"{A3}/outputs/analysis3_nlo_O0_FF_FK.pdf"], ["analysis3_NLO_FFFK.pdf"]),
    3: (f"{A3}/plot_analysis3_cl_decomposition.py", "pyccl", [], False,
        [f"{A3}/outputs/analysis3_cl_O0_FF_FK.pdf"], ["analysis3_cl_full_vs_O0.pdf"]),
    4: (f"{A3}/plot_cl_EB_polarization.py", "pyccl", [], False,
        [f"{A3}/outputs/cl_EB_polarization.pdf"], ["cl_EB_polarization.pdf"]),
    5: (f"{A3}/plot_multiz_kappa_2x5.py", "pyccl", [], False,
        [f"{A3}/outputs/multiz_kappa_xi_cl_2x5.pdf"], ["multiz_kappa_xi_cl.pdf"]),
    6: (f"{MC}/fig_xi_channels.py", "pyccl", ["--from-cache"], True,
        [f"{MC}/figures/xi_kappa_channels.pdf"], ["appendix_mc_workflow.pdf"]),
    7: (f"{ETL}/plot_zeta_figures.py", "pyccl", [], False,
        [f"{ETL}/outputs/zeta_driving_field_slices_draft.pdf"], ["zeta_driving_field_slices_draft.pdf"]),
    8: ("figures/make_response_corr_operator_figures.py", "pyccl", [], True,
        ["figures/outputs/corr_operator_slices_draft.pdf", "figures/outputs/response_operator_draft.pdf"],
        ["corr_operator_slices_draft.pdf", "response_operator_draft.pdf"]),
    10: ("figures/make_screen_basis_spin2_figure.py", "pyccl", [], False,
         ["figures/outputs/screen_basis_spin2_pattern.pdf"], ["screen_basis_spin2_pattern.pdf"]),
    11: ("figures/plot_driving_field_portrait.py", "pyccl", ["--seed", "42"], False,
         ["figures/outputs/driving_field_portrait.pdf"], ["driving_field_portrait.pdf"]),
    12: ("figures/make_selection_rule_schematic.py", "pyccl", ["--seed", "42"], False,
         ["figures/outputs/selection_rule_schematic.pdf"], ["selection_rule_schematic.pdf"]),
    13: ("figures/make_null_congruence_schematic.py", "pyccl", [], False,
         ["figures/outputs/null_geodesic_congruence.pdf"], ["null_geodesic_congruence.pdf"]),
    14: ("figures/make_cosmology_coord_maps.py", "pyccl", [], False,
         ["figures/outputs/cosmology_coord_maps.pdf"], ["cosmology_coord_maps.pdf"]),
    15: ("figures/make_feyndiag.py", "sftw", [], False,
         ["figures/outputs/FeynDiag_2pt_order2.pdf", "figures/outputs/FeynDiag_3pt_order1.pdf"],
         ["FeynDiag_2pt_order2.pdf", "FeynDiag_3pt_order1.pdf"]),
    17: (f"{RB}/fig_cutoff_paper.py", "pyccl", [], False,
         [f"{RB}/outputs/fk_cutoff_convergence.pdf"], ["fk_cutoff_convergence.pdf"]),
}
#: figures 9 and 16 are the second outputs of rows 8 and 15.
ALIASES = {9: 8, 16: 15}

if active_manifest() is not None:
    aligned = "sachs_sft/analyses/analysis1/r1_aligned"
    FIGURES[1] = (
        f"{aligned}/plot_comparison.py", "pyccl",
        ["--sft", str(resolve_product("order0_comparison", Path())),
         "--reference", str(resolve_product("pyccl_reference", Path())),
         "--output-dir", str(PKG / aligned / "outputs/comparison_final"),
         "--archive-existing"], False,
        [f"{aligned}/outputs/comparison_final/analysis1_O0_vs_pyccl_aligned.pdf"],
        ["analysis1_O0_vs_pyccl.pdf"],
    )


def run_one(n: int) -> tuple[int, float]:
    n = ALIASES.get(n, n)
    gen, key, args, canoes, outs, _ = FIGURES[n]
    py = PYCCL if key == "pyccl" else SFTW
    env = dict(os.environ, MPLBACKEND="Agg")
    if canoes:
        env["PYTHONPATH"] = f"{CANOES}/src" + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    script = PKG / gen
    for o in outs:
        (PKG / o).parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    r = subprocess.run([py, str(script), *args], cwd=str(script.parent), env=env,
                       capture_output=True, text=True)
    dt = time.perf_counter() - t0
    log = PKG / "reproduce" / "logs"
    log.mkdir(exist_ok=True)
    (log / f"fig{n:02d}.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    return r.returncode, dt


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", type=int, nargs="+")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for n, (gen, key, args, canoes, outs, papers) in FIGURES.items():
            print(f"{n:2d}  {gen} {' '.join(args)}  [{key}{' +canoes' if canoes else ''}]  -> {', '.join(papers)}")
        return 0
    todo = sorted(FIGURES) if a.all else sorted({ALIASES.get(n, n) for n in (a.figure or [])})
    if not todo:
        ap.error("give --figure N ... or --all")
    bad = 0
    for n in todo:
        rc, dt = run_one(n)
        status = "ok" if rc == 0 else f"FAILED rc={rc}"
        print(f"[fig {n:2d}] {status} in {dt:6.1f} s  ({FIGURES[n][0]})", flush=True)
        bad += rc != 0
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
