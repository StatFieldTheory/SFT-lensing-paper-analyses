#!/usr/bin/env python
"""Phase-0 figure regeneration harness: import each paper generator, rebind its
output constants into the scratch dir (never into the repo), optionally rebind
its FK input to the converged product the deployed figure used, and run it.

One figure per process (the two _plot_style.py copies must not share a process).

    python p0_regen.py --figure KEY [--mode standalone|rebound|cache|valfig]
"""
import argparse, importlib.util, os, shutil, sys
from pathlib import Path
os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np

REPO = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing")
SCR = Path(os.environ.get("P0_SCRATCH", str(Path(__file__).resolve().parent / "p0")))
SSA = REPO / "SFT-lensing-paper-analyses"
SACHS = SSA / "sachs_sft"
A1 = SACHS / "analyses" / "analysis1"
A3 = SACHS / "analyses" / "analysis3"
MC = SACHS / "analyses" / "mc_sachs_2pt"
ETL = SACHS / "callables" / "kappa3_vertex" / "equal_time_limber"
FIGGEN = SSA / "figures"
DFE = REPO / "driver_field_emulators"
FK_CONV = DFE / "products" / "table_permclosed_cut15360_permfix_xi.npz"
ZETA_CONV = DFE / "products" / "zeta_bands_cut15360_gmax85.npz"
FK_MARKERS = DFE / "code" / "mc_fk_complete" / "_markers_pooled.npz"


def load(directory, filename, name=None):
    directory = Path(directory)
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
    name = name or Path(filename).stem
    spec = importlib.util.spec_from_file_location(name, directory / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def out(sub):
    p = SCR / sub
    p.mkdir(parents=True, exist_ok=True)
    return p


def fig_analysis1(mode):
    m = load(A1, "plot_analysis1_O0_vs_pyccl_fkem.py")
    m.OUT_STEM = out("analysis1") / "analysis1_O0_vs_pyccl_fkem"
    m.main()
    return [m.OUT_STEM.with_suffix(".pdf")]


def fig_nlo(mode):
    m = load(A3, "plot_analysis3_nlo_decomposition.py")
    if mode == "rebound":
        m.FK_NPZ = FK_CONV
    m.OUT_STEM = out(f"nlo_{mode}") / "analysis3_nlo_O0_FF_FK"
    print(f"[harness] FK_NPZ = {m.FK_NPZ}")
    m.main()
    return [m.OUT_STEM.with_suffix(".pdf")]


def fig_cl(mode):
    m = load(A3, "plot_analysis3_cl_decomposition.py")
    if mode == "rebound":
        m.FK_NPZ = FK_CONV
    m.OUT_STEM = out(f"cl_{mode}") / "analysis3_cl_O0_FF_FK"
    print(f"[harness] FK_NPZ = {m.FK_NPZ}")
    m.main()
    return [m.OUT_STEM.with_suffix(".pdf")]


def fig_eb(mode):
    C = load(A3, "plot_analysis3_cl_decomposition.py")
    if mode == "rebound":
        C.FK_NPZ = FK_CONV
    m = load(A3, "plot_cl_EB_polarization.py")
    assert m.C is C
    m.OUT_STEM = out(f"eb_{mode}") / "cl_EB_polarization"
    print(f"[harness] FK_NPZ = {C.FK_NPZ}")
    m.main()
    return [m.OUT_STEM.with_suffix(".pdf")]


def fig_multiz(mode):
    m = load(A3, "plot_multiz_kappa_2x5.py")
    m.OUT_STEM = out("multiz") / "multiz_kappa_xi_cl_2x5"
    print(f"[harness] NPZ = {m.NPZ}")
    m.main()
    return [m.OUT_STEM.with_suffix(".pdf")]


def fig_mc(mode):
    """cache: fig_xi_channels --from-cache (replot). valfig: the make_val_figure.py
    route (FK line rebound to the converged fold, FK markers from mc_fk_complete),
    replicated here so the cache rewrite and the figure land in scratch."""
    d = out(f"mc_{mode}")
    (d / "outputs").mkdir(exist_ok=True)
    (d / "figures").mkdir(exist_ok=True)
    shutil.copy(MC / "outputs" / "appendix_mc_curve.npz", d / "outputs" / "appendix_mc_curve.npz")
    m = load(MC, "fig_xi_channels.py")
    m.__file__ = str(d / "fig_xi_channels.py")   # redirects _cache_path() and _plot() output
    if mode == "cache":
        m.main(48_000, 4000, True)
        return [d / "figures" / "xi_kappa_channels.pdf"]
    import math
    original = m._load_kk

    def patched(name):
        if name != "C_corr_op_K_limber_FK":
            return original(name)
        dd = np.load(FK_CONV, allow_pickle=True)
        mask = (dd["a"] == 0) & (dd["b"] == 0) & (dd["order"] == 2)
        gs = []
        for xi, yi in zip(dd["x"][mask], dd["y"][mask]):
            xi = np.asarray(xi, float); yi = np.asarray(yi, float)
            gs.append(math.degrees(math.acos(float(np.clip(
                np.dot(xi / np.linalg.norm(xi), yi / np.linalg.norm(yi)), -1, 1)))) * 60)
        order = np.argsort(gs)
        return np.asarray(gs)[order], np.asarray(dd["value"], float)[mask][order]

    m._load_kk = patched
    g, o0 = m._load_kk("C_corr_op_O0")
    _, ff_mom = m._load_kk("C_corr_op_K_limber_FF")
    _, fk3 = m._load_kk("C_corr_op_K_limber_FK")
    cache = dict(np.load(m._cache_path(), allow_pickle=True))
    g_c = np.asarray(cache["g_mc"], float); ff_c = np.asarray(cache["ffc_mc"], float); se_c = np.asarray(cache["ffc_se"], float)
    g_mc, ffc_mc, ffc_se = [], [], []
    for gv in m.GAMMA_MC:
        j = int(np.argmin(np.abs(g_c - gv)))
        assert abs(g_c[j] - gv) < 1e-6
        g_mc.append(gv); ffc_mc.append(ff_c[j]); ffc_se.append(se_c[j])
    D = dict(g=g, o0=o0, ff_mom=ff_mom, fk3=fk3, g_mc=np.array(g_mc),
             ffc_mc=np.array(ffc_mc), ffc_se=np.array(ffc_se))
    mk = np.load(FK_MARKERS)
    D.update(g_fk=np.asarray(mk["gamma"], float), fk_mc=np.asarray(mk["fk"], float),
             fk_se=np.asarray(mk["err"], float))
    np.savez(m._cache_path(), **D)
    m._plot(D)
    return [d / "figures" / "xi_kappa_channels.pdf"]


def fig_zeta(mode):
    m = load(ETL, "plot_zeta_figures.py")
    m._FIGDIR = out(f"zeta_{mode}")
    if mode == "rebound":
        m._NPZ = ZETA_CONV
    print(f"[harness] _NPZ = {m._NPZ}")
    m.main()
    return [m._FIGDIR / "zeta_driving_field_slices_draft.pdf"]


def fig_ops(mode):
    m = load(FIGGEN, "make_response_corr_operator_figures.py")
    m.FIG_DIR = out("ops")
    m.main()
    return [m.FIG_DIR / "response_operator_draft.pdf", m.FIG_DIR / "corr_operator_slices_draft.pdf"]


def fig_screen(mode):
    m = load(FIGGEN, "make_screen_basis_spin2_figure.py")
    m.FIG_DIR = out("screen")
    m.render()
    return [m.FIG_DIR / "screen_basis_spin2_pattern.pdf"]


def fig_null(mode):
    m = load(FIGGEN, "make_null_congruence_schematic.py")
    d = out("null")
    m.OUT_PDF = d / "null_geodesic_congruence.pdf"
    m.OUT_PNG = d / "null_geodesic_congruence.png"
    m.main()
    return [m.OUT_PDF]


FIGS = {"analysis1": fig_analysis1, "nlo": fig_nlo, "cl": fig_cl, "eb": fig_eb,
        "multiz": fig_multiz, "mc": fig_mc, "zeta": fig_zeta, "ops": fig_ops,
        "screen": fig_screen, "null": fig_null}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True, choices=sorted(FIGS))
    ap.add_argument("--mode", default="standalone")
    a = ap.parse_args()
    outs = FIGS[a.figure](a.mode)
    for p in outs:
        print(f"[harness] OUTPUT {p} exists={Path(p).exists()}")
