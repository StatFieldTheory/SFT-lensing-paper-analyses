#!/usr/bin/env python
"""Fail-loud path patcher for the Phase-2 moves.

Usage: p2_patch.py <group>

Each entry is (relative file, [(old, new), ...]) applied after the git mv of
that group; a missing `old` string aborts the run so nothing rots silently.
Regex entries are marked with a leading "re:" on the old string.
"""
import re
import sys
from pathlib import Path

PKG = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses")
KV = "sachs_sft/callables/kappa3_vertex"
RB = f"{KV}/rebuild"
PA = f"{KV}/equal_time_limber_cut15360_permaware"
ETL = f"{KV}/equal_time_limber"
ROOT5 = "Path(__file__).resolve().parents[5]"

OLD_PREFIX = "/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products"
NEW_PREFIX = f"/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/{RB}"

GROUP_A = {
    f"{RB}/_common.py": [
        ("#: driver_field_emulators/code, which holds the b_model package.\n_CODE_ROOT = _HERE.parent\n"
         "#: repository root, STF_lensing.\n_REPO = _CODE_ROOT.parents[1]\n",
         "#: sachs_sft/callables/kappa3_vertex, which holds the b_model package.\n_CODE_ROOT = _HERE.parent\n"
         "#: repository root, STF_lensing.\n_REPO = _CODE_ROOT.parents[3]\n"),
        ("re:ARCHIVED_L2 = \\(\n.*?\n\\)\n",
         "ARCHIVED_L2 = _ETL / \"inputs\" / \"l2_lambda_grid_z5covgrid16.npz\"\n"),
    ],
    f"{RB}/run_fk_variant.py": [
        ("re:SACHS_SFT = \\(\n    Path\\(__file__\\)\\.resolve\\(\\)\\.parents\\[2\\]\n    / \"SFT-lensing-paper-analyses\" / \"sachs_sft\"\n\\)\n",
         "SACHS_SFT = Path(__file__).resolve().parents[3]\n"),
        ("Usage (from ``driver_field_emulators/code``)::", "Usage (from ``sachs_sft/callables/kappa3_vertex/rebuild``)::"),
        ("$SFT run_fk_variant.py ../products/collapsed_dense_tree_cut1000.npz", "$SFT run_fk_variant.py products/<table>.npz"),
    ],
    f"{RB}/sweep_cutoff.sh": [
        ("P=/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/products\n",
         "P=\"$(cd \"$(dirname \"$0\")\" && pwd)/products\"\n"),
        ("\"$SFT\" ../run_fk_variant.py \"$OUT\"", "\"$SFT\" \"$(dirname \"$0\")/run_fk_variant.py\" \"$OUT\""),
    ],
    f"{RB}/assemble.py": [
        ("producer=\"driver_field_emulators/code/rebuild/assemble.py\"",
         "producer=\"SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/assemble.py\""),
    ],
    f"{PA}/tests/test_perm_aware.py": [
        ("DEPLOYED = (Path(__file__).resolve().parents[4]\n"
         "            / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"callables\"\n"
         "            / \"kappa3_vertex\" / \"equal_time_limber\"\n",
         "DEPLOYED = (Path(__file__).resolve().parents[2]\n"
         "            / \"equal_time_limber\"\n"),
    ],
    f"{PA}/compare_callables.py": [
        ("O0 = (\"/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/\"\n"
         "      \"sachs_sft/sftwick_outputs/2PCF/C_corr_op_O0/xi_C_corr_op_O0.npz\")",
         "O0 = str(Path(__file__).resolve().parents[3] / \"sftwick_outputs\" / \"2PCF\"\n"
         "         / \"C_corr_op_O0\" / \"xi_C_corr_op_O0.npz\")"),
    ],
    f"{ETL}/build_equal_time_limber_table.py": [
        ("re:_ARCHIVED_L2 = \\(\n.*?\n\\)\n",
         "_ARCHIVED_L2 = Path(__file__).resolve().parent / \"inputs\" / \"l2_lambda_grid_z5covgrid16.npz\"\n"),
    ],
    f"{RB}/fig_cutoff_paper.py": [
        ("    ap.add_argument(\"--products\", type=Path, required=True)\n    ap.add_argument(\"--out\", type=Path, required=True)\n",
         "    ap.add_argument(\"--products\", type=Path, default=_LADDER,\n"
         "                    help=\"folder holding the folded cutoff ladder\")\n"
         "    ap.add_argument(\"--out\", type=Path,\n"
         "                    default=Path(__file__).resolve().parent / \"outputs\" / \"fk_cutoff_convergence.pdf\")\n"),
        ("_MC = _REPO / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"analyses\" / \"mc_sachs_2pt\"\n",
         "_MC = _REPO / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"analyses\" / \"mc_sachs_2pt\"\n"
         "_LADDER = (_REPO / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"sftwick_outputs\"\n"
         "           / \"2PCF\" / \"cutoff_ladder\")\n"),
        ("    <PyCCL python> fig_cutoff_paper.py --products ../../products \\\n        --out ../../../figures/fk_cutoff_convergence.pdf\n",
         "    <PyCCL python> fig_cutoff_paper.py            # defaults: the cutoff_ladder run folder, outputs/\n"),
    ],
}

# Plain-string passes applied to every .py/.sh in a folder (after the explicit
# entries). A third element lists file names the pattern must skip.
GROUP_A_GLOBAL = {
    RB: [
        ("Path(__file__).resolve().parents[3]", ROOT5, ("run_fk_variant.py",)),
        ("_REPO / \"driver_field_emulators\" / \"code\" / \"callable_fixed\"",
         "Path(__file__).resolve().parents[1] / \"equal_time_limber_cut15360_permaware\""),
        ("_REPO / \"driver_field_emulators\" / \"products\"", "Path(__file__).resolve().parent / \"products\""),
        ("../../products/", "products/"),
        ("../products/", "products/"),
        ("driver_field_emulators/code/rebuild", "sachs_sft/callables/kappa3_vertex/rebuild"),
        ("driver_field_emulators/code", "sachs_sft/callables/kappa3_vertex"),
    ],
    f"{KV}/b_model": [
        ("driver_field_emulators/code", "sachs_sft/callables/kappa3_vertex"),
    ],
}

GROUP_A_JSON = {  # jobs_*.json: absolute paths of the August builds, prefix-mapped
    f"{RB}": [
        (OLD_PREFIX + "/pieces/", NEW_PREFIX + "/products/pieces/"),
        (OLD_PREFIX + "/logs/", NEW_PREFIX + "/logs/"),
        (OLD_PREFIX + "/", NEW_PREFIX + "/products/"),
    ],
}

GROUP_B = {
    "sachs_sft/scripts/run_FF_single.py": [
        ("    cfg_path = here / cfg_yaml_name\n    os.chdir(here)\n    cfg = load_workflow_config(cfg_path)\n"
         "    _clear_cache(cfg.expand.cache_path, here)\n    _clear_cache(cfg.propagators.cache_path, here)\n",
         "    cfg_path = (here / cfg_yaml_name).resolve()\n"
         "    # Everything relative in the YAML resolves against the YAML's own folder:\n"
         "    # sft-wick resolves module and data paths that way at load time, and its\n"
         "    # output and cache paths are CWD-relative, so the CWD is set to that folder.\n"
         "    base = cfg_path.parent\n    os.chdir(base)\n    cfg = load_workflow_config(cfg_path)\n"
         "    _clear_cache(cfg.expand.cache_path, base)\n    _clear_cache(cfg.propagators.cache_path, base)\n"),
        ("    default_out = yaml_out if yaml_out.is_absolute() else (here / yaml_out).resolve()\n"
         "    if default_out.exists():\n        default_out.unlink()\n",
         "    default_out = yaml_out if yaml_out.is_absolute() else (base / yaml_out).resolve()\n"
         "    if default_out.exists():\n"
         "        # A run folder carrying a PRODUCTION marker holds a deployed result:\n"
         "        # refuse to unlink it unless the caller says so explicitly. A variant\n"
         "        # config that inherited a production output path destroyed the\n"
         "        # deployed FK sweep on 2026-08-25.\n"
         "        if (default_out.parent / \"PRODUCTION\").exists() and not os.environ.get(\"SFT_WICK_FORCE_OVERWRITE\"):\n"
         "            sys.exit(f\"[run_FF_single] refusing to overwrite {default_out}: its folder is marked \"\n"
         "                     \"PRODUCTION; set SFT_WICK_FORCE_OVERWRITE=1 to override\")\n"
         "        default_out.unlink()\n"),
    ],
}


AN = "sachs_sft/analyses"
MCF = f"{AN}/mc_fk_complete"
REV = f"{AN}/revision_2026-08"
AUD = f"{AN}/fk_audit_2026-08"
A3 = f"{AN}/analysis3"
PERM = "callables/kappa3_vertex/equal_time_limber_cut15360_permaware"
FKRUN = "sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"
MULTIZ_EXPR = "_REPO / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"sftwick_outputs\" / \"2PCF\" / \"multiz\""

GROUP_C = {
    f"{MCF}/_bootstrap.py": [
        ("REPO = _HERE.parents[2]\n", "REPO = _HERE.parents[4]\n"),
        ("FIX_DIR = REPO / \"driver_field_emulators\" / \"code\" / \"callable_fixed\"\n",
         "FIX_DIR = _HERE.parents[2] / \"callables\" / \"kappa3_vertex\" / \"equal_time_limber_cut15360_permaware\"\n"),
        ("PRODUCTS = REPO / \"driver_field_emulators\" / \"products\"\n",
         "PRODUCTS = _HERE.parents[2] / \"sftwick_outputs\" / \"2PCF\" / \"C_corr_op_K_limber_FK_cut15360_permfix\"\n"),
        ("TABLE = PRODUCTS / \"table_permclosed_cut15360.npz\"\n", "TABLE = FIX_DIR / \"table_permclosed.npz\"\n"),
        ("FOLD = PRODUCTS / \"table_permclosed_cut15360_permfix_xi.npz\"\n",
         "FOLD = PRODUCTS / \"xi_C_corr_op_K_limber_FK_cut15360_permfix.npz\"\n"),
    ],
    f"{REV}/regenerate.py": [
        ("_OUT = _REPO / \"driver_field_emulators\" / \"figures\" / \"corrected\"\n", "_OUT = _HERE / \"outputs\"\n"),
        ("    python regenerate.py --fk-npz ../../products/<corrected>_xi.npz \\\n",
         "    python regenerate.py --fk-npz ../../" + FKRUN + " \\\n"),
    ],
    f"{REV}/regenerate_multiz.py": [
        ("_TALK = _REPO / \"talk\" / \"scripts\"\n",
         "_TALK = _A3  # multiz_sweep.py lives beside the compute script since 2026-09-02\n"),
        ("                    default=_REPO / \"driver_field_emulators\" / \"products\" / \"_variant_multiz_conv\")\n",
         "                    default=_HERE / \"outputs\" / \"_variant_multiz_conv\")\n"),
        ("    sys.path.insert(0, str(_REPO / \"driver_field_emulators\" / \"code\"))\n",
         "    sys.path.insert(0, str(_REPO / \"SFT-lensing-paper-analyses\" / \"sachs_sft\" / \"callables\" / \"kappa3_vertex\" / \"rebuild\"))\n"),
    ],
    f"{REV}/regenerate_zeta_slices.py": [
        ("_OUT_NPZ = _REPO / \"driver_field_emulators\" / \"products\" / \"zeta_bands_cut15360.npz\"\n",
         "_OUT_NPZ = _ETL / \"outputs\" / \"zeta_bands_cut15360.npz\"\n"),
        ("_FIGDIR = _REPO / \"driver_field_emulators\" / \"figures\" / \"corrected\"\n", "_FIGDIR = _HERE / \"outputs\"\n"),
    ],
    f"{REV}/run_ff_one_lambda.py": [
        ("PRODUCTS = _REPO / \"driver_field_emulators\" / \"products\"\n", "PRODUCTS = " + MULTIZ_EXPR + "\n"),
    ],
    f"{REV}/run_ff_multiz_serial.py": [
        ("OUT = _REPO / \"driver_field_emulators\" / \"products\" / \"multiz_ff_real.npz\"\n",
         "OUT = " + MULTIZ_EXPR + " / \"multiz_ff_real.npz\"\n"),
    ],
    f"{AUD}/redshift_fold_trends.py": [
        ("PRODUCTS = HERE.parent / \"products\"\n", "PRODUCTS = HERE / \"products\"\n"),
        ("_SCRIPTS = (HERE.parent.parent / \"SFT-lensing-paper-analyses\" / \"sachs_sft\"\n", "_SCRIPTS = (HERE.parents[1]\n"),
        ("MULTIZ = (HERE.parent.parent / \"talk\" / \"assets\" / \"figures\" / \"_data\"\n", "MULTIZ = (HERE.parents[1] / \"analysis3\" / \"inputs\"\n"),
        ("\"multiz_components.npz\"", "\"multiz_components_talk_2026-06-10.npz\""),
    ],
    f"{AUD}/redshift_pershell.py": [
        ("PRODUCTS = HERE.parent / \"products\"\n", "PRODUCTS = HERE / \"products\"\n"),
        ("_SCRIPTS = (HERE.parent.parent / \"SFT-lensing-paper-analyses\" / \"sachs_sft\"\n", "_SCRIPTS = (HERE.parents[1]\n"),
    ],
    f"{AUD}/plot_fk_variants.py": [
        ("PRODUCTS = HERE.parent / \"products\"\n", "PRODUCTS = HERE / \"products\"\n"),
        ("FIGURES = HERE.parent / \"figures\"\n", "FIGURES = HERE / \"figures\"\n"),
    ],
    f"{AUD}/mc_crosscheck_dense.py": [
        ("re:Path\\(__file__\\)\\.resolve\\(\\)\\.parents\\[2\\]\n\\s*/ \"SFT-lensing-paper-analyses\" / \"sachs_sft\"\\)",
         "Path(__file__).resolve().parents[1])"),
    ],
    f"{AUD}/reduced_shear/reduced_shear.py": [
        ("sys.path.insert(0, str(Path(__file__).resolve().parents[1]))\n",
         "sys.path.insert(0, str(Path(__file__).resolve().parents[3] / \"callables\" / \"kappa3_vertex\"))\n"),
    ],
    f"{A3}/compute_multiz_kappa_2pcf.py": [
        ("TALK_SCRIPTS = \"/Users/zzhang/Documents/MyDrafts/STF_lensing/talk/scripts\"\nsys.path.insert(0, TALK_SCRIPTS)\nimport run_multiz_components as R  # noqa: E402\n",
         "sys.path.insert(0, str(HERE))\nimport multiz_sweep as R  # noqa: E402  (moved here from talk/scripts/run_multiz_components.py, 2026-09-02)\n"),
        ("lifting (run_vertex_sweep, z->lambda mapping) is imported from the talk's\n``run_multiz_components.py``.",
         "lifting (run_vertex_sweep, z->lambda mapping) is imported from ``multiz_sweep.py``\nin this folder (the talk repository's ``run_multiz_components.py`` until 2026-09-02)."),
    ],
    f"{A3}/plot_multiz_kappa_2x5.py": [
        ("    _FF_REAL = (Path(__file__).resolve().parents[4] / \"driver_field_emulators\"\n                / \"products\" / \"multiz_ff_real_all5.npz\")\n",
         "    _FF_REAL = (Path(__file__).resolve().parents[2] / \"sftwick_outputs\" / \"2PCF\"\n                / \"multiz\" / \"multiz_ff_real_all5.npz\")\n"),
    ],
    f"{A3}/multiz_sweep.py": [
        ("REPO = Path(\"/Users/zzhang/Documents/MyDrafts/STF_lensing\")\n", "REPO = Path(__file__).resolve().parents[4]\n"),
        ("OUT_NPZ = REPO / \"talk\" / \"assets\" / \"figures\" / \"_data\" / \"multiz_components.npz\"\n",
         "OUT_NPZ = Path(__file__).resolve().parent / \"outputs\" / \"multiz_components_20planes.npz\"\n"),
        ("SCRATCH = Path(\"/tmp/multiz_components_scratch\")\n",
         "SCRATCH = Path(__file__).resolve().parent / \"outputs\" / \"_multiz_scratch\"\n"),
        ("\"\"\"Compute multi-redshift xi_kappa(gamma) decomposition for the conference talk gif.\n",
         "\"\"\"Multi-redshift xi_kappa(gamma) sweep helper: z -> lambda mapping and the per-plane\nsft-wick sweep used by compute_multiz_kappa_2pcf.py.\n\nMoved here 2026-09-02 from talk/scripts/run_multiz_components.py (written 2026-06-10 for the\nconference-talk animation); the talk keeps a thin wrapper importing this module. Original docstring:\n\nCompute multi-redshift xi_kappa(gamma) decomposition for the conference talk gif.\n"),
    ],
}

GROUP_C_GLOBAL = {
    MCF: [
        ("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/code/mc_fk_complete",
         "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"),
        ("driver_field_emulators/products/table_permclosed_cut15360_permfix_xi.npz", "SFT-lensing-paper-analyses/sachs_sft/" + FKRUN),
        ("driver_field_emulators/products/table_permclosed_cut15360.npz", "SFT-lensing-paper-analyses/sachs_sft/" + PERM + "/table_permclosed.npz"),
        ("products/table_permclosed_cut15360_permfix_xi.npz", FKRUN),
        ("products/table_permclosed_cut15360.npz", PERM + "/table_permclosed.npz"),
        ("driver_field_emulators/code/mc_fk_complete", "SFT-lensing-paper-analyses/sachs_sft/analyses/mc_fk_complete"),
        ("driver_field_emulators/code/callable_fixed", "SFT-lensing-paper-analyses/sachs_sft/" + PERM),
        ("driver_field_emulators/code/rebuild", "SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild"),
        ("driver_field_emulators/code/figures_corrected", "SFT-lensing-paper-analyses/sachs_sft/analyses/revision_2026-08"),
        ("code/figures_corrected/regenerate_mc.py", "analyses/revision_2026-08/regenerate_mc.py"),
    ],
    REV: [
        ("_HERE.parents[2]", "_HERE.parents[4]"),
        ("Path(__file__).resolve().parents[3]", "Path(__file__).resolve().parents[4]"),
        ("../../products/table_permclosed_cut15360.npz", "../../" + PERM + "/table_permclosed.npz"),
        ("../callable_fixed/perm_aware_kappa3_callable.py", "../../" + PERM + "/perm_aware_kappa3_callable.py"),
        ("driver_field_emulators/figures/corrected/", "outputs/"),
        ("driver_field_emulators/products/multiz_ff_real.npz", "sftwick_outputs/2PCF/multiz/multiz_ff_real.npz"),
        ("code/rebuild/fk_kernel_crosscheck.py", "callables/kappa3_vertex/rebuild/fk_kernel_crosscheck.py"),
        ("code/mc_fk_complete/_markers_pooled.npz", "analyses/mc_fk_complete/_markers_pooled.npz"),
    ],
    AUD: [
        ("HERE.parents[1] / \"SFT-lensing-paper-analyses\" / \"sachs_sft\"", "HERE.parents[1]"),
        ("driver_field_emulators/code", "sachs_sft/analyses/fk_audit_2026-08"),
        ("driver_field_emulators/figures", "figures"),
        ("driver_field_emulators/products", "products"),
    ],
    f"{A3}": [
        ("driver_field_emulators/code/figures_corrected/run_ff_one_lambda.py", "analyses/revision_2026-08/run_ff_one_lambda.py"),
        ("code/rebuild/mean_kappa_z.py", "callables/kappa3_vertex/rebuild/mean_kappa_z.py"),
    ],
    f"{AN}/mc_sachs_2pt": [
        ("driver_field_emulators/code/mc_fk_complete", "analyses/mc_fk_complete"),
        ("driver_field_emulators/code/figures_corrected/make_val_figure.py", "analyses/revision_2026-08/make_val_figure.py"),
        ("driver_field_emulators/code/rebuild/fk_kernel_crosscheck.py", "callables/kappa3_vertex/rebuild/fk_kernel_crosscheck.py"),
    ],
}


def apply(entries: dict, label: str) -> None:
    for rel, pairs in entries.items():
        p = PKG / rel
        s = p.read_text()
        for old, new in pairs:
            if old.startswith("re:"):
                pat = re.compile(old[3:], re.DOTALL)
                if not pat.search(s):
                    sys.exit(f"[{label}] {rel}: regex not found: {old[3:60]!r}")
                s = pat.sub(new, s, count=1)
            else:
                if old not in s:
                    sys.exit(f"[{label}] {rel}: string not found: {old[:60]!r}")
                s = s.replace(old, new)
        p.write_text(s)
        print(f"[{label}] patched {rel}")


def apply_global(entries: dict, suffixes=(".py",), label="global") -> None:
    for rel, pairs in entries.items():
        for p in sorted((PKG / rel).rglob("*")):
            if p.suffix not in suffixes or "__pycache__" in p.parts:
                continue
            s = p.read_text()
            n = 0
            for entry in pairs:
                old, new = entry[0], entry[1]
                skip = entry[2] if len(entry) > 2 else ()
                if p.name in skip:
                    continue
                if old in s:
                    s = s.replace(old, new)
                    n += 1
            if n:
                p.write_text(s)
                print(f"[{label}] {p.relative_to(PKG)}: {n} pattern(s)")


if __name__ == "__main__":
    group = sys.argv[1]
    if group == "a":
        apply(GROUP_A, "a")
        apply_global(GROUP_A_GLOBAL, (".py", ".sh"), "a/global")
        apply_global(GROUP_A_JSON, (".json",), "a/json")
    elif group == "b":
        apply(GROUP_B, "b")
    elif group == "c":
        apply(GROUP_C, "c")
        apply_global(GROUP_C_GLOBAL, (".py", ".md", ".wl", ".sh"), "c/global")
    else:
        sys.exit(f"unknown group {group}")
