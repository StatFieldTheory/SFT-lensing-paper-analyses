#!/bin/zsh
# Phase 2, group (e) leftovers and group (f) default paths:
#  (e) drop the 12 materialised variant callables still tracked under the imported
#      driver_field_emulators copy and move its untracked remainder to the Trash;
#  (f) every paper generator defaults to the corrected inputs and writes into its own
#      outputs/ folder; the three June run configs become folder-relative; the variant
#      runner absolutises inherited relative paths; the FK Monte-Carlo markers become a
#      generator option with the pooled markers as default.
set -e
R=/Users/zzhang/Documents/MyDrafts/STF_lensing
PKG=$R/SFT-lensing-paper-analyses
T=$HOME/.Trash/STF_lensing_reorg_2026-09-02
PY=/opt/homebrew/Caskroom/miniconda/base/envs/PyCCL/bin/python
CO="Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
cd "$PKG"

echo "=== (e) imported-copy leftovers ==="
if [ -d driver_field_emulators ]; then
  git rm -r -q driver_field_emulators
  git commit -q -m "chore: drop the materialised variant callables left under the imported driver_field_emulators copy; the promoted configs and tables carry the same information" -m "$CO"
  echo "committed: $(git rev-parse --short HEAD)"
  mkdir -p "$T/SFT-lensing-paper-analyses"
  mv driver_field_emulators "$T/SFT-lensing-paper-analyses/driver_field_emulators_imported_copy_remainder"
  echo "untracked remainder moved to the Trash"
fi

echo "=== (f1) generator defaults ==="
$PY - <<'PYEOF'
from pathlib import Path
A3 = Path("sachs_sft/analyses/analysis3"); MC = Path("sachs_sft/analyses/mc_sachs_2pt")
ETL = Path("sachs_sft/callables/kappa3_vertex/equal_time_limber"); FG = Path("figures")
FKLINE = 'FK_NPZ = RUNS / "C_corr_op_K_limber_FK" / "xi_C_corr_op_K_limber_FK.npz"\n'
FKNEW = ('# The manuscript\'s FK: ell_max = 15360 with the permutation-aware vertex (2026-08-26).\n'
         '# The June cut1000 sweep in C_corr_op_K_limber_FK/ is superseded and about 4.5x low at 0.5\'.\n'
         'FK_NPZ = RUNS / "C_corr_op_K_limber_FK_cut15360_permfix" / "xi_C_corr_op_K_limber_FK_cut15360_permfix.npz"\n')
edits = {
    A3 / "plot_analysis3_nlo_decomposition.py": [(FKLINE, FKNEW)],
    A3 / "plot_analysis3_cl_decomposition.py": [(FKLINE, FKNEW)],
    MC / "fig_xi_channels.py": [
        ('    _, fk3 = _load_kk("C_corr_op_K_limber_FK")\n', '    _, fk3 = _load_kk("C_corr_op_K_limber_FK_cut15360_permfix")\n'),
        ('def _compute(ff_real: int, ff_nl: int) -> dict:\n', 'def _compute(ff_real: int, ff_nl: int, fk_markers: Path | None = None) -> dict:\n'),
        ('             ffc_mc=ffc_mc, ffc_se=ffc_se)\n    np.savez(_cache_path(), **D)\n',
         '             ffc_mc=ffc_mc, ffc_se=ffc_se)\n'
         '    if fk_markers is not None and Path(fk_markers).exists():\n'
         '        # Placement-complete FK Monte-Carlo markers (analyses/mc_fk_complete/).\n'
         '        m = np.load(fk_markers)\n'
         '        D.update(g_fk=np.asarray(m["gamma"], float), fk_mc=np.asarray(m["fk"], float),\n'
         '                 fk_se=np.asarray(m["err"], float))\n'
         '        print(f"[fk markers <- {fk_markers}]")\n'
         '    np.savez(_cache_path(), **D)\n'),
        ('def main(ff_real: int, ff_nl: int, from_cache: bool) -> None:\n', 'def main(ff_real: int, ff_nl: int, from_cache: bool, fk_markers=None) -> None:\n'),
        ('        D = _compute(ff_real, ff_nl)\n', '        D = _compute(ff_real, ff_nl, fk_markers)\n'),
        ('    a = ap.parse_args()\n    main(a.ff_real, a.ff_nl, a.from_cache)\n',
         '    ap.add_argument("--fk-markers", type=Path, default=FK_MARKERS_DEFAULT,\n'
         '                    help="pooled FK Monte-Carlo markers npz (gamma, fk, err); "\n'
         '                         "default: analyses/mc_fk_complete/_markers_pooled.npz")\n'
         '    a = ap.parse_args()\n    main(a.ff_real, a.ff_nl, a.from_cache, a.fk_markers)\n'),
        ('GAMMA_MC = (1.0, 2.6, 6.7, 17.3, 44.4, 114.3, 293.9,\n',
         'FK_MARKERS_DEFAULT = Path(__file__).resolve().parents[1] / "mc_fk_complete" / "_markers_pooled.npz"\n'
         'GAMMA_MC = (1.0, 2.6, 6.7, 17.3, 44.4, 114.3, 293.9,\n'),
        ('the caller supplies g_fk/fk_mc/fk_se; `_compute` does not produce them, so the\npaper figure is built through\n',
         'the caller supplies g_fk/fk_mc/fk_se; since 2026-09-02 `_compute` takes\n`--fk-markers` (default: the pooled markers) so a fresh run reproduces the deployed\nfigure; the August deployment went through\n'),
    ],
    MC / "driver_stats.py": [
        ('_KAPPA3_DIR = _SACHS_SFT / "callables" / "kappa3_vertex" / "equal_time_limber"\n',
         '_KAPPA3_DIR = _SACHS_SFT / "callables" / "kappa3_vertex" / "equal_time_limber_cut15360_permaware"\n'),
        ('import equal_time_limber_kappa3_callable as _k3  # noqa: E402  coupling_fn -> (3,3,3)\n',
         '# Permutation-aware vertex at ell_max = 15360 since 2026-09-02 (was the cut1000\n'
         '# equal_time_limber_kappa3_callable); only a fresh Monte-Carlo run is affected.\n'
         'import perm_aware_kappa3_callable as _k3  # noqa: E402  coupling_fn -> (3,3,3)\n'),
    ],
    ETL / "plot_zeta_figures.py": [
        ('_NPZ = _HERE / "ell_band_decomp_results.npz"\n',
         '# Converged cutoff (eight HIGH windows to ell = 15360, clipped at 85 arcmin), 2026-08-27;\n'
         '# the June cut1000 ell_band_decomp_results.npz next to this script is superseded.\n'
         '_NPZ = _HERE / "outputs" / "zeta_bands_cut15360_gmax85.npz"\n'),
        ('_FIGDIR = Path("/Users/zzhang/Documents/MyDrafts/STF_lensing/figures")\n',
         '_FIGDIR = _HERE / "outputs"   # the paper copy is deployed by reproduce/deploy.py\n_FIGDIR.mkdir(parents=True, exist_ok=True)\n'),
    ],
    ETL / "ell_band_decomp.py": [
        ('_HIGH_WINDOWS = [(60, 125), (125, 250), (250, 500), (500, 1000)]\n',
         '# Eight disjoint Born-Limber windows to the converged cutoff (2026-08-27); the first\n'
         '# four alone reproduce the June cut1000 decomposition.\n'
         '_HIGH_WINDOWS = [(60, 125), (125, 250), (250, 500), (500, 1000),\n'
         '                 (1000, 2000), (2000, 4000), (4000, 8000), (8000, 15360)]\n'),
        ('    ap.add_argument("--out", type=Path, default=_HERE / "ell_band_decomp_results.npz")\n',
         '    ap.add_argument("--out", type=Path, default=_HERE / "outputs" / "zeta_bands_cut15360.npz")\n'),
    ],
    FG / "make_screen_basis_spin2_figure.py": [('FIG_DIR = ROOT / "figures"\n', 'FIG_DIR = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py\n')],
    FG / "make_response_corr_operator_figures.py": [('FIG_DIR = ROOT / "figures"\n', 'FIG_DIR = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py\n')],
    FG / "make_null_congruence_schematic.py": [('FIG_DIR = ROOT / "figures"\n', 'FIG_DIR = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py\n')],
    FG / "make_selection_rule_schematic.py": [('FIGURES = ROOT / "figures"\n', 'FIGURES = Path(__file__).resolve().parent / "outputs"   # deployed by reproduce/deploy.py\n')],
    FG / "make_cumulant_hierarchy.py": [('FIGURES = ROOT / "figures"\n', 'FIGURES = Path(__file__).resolve().parent / "outputs"   # orphan schematic, not in the paper\n')],
    FG / "plot_driving_field_portrait.py": [('FIGURES = ROOT / "figures"\n', 'FIGURES = HERE / "outputs"   # deployed by reproduce/deploy.py\n')],
}
for p, pairs in edits.items():
    s = p.read_text()
    for old, new in pairs:
        assert old in s, (str(p), old[:60])
        s = s.replace(old, new)
    p.write_text(s); print("patched", p)
PYEOF

echo "=== (f2) run configs relative to their folder; variant runner absolutises inherited paths ==="
$PY - <<'PYEOF'
from pathlib import Path
RUNS = Path("sachs_sft/sftwick_outputs/2PCF")
ABS = "/Users/zzhang/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/sachs_sft/"
for run in ("C_corr_op_O0", "C_corr_op_K_limber_FF", "C_corr_op_K_limber_FK"):
    p = RUNS / run / "config_L2.yaml"; s = p.read_text()
    if ABS not in s:
        print("already relative:", p); continue
    s = s.replace(ABS + f"sftwick_outputs/2PCF/{run}/", "").replace(ABS, "../../../")
    header = ("# Paths are relative to this folder since 2026-09-02: sft-wick resolves module and data\n"
              "# paths against the YAML directory, and scripts/run_FF_single.py changes into it before\n"
              "# running, so output and cache paths resolve here too. Values are otherwise unchanged\n"
              "# from the 2026-06 production run (absolute paths of that machine before this date).\n")
    p.write_text(header + s); print("relative:", p)
rv = Path("sachs_sft/callables/kappa3_vertex/rebuild/run_fk_variant.py"); s = rv.read_text()
old = "    config = yaml.safe_load(base.read_text())\n"
new = old + (
    "    # The run folders write their configs with paths relative to the YAML folder\n"
    "    # (since 2026-09-02); the variant config lives elsewhere, so every inherited\n"
    "    # module and data path is made absolute here before anything is rewritten.\n"
    "    def _abs(p):\n"
    "        return p if Path(p).is_absolute() else str((base.parent / p).resolve())\n"
    "    sysd = config[\"system\"]\n"
    "    sysd[\"linear\"][\"R_time_module\"] = _abs(sysd[\"linear\"][\"R_time_module\"])\n"
    "    sysd[\"noise\"][\"kappa2\"][\"module\"] = _abs(sysd[\"noise\"][\"kappa2\"][\"module\"])\n"
    "    for v in sysd.get(\"vertices\", []):\n"
    "        if \"coupling_path\" in v:\n"
    "            v[\"coupling_path\"] = _abs(v[\"coupling_path\"])\n"
    "    for v in sysd.get(\"nonlocal_vertices\", []):\n"
    "        if \"coupling_module\" in v:\n"
    "            v[\"coupling_module\"] = _abs(v[\"coupling_module\"])\n"
    "    if config.get(\"propagators\", {}).get(\"c_closed_form_module\"):\n"
    "        config[\"propagators\"][\"c_closed_form_module\"] = _abs(config[\"propagators\"][\"c_closed_form_module\"])\n")
assert old in s and "_abs(" not in s
rv.write_text(s.replace(old, new)); print("patched run_fk_variant.materialise_variant")
PYEOF

echo "=== (f3) reproduce/ tools ==="
mkdir -p reproduce
cp reorg_2026-09/phase0/harness/pdfcmp.py reproduce/pdfcmp.py
git add sachs_sft figures reproduce
git commit -q -m "refactor: generators default to the corrected inputs and write into their own outputs folders; run configs relative to their folder; variant runner absolutises inherited paths; FK Monte-Carlo markers become a generator option" -m "$CO"
echo "committed group f: $(git rev-parse --short HEAD)"
echo "=== status ==="; git status --short | head -5
