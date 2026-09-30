# Figure manifest

This mapping describes the R1 revision selected by
[reproduce/active_products.json](reproduce/active_products.json). The manifest
pins the numerical products; [REPRODUCE.md](REPRODUCE.md) explains environment
setup and validation limits. Run `python reproduce/regen_figure.py --list` for
the exact current commands, interpreter choices and output names.

The IDs below are reproduction IDs, not the figure numbering in the manuscript.
The deployed PDFs live in the separate manuscript project's `../figures/`.
Generator paths are relative to this analysis repository. Product keys refer to
the active manifest; other dependencies are stated explicitly.

| ID | Manuscript PDF | Generator | Main inputs | Environment |
|---|---|---|---|---|
| 1 | `analysis1_O0_vs_pyccl.pdf` | `sachs_sft/analyses/analysis1/r1_aligned/plot_comparison.py` | `order0_comparison`, `pyccl_reference` | PyCCL |
| 2 | `analysis3_NLO_FFFK.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_nlo_decomposition.py` | `order0`, `ff`, `fk` | PyCCL |
| 3 | `analysis3_cl_full_vs_O0.pdf` | `sachs_sft/analyses/analysis3/plot_analysis3_cl_decomposition.py` | `order0`, `ff`, `fk` | PyCCL |
| 4 | `cl_EB_polarization.pdf` | `sachs_sft/analyses/analysis3/plot_cl_EB_polarization.py` | `ff`, `fk`; shared curved-sky transform | PyCCL |
| 5 | `multiz_kappa_xi_cl.pdf` | `sachs_sft/analyses/analysis3/plot_multiz_kappa_2x5.py` | `multiz`, `multiz_ff` | PyCCL |
| 6 | `appendix_mc_workflow.pdf` | `sachs_sft/analyses/mc_sachs_2pt/fig_xi_channels.py --from-cache` | `mc_cache`, `mc_markers`; optical callable dependencies | PyCCL + canoes |
| 7 | `zeta_driving_field_slices_draft.pdf` | `sachs_sft/callables/kappa3_vertex/equal_time_limber/plot_zeta_figures.py` | `zeta_slices` | PyCCL |
| 8, 9 | `corr_operator_slices_draft.pdf`, `response_operator_draft.pdf` | `figures/make_response_corr_operator_figures.py` | Stored correlation-propagator table and response callable | PyCCL + canoes |
| 10 | `screen_basis_spin2_pattern.pdf` | `figures/make_screen_basis_spin2_figure.py` | Procedural geometry | PyCCL |
| 11 | `driving_field_portrait.pdf` | `figures/plot_driving_field_portrait.py --seed 42` | CAMB linear power spectrum evaluated at runtime | PyCCL |
| 12 | `selection_rule_schematic.pdf` | `figures/make_selection_rule_schematic.py --seed 42` | Procedural schematic | PyCCL |
| 13 | `null_geodesic_congruence.pdf` | `figures/make_null_congruence_schematic.py` | Procedural schematic | PyCCL |
| 14 | `cosmology_coord_maps.pdf` | `figures/make_cosmology_coord_maps.py` | Background mapping in the generator | PyCCL |
| 15, 16 | `FeynDiag_2pt_order2.pdf`, `FeynDiag_3pt_order1.pdf` | `figures/make_feyndiag.py` | SFT-Wick action and diagram renderer | sft-wick |
| 17 | `fk_cutoff_convergence.pdf` | `sachs_sft/callables/kappa3_vertex/rebuild/fig_cutoff_paper.py` | `cutoff_ladder`, `order0` | PyCCL |

IDs 8/9 and 15/16 each run one shared generator. Regeneration writes local
outputs and never deploys them to the manuscript. The separate
`reproduce/deploy.py` step records deployed-byte checksums in
`reproduce/deployed_md5.txt`.

The [2026-09-29 comparison](sachs_sft/analyses/r1_sft061/closeout_figure_check.txt)
found identical content and zero raster differences for all 17 PDFs after
regeneration. PDF date metadata can prevent raw-byte equality on another run;
font metrics and external environment changes can affect rendering.

Historical arrays under `sachs_sft/sftwick_outputs/2PCF/` remain as comparison
inputs and dependencies of numerical replay. The active manifest replaces their
role in the affected figure readers; their presence does not mean they are the
current plotted results. Do not remove historical inputs without checking the
replay configurations and recorded provenance.
