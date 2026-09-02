# Phase-0 figure checks (2026-09-02)

Deployed `figures/*.pdf` against candidates regenerated into the session scratch dir by `harness/p0_regen.py` and compared with `harness/pdfcmp.py`.

| # | deployed figure | candidate (how produced) | verdict | raster diff frac (60 dpi) | page size deployed / candidate |
|---|---|---|---|---|---|
| 1 | analysis1_O0_vs_pyccl.pdf | standalone, sft-wick env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 903.887 x 654.524 pts / 903.887 x 654.524 pts |
| 1 | analysis1_O0_vs_pyccl.pdf | standalone, PyCCL env | DIFFERS | SHAPE DIFFERS | 903.887 x 654.524 pts / 903.88 x 653.555 pts |
| 2 | analysis3_NLO_FFFK.pdf | standalone (default FK_NPZ = cut1000 sweep), PyCCL env | DIFFERS | 0.02316 | 903.88 x 653.555 pts / 903.88 x 653.555 pts |
| 2 | analysis3_NLO_FFFK.pdf | FK_NPZ rebound to cut15360 permfix, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 903.88 x 653.555 pts / 903.88 x 653.555 pts |
| 2 | analysis3_NLO_FFFK.pdf | FK_NPZ rebound, sft-wick env | DIFFERS | SHAPE DIFFERS | 903.88 x 653.555 pts / 903.887 x 654.524 pts |
| 3 | analysis3_cl_full_vs_O0.pdf | standalone (cut1000 FK), PyCCL env | DIFFERS | 0.016166 | 897.002 x 833.912 pts / 897.002 x 833.912 pts |
| 3 | analysis3_cl_full_vs_O0.pdf | FK_NPZ rebound to cut15360 permfix, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 897.002 x 833.912 pts / 897.002 x 833.912 pts |
| 4 | cl_EB_polarization.pdf | standalone (cut1000 FK), PyCCL env | DIFFERS | 0.010516 | 1061.64 x 378.328 pts / 1061.64 x 378.328 pts |
| 4 | cl_EB_polarization.pdf | FK_NPZ rebound to cut15360 permfix, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 1061.64 x 378.328 pts / 1061.64 x 378.328 pts |
| 5 | multiz_kappa_xi_cl.pdf | plot stage from cached 5z npz + multiz_ff_real_all5.npz, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 1095.09 x 494.153 pts / 1095.09 x 494.153 pts |
| 6 | appendix_mc_workflow.pdf | fig_xi_channels --from-cache, PyCCL env + canoes on path | IDENTICAL (content; only /CreationDate differs) | 0.0 | 463.44 x 348.24 pts / 463.44 x 348.24 pts |
| 6 | appendix_mc_workflow.pdf | make_val_figure route (FK line cut15360 permfix + mc_fk_complete markers), PyCCL env + canoes | IDENTICAL (content; only /CreationDate differs) | 0.0 | 463.44 x 348.24 pts / 463.44 x 348.24 pts |
| 7 | zeta_driving_field_slices_draft.pdf | standalone (default June ell_band_decomp_results.npz), PyCCL env | DIFFERS | SHAPE DIFFERS | 987.6 x 372.931 pts / 987.6 x 369.336 pts |
| 7 | zeta_driving_field_slices_draft.pdf | _NPZ rebound to zeta_bands_cut15360_gmax85.npz, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 987.6 x 372.931 pts / 987.6 x 372.931 pts |
| 8 | corr_operator_slices_draft.pdf | standalone, PyCCL env + canoes | DIFFERS | 0.007744 | 661.242 x 601.604 pts / 661.445 x 601.815 pts |
| 9 | response_operator_draft.pdf | standalone, PyCCL env + canoes | DIFFERS | 0.02133 | 407.76 x 256.551 pts / 408 x 256.791 pts |
| 10 | screen_basis_spin2_pattern.pdf | standalone, PyCCL env | DIFFERS | 0.017811 | 509.686 x 275.218 pts / 509.966 x 275.666 pts |
| 10 | screen_basis_spin2_pattern.pdf | standalone, sft-wick env | DIFFERS | 0.018015 | 509.686 x 275.218 pts / 509.966 x 275.666 pts |
| 11 | driving_field_portrait.pdf | --seed 42, PyCCL env (camb 1.6.5) | IDENTICAL (content; only /CreationDate differs) | 0.0 | 1022.2 x 316.842 pts / 1022.2 x 316.842 pts |
| 12 | selection_rule_schematic.pdf | --seed 42, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 461.88 x 165.096 pts / 461.88 x 165.096 pts |
| 12 | selection_rule_schematic.pdf | --seed 42, sft-wick env | DIFFERS | 0.001412 | 461.88 x 165.096 pts / 461.88 x 165.096 pts |
| 13 | null_geodesic_congruence.pdf | standalone, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 394.2 x 200.412 pts / 394.2 x 200.412 pts |
| 13 | null_geodesic_congruence.pdf | standalone, sft-wick env | DIFFERS | 0.000887 | 394.2 x 200.412 pts / 394.2 x 200.412 pts |
| 14 | cosmology_coord_maps.pdf | archive generator + scripts.canoes_pipeline shim, PyCCL env + canoes | IDENTICAL (content; only /CreationDate differs) | 0.0 | 837.84 x 333.84 pts / 837.84 x 333.84 pts |
| 14 | cosmology_coord_maps.pdf | same, canoes .venv (matplotlib 3.10.9) | DIFFERS | 0.0 | 837.84 x 333.84 pts / 837.84 x 333.84 pts |
| 15 | FeynDiag_2pt_order2.pdf | sft-wick/examples/FeynDiag_2pt_order2.pdf (not rerun) | IDENTICAL (raw bytes) | 0.0 | 841.18 x 624.309 pts / 841.18 x 624.309 pts |
| 16 | FeynDiag_3pt_order1.pdf | sft-wick/examples/FeynDiag_3pt_order1.pdf (not rerun) | IDENTICAL (raw bytes) | 0.0 | 853.338 x 624.309 pts / 853.338 x 624.309 pts |
| 17 | fk_cutoff_convergence.pdf | fig_cutoff_paper.py --products, PyCCL env | IDENTICAL (content; only /CreationDate differs) | 0.0 | 463.44 x 319.44 pts / 463.44 x 319.44 pts |
| 2* | analysis3_NLO_FFFK.pdf | in-tree twin: dfe/figures/corrected/analysis3_nlo_O0_FF_FK_corrected.pdf | IDENTICAL (raw bytes) | 0.0 | 903.88 x 653.555 pts / 903.88 x 653.555 pts |
| 2* | analysis3_NLO_FFFK.pdf | in-tree STALE twin: analysis3/outputs/analysis3_nlo_O0_FF_FK.pdf (June) | DIFFERS | 0.161427 | 903.88 x 653.555 pts / 903.88 x 653.555 pts |
| 3* | analysis3_cl_full_vs_O0.pdf | in-tree STALE twin: analysis3/outputs/analysis3_cl_O0_FF_FK.pdf (June) | DIFFERS | 0.024676 | 897.002 x 833.912 pts / 897.002 x 833.912 pts |
| 4* | cl_EB_polarization.pdf | in-tree STALE twin: analysis3/outputs/cl_EB_polarization.pdf (June) | DIFFERS | 0.013477 | 1061.64 x 378.328 pts / 1061.64 x 378.328 pts |
| 5* | multiz_kappa_xi_cl.pdf | in-tree twin: analysis3/outputs/multiz_kappa_xi_cl_2x5.pdf | IDENTICAL (raw bytes) | 0.0 | 1095.09 x 494.153 pts / 1095.09 x 494.153 pts |
| 6* | appendix_mc_workflow.pdf | in-tree twin: mc_sachs_2pt/figures/xi_kappa_channels.pdf | IDENTICAL (raw bytes) | 0.0 | 463.44 x 348.24 pts / 463.44 x 348.24 pts |
| 7* | zeta_driving_field_slices_draft.pdf | in-tree twin: dfe/figures/corrected/zeta_driving_field_slices_draft.pdf | IDENTICAL (raw bytes) | 0.0 | 987.6 x 372.931 pts / 987.6 x 372.931 pts |
| 1* | analysis1_O0_vs_pyccl.pdf | in-tree twin: analysis1/outputs/analysis1_O0_vs_pyccl_fkem.pdf | IDENTICAL (raw bytes) | 0.0 | 903.887 x 654.524 pts / 903.887 x 654.524 pts |
