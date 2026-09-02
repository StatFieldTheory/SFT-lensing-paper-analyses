# 2PCF / cutoff_ladder

The ten folded FK sweeps behind figure 17 (`fk_cutoff_convergence.pdf`) and the convergence numbers of the
cutoff subsection: `table_<model>_cut<N>_r4_xi.npz` for model = tree, bihalofit and N = 960, 1920, 3840,
7680, 15360. Each is the production FK fold (sorting callable, 40-point gamma grid, t_final = 2313.029,
n_gauss = 24, n_jobs = 4) of the vertex table `callables/kappa3_vertex/rebuild/products/table_<model>_cut<N>_r4.npz`,
assembled from the octave-band pieces by `rebuild/assemble.py` (2026-08-26; logs in `rebuild/logs/`).

`configs/` holds the verbatim run configs of 2026-08-26 (absolute paths of that machine, callable =
a materialised copy of the production callable with `TABLE_PATH` rewritten to the table). To regenerate:

    cd sachs_sft/callables/kappa3_vertex/rebuild
    sh sweep_cutoff.sh tree 960 1920 3840 7680 15360        # assemble (seconds) + fold (minutes each)
    sh sweep_cutoff.sh bihalofit 960 1920 3840 7680 15360

which writes `products/table_<model>_cut<N>_r4_xi.npz`; copy the folds here. Tree at 0.5' as a
percentage of Order-0: 0.48, 0.97, 1.53, 2.02, 2.31; BiHalofit: 0.78, 2.56, 8.01, 21.5, 41.1.
