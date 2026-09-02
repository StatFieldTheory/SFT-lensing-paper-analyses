# analysis3/inputs

* `multiz_components_talk_2026-06-10.npz`: byte-identical copy of
  `talk/assets/figures/_data/multiz_components.npz` (talk repository, 2026-06-10 06:15), the 20-plane
  sweep (z_s = 0.5 to 5, keys z, gamma, o0, ff, fk, lam) that `run_multiz_components.py` produced for
  the conference-talk animation with the production configs. The Order-0 planes of
  `outputs/multiz_kappa_2pcf_5z.npz`, and hence of the deployed figure 5, are PCHIP-in-z interpolants
  of its `o0` array at z_s = 1, 1.7, 2.5, 3.2, 4 (maximum relative difference 0.0, checked 2026-09-02).
  Its `ff` planes are the talk-era placeholder (z = 5 FF scaled per gamma) and its `fk` planes are the
  corner-frozen cut1000 vertex; neither is used by the paper any more.
