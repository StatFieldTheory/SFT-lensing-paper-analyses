# Historical audit context

This note preserves pre-acceptance audit context. Pending/candidate wording below
describes that earlier stage. The current public selection and finite acceptance
boundaries are recorded in active_products.json and public_product_lineage.json.

# Historical FK documentation errata, current-facing addendum

The unchanged June sweep README at sachs_sft/sftwick_outputs/2PCF/C_corr_op_K_limber_FK/README.md:55-58 calls the August +1.9480e-5 value the manuscript's converged value and points to the August permfix folder. Those statements describe historical products. They do not select current inputs or establish convergence of a rebuilt candidate. Preserve that historical README, its numbers and archived products. Current consumers use reproduce/active_products.json with its explicit hashes; candidate activation is a separate accepted-bundle action. Do not replace the August number with a running-fold amplitude.

The September vertex rebuild note is located at reproduce/provenance/vertex_rebuild_20260906.md, with a root HANDOFF_2026-09-06_nphi512_vertex_rebuild.md historical copy. Neither should be rewritten as T-001 acceptance. This resolves the earlier missing docs/design locator.

BRIEF C5's 33-row transcript and other-30-rows-unchanged wording describe the historical flip-control procedure. FIX_LIST's later corrections specify a current 32-row checker. Static source inspection finds 28 add call sites, including loop-expanded and conditional entries; call-site count is not output-row count. The 32-row statement is from FIX_LIST, not a newly measured checker output. A blanket 30-unchanged criterion must not be imposed on rebuilt FK/FF products: N2/N3/N5/N7/N8 can change entries, and the author significance rule governs manuscript changes separately from corrected package products. For B1 use only its dedicated C3 invariants and parent_b1_c3_control_array_review_20261004_v2.json; do not infer an unchanged checker roster by subtracting two from 32. For a future checker comparison, match explicit claim labels and classify changed, unchanged, removed or newly added claims. Archive existing JSON/stdout before running the checker, which writes reproduce/numbers.json.

Odd-spin scientific closure against an accepted candidate remains pending. Operational sign-seam tests and the saved B1 fold are controls, not an external check of accepted-candidate TTP/PPP. Retain DECISION option A and local E as primitive; do not substitute E=E0/a or identify FK with reduced shear.
