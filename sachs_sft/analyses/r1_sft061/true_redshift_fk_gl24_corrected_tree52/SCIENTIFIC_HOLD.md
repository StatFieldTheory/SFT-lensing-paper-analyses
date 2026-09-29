# Scientific hold before any fold

The prepared restricted ladder uses the historical geometry-sorting callable.
A source review identified that this callable can misassign spin-dependent K
components. The kappa-kappa FK contraction includes Ricci-Weyl modulus terms,
so scalar external indices do not make this sorting harmless. The restriction
check proved equality to the historical sorting callable only.

No fold has started in this directory. The original plan is preserved as
`prepared_blocked_sorting_callable.json`. Its missing active filename makes
the existing serial coordinator fail before starting this batch. Use a fresh
batch with a verified permutation-aware contraction for final science.
