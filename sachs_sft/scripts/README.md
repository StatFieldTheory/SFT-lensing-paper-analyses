# scripts/

ONLY genuinely-shared, necessary utility code (cosmology helpers, the sft-wick
run harness, common plotting). NOT analysis logic, NOT callable implementations
(those live in their own folders under `callables/`).

Keep this small. If something is used by only one callable or one analysis, it
belongs in that folder, not here.
