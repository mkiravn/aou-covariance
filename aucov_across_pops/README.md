# aucov_across_pops

Phenotypic covariance across relatedness, per target population. One notebook per stage, parametrised by `RUN_ID` (a config in `runs/`).

Order: `ancestry_round1` → `hm3_panel` → (`hm3_overlap`, optional diagnostic) → `ancestry_round2` → `covariate_pcs` → `grm` → `phenotypes` → `crossproducts` → `curves`.

Design: `docs/pipeline_architecture.md`. Local runs on simulated data: `docs/local_walkthrough.md`.
