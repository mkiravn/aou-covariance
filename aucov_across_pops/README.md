# aucov_across_pops

Phenotypic covariance across relatedness, per target population. One notebook per stage, parametrised by `RUN_ID` (a config in `runs/`).

Order: `hm3_overlap` (writes `hm3_ids`; once per CDR) → `ancestry_round1` → `hm3_panel` → `ancestry_round2` → `covariate_pcs` → `grm` → `phenotypes` → `crossproducts` → `curves`.

`uniform_sample` replaces round 1 for a random sample of everyone (`runs/uniform_10k.json`); `even_sample` for one spread evenly over PC space (`runs/even_10k.json`).

Design: `docs/pipeline_architecture.md`. Local runs on simulated data: `docs/local_walkthrough.md`. On the Workbench: `docs/workbench_guide.md`.
