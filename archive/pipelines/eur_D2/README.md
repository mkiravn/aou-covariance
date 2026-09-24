# eur_D2 pipeline

Phenotypic covariance by relatedness class for the **eur_D2** sample set:
European-ancestry individuals selected reference-free by a second round of PCA
on `eur_loose` followed by a centroid filter, then a sharded GRM, a
relatedness screen that separates parent-offspring from full-sib pairs, binned
phenotype cross-products, and relatedness-regression estimators.

Notebooks live together under `eur_D2/notebooks/` in the bucket. **Data stays in
the stage folders** the rest of the project uses; this file is the index that
ties them together. Everything below is relative to

```
B = gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9
```

## At a glance

| | |
|---|---|
| individuals | 221,992 |
| GRM variants | 1,179,292 |
| GRM | 24,640,335,028 entries, 91.8 GiB, 16 shards |
| phenotypes | 36 × 2 transforms × 4 covariate sets = 288 models |
| relatedness screen | 14,004 pairs with a_ij > 0.2, 23,503 individuals |
| classified first-degree pairs | 8,887 (PO / FS) |
| bins | 153 (`make_bins.py --wide`) |

## Steps

| # | notebook | original location | what it does |
|---|---|---|---|
| 1 | `01_d2_ancestry_selection.ipynb` | `B/extra_notebooks/EUR-loose-r2-plot-filter-finalPCA.ipynb` | second PCA on eur_loose, centroid filter, PCA on the kept set |
| 1a | `01a_eur_loose_r2_plots.ipynb` | `B/extra_notebooks/EUR-loose-r2-plot.ipynb` | exploratory plots of the eur_loose round-2 PCA |
| 2 | `02_grm_panel_qc.ipynb` | `B/extra_notebooks/EUR_D2_GRM_QC.ipynb` | variant QC of the GRM panel for the D2 keep list |
| 3 | `03_residualize_phenotypes.ipynb` | `B/extra_notebooks/EUR_D2_resid.ipynb` | pull, filter, residualize, standardize within sex |
| 4 | `04_grm_shards_batch.ipynb` | `~/notebooks/EUR_D2_submitGRM.ipynb` | sharded GRM on Google Batch |
| 5 | `05_relatedness_screen_king.ipynb` | `~/notebooks/EUR_D2_POscan.ipynb` | GRM screen, KING on candidates, PO/FS classification |
| 6 | `06_accumulate.ipynb` | `~/notebooks/EUR_D2_accumulate.ipynb` | classed cross-product accumulation, all phenotypes |
| 7 | `07_merge_plots_meta.ipynb` | `~/notebooks/EUR_D2_plot.ipynb` | merge, per-phenotype plots, estimators, meta-analysis |
| 8 | `08_export_results.ipynb` | `~/notebooks/Untitled.ipynb` | small-count-suppressed export zip |

`notebooks/dev/` holds superseded attempts and diagnostics kept for the record;
see the end of this file.

### 1 — D2 ancestry selection

- **input** `B/01_ancestry_filtering/r2_eur_gbr_ceu_projection/eur_loose.{pgen,pvar,psam}`
  — `eur_loose` is the loose Mahalanobis filter on AoU's premade ancestry PCs
- **output** `B/01_ancestry_filtering/r2_eur_gbr_ceu_projection/`
  - `aou_D2_keep_ids.txt` — the eur_D2 membership
  - `aou_eur_loose_D2_filtered.{eigenvec,eigenval,eigenvec.allele}` — PCA fitted on the kept set
  - `aou_eur_loose_D2_filtered.prune.{in,out}` — variants used for that PCA
- **PCA command** (from the plink log): `plink2 --pfile eur_loose --keep aou_D2_keep_ids.txt
  --extract aou_eur_loose_D2_filtered.prune.in --pca approx 20 allele-wts`
- Centroid-filter threshold and pruning parameters: see the notebook.

> The folder name `r2_eur_gbr_ceu_projection` is historical. eur_D2 is
> reference-free; nothing in this step projects onto GBR/CEU.

### 2 — GRM panel QC

- **input** unified panel restricted to `aou_D2_keep_ids.txt`
- **output** `B/03_grm_shards/eur_D2/eur_D2_GRM_QC.{bed,bim,fam}` plus
  `.afreq`, `_freq.frq`, QC logs, and `eur_D2_GRM_LD_check.prune.{in,out}`
- The same panel is staged, with the precomputed `.frq`, for step 4 at
  `B/03_grm_shards/eur_D2/grm_input/`. **`grm_input/` is the canonical copy**;
  see *Housekeeping*.

### 3 — Phenotype residualization

- **inputs** D2 keep list and PCs from step 1; `02_phenotype/docs/phenotype_list.tsv`
- **outputs** `B/02_phenotype/eur_D2/`
  - `raw_pheno_cache/` — per-phenotype pulls
  - `modeling_tables/` — joined phenotype + covariates
  - `residualized/` — `<phenotype>__<raw|invnorm>__<covset>.pheno`, 288 files
- Covariate sets are a nested staircase: `base` → `base_pcs` → `base_pcs_zip3` → `base_pcs_zip3_ses`.
- `.pheno` files carry `FID = IID = person_id`; step 6 keys on IID so the GRM's
  `FID = 0` does not matter.

### 4 — Sharded GRM (Google Batch)

- **input** `B/03_grm_shards/eur_D2/grm_input/` and `B/03_grm_shards/eur_D2/bin/plink`
- **output** `B/03_grm_shards/eur_D2/shards/grm_shard_{k}_of_16.grm.bin.{k}` (k = 1..16),
  `grm_shard_1_of_16.grm.id`
- plink 1.9 `--make-grm-bin --parallel k 16 --read-freq`; `.grm.N.bin` dropped per shard
- worker `n1-custom-16-138240-ext`, `--disk-size 150`, `--use-private-address`
- **dsub must be ≥ 0.5.3 with `CLOUD_SDK_IMAGE` patched to a published tag**
  (`581.0.0-slim` worked). Every released dsub pins a purged cloud-sdk tag, and
  0.5.2's wrappers call `gsutil`, which current images no longer ship. Patch
  and submit in the same shell; any `pip install dsub` reverts the patch.
- `.grm.id` has `FID = 0` for every row.
- Logs: `B/03_grm_shards/eur_D2/dsub_logs/`
- Every shard's size is exactly computable; verify before use (a truncated
  copy is the one failure that yields plausible wrong numbers).

### 5 — Relatedness screen and PO/FS classification

- **inputs** step 4 shards; `grm_input/eur_D2_GRM_QC` panel
- **outputs** `B/03_grm_shards/eur_D2/relatedness/`
  - `candidate_pairs_gt0.2.tsv`, `candidate_individuals_gt0.2.keep` — GRM screen
  - `king_subset.kin0`, `king_panel.log` — kinship and IBS0 on the candidates
  - `deg1_classified.tsv`, `po_pairs_exclude.tsv`, `fs_pairs_confirmed.tsv`
  - `screen_provenance.txt` — every parameter below, as run
- Screen: every pair with a_ij > 0.2 across all 16 shards; float count
  reconciled against N(N+1)/2.
- KING on the 23,503 candidates only: `plink2 --keep … --thin → 50,032 variants --seed 1`,
  then `--make-king-table --king-table-filter 0.15`.
- First degree: `KINSHIP ∈ [0.177, 0.354)`. PO vs FS split on IBS0, **which
  plink2 reports as a proportion, not a count**. The cut is in `screen_provenance.txt`.
- Validated by GRM relatedness: PO pairs are markedly narrower around 0.5 than FS.

### 6 — Classed accumulation

- **inputs** step 4 shards (copied locally and size-verified), step 3
  `residualized/`, step 5 `deg1_classified.tsv`, bins from `make_bins.py --wide`
- **tool** `04_process_shards/src/grm_class_tool` (branch `feat/grm-class-tool`);
  PO and FS pairs are routed into their own bin sets rather than excluded
- `accumulate --pheno-list` per phenotype, all 8 models in one pass per shard;
  576 invocations; `--nblocks 50 --seed 1`
- **output** accumulators on the VM at `~/grm_pheno_cov_eur_D2/*_shard{k}.acc.tsv`
  (not in the bucket; regenerable)
- Filenames carry `__deg1__wide`. Bins are stored as indices, so accumulators
  must be merged with the same bins file.

### 7 — Merge, plots, estimators, meta-analysis

- **outputs** `B/03_grm_shards/eur_D2/crossproducts/`
  - `merged/` — `*_merged.full.tsv`, `*_merged.jk.tsv` per model, including a `pooled` class
  - `summaries/` — per-model summary table
  - `plots/` — `<phenotype>_models.png`
  - `meta/estimators_deg1_wide.tsv`, `meta/plots/`
- Estimators: `04_process_shards/src/he_estimators.py` — h2_Unrel, h2_FS,
  h2_PedW25, h2_Pedf, b2_FS, b2_step, FS/PO excess; delete-block jackknife SEs;
  each on noPO / pooled / PO-only pair sets. Tested against synthetic data by
  `make test`.

### 8 — Export

Allow-listed plots and summary tables only, with every count of 1–20
suppressed, zipped on the VM for download. Bin-level merged files, relatedness
files, accumulators and `.pheno` files are participant-level or carry small
counts and **must not leave the workbench**.

## Code

- `04_process_shards/src/` — `grm_class_tool`, `make_bins.py`, `he_estimators.py`, tests
- Cell sets the notebooks were built from:
  `03_grm_shards/notebooks/remote/eur_D2_grm_batch_cells.md`,
  `04_process_shards/notebooks/remote/00_relatedness_screen_cells.md`,
  `eur_D2_accumulate_cells.md`, `eur_D2_monitor_merge_cells.md`,
  `eur_D2_meta_analysis_cells.md`

## Housekeeping

| item | status |
|---|---|
| `B/03_grm_shards/eur_D2/eur_D2_GRM_QC.bed` (65.4 GB) | same size as `grm_input/` copy; delete once hashes match |
| `B/03_grm_shards/eur_D2/crossproducts/plots/*__deg1__wide.png` | early per-model plots with per-bin pair-count histograms; fine in the bucket, never export |
| `B/01_ancestry_filtering/king_po_exclusion/aou_D2_king.*`, `aou_D2allsegs.txt` | earlier full-cohort KING / `--ibdseg` attempt on eur_D2, superseded by step 5 |
| `B/extra_notebooks/` | originals of steps 1–3; remove after the copies in `eur_D2/notebooks/` are verified |
| stage placement | eur_D2 cross-products live under `03_grm_shards/eur_D2/`, not `04_process_shards/eur_D2/` as for other sample sets |

## notebooks/dev/

| notebook | original | contents |
|---|---|---|
| `grm_local_attempt.ipynb` | `~/notebooks/EUR_D2_GRM.ipynb` | local-VM GRM run, superseded by Batch |
| `phenotype_coverage_checks.ipynb` | `~/scratch_grm/relatedness_screen/Untitled.ipynb` | residualized files vs GRM overlap |
| `classed_accumulate_single_shard.ipynb` | Trash `Untitled 1.ipynb` | first real-shard test of grm_class_tool |
| `accumulate_preflight.ipynb` | Trash `Untitled1.ipynb` | input checks before accumulation |
| `classed_accumulate_concurrent_tests.ipynb` | Trash `Untitled2.ipynb` | 16-shard concurrent test, bins, plots |
| `meta_analysis_first_pass.ipynb` | Trash `Untitled5.ipynb` | first estimator run |

Notebooks contain outputs computed on AoU data. They belong in the workspace
bucket; strip outputs before any of them goes into git.
