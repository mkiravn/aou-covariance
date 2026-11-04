# aucov_across_pops: architecture

Draft, 2026-11-04. Nothing is built yet. This describes one parametrised pipeline
that runs once per target population (3–4 runs) and shares code and conventions
with the existing CEU/GBR pipeline (`notebooks/01–12`, run id `1kg_eur`).

## What each requirement maps to

| requirement | mechanism |
|---|---|
| reports results and figures | every stage writes figures + an aggregate summary to `{root}/report/{stage}/` through one helper that enforces AoU's no-counts-of-1–20 rule; `python -m aoucov --run X report` assembles `RESULTS.md` |
| code visible and editable | each stage notebook holds all of its analysis code: Python functions at the top, plink as `%%bash` cells, bigsnpr as `%%script Rscript` cells. Only plumbing (`aoucov`) and code shared with CEU/GBR (`residualize_lib.R`, `he_estimators.py`, `grm_class_tool`) live elsewhere |
| aware of input and output | `aoucov` (already drafted, untracked): a committed run config declares each stage's outputs and `needs`; `python -m aoucov --run X status` shows done/missing/blocked; each output gets a `.meta.json` sidecar with git SHA, config hash, machine |
| documented on GitHub | repo holds code, configs, this doc, per-stage READMEs; results reach GitHub only after you review the downloaded `report/` |
| AoU best practice | compute stays in the workbench; Batch via dsub with the known pet-SA/private-network config; CDR-provided resources used where they exist (premade PCs, ACAF callset, relatedness table); CDR version pinned in config; count suppression built into reporting |
| shared with CEU/GBR | one gate module, one LD-region file, one phenotype library, one Batch helper; `1kg_eur` gets a run config pointing at its existing outputs so both runs show in `status` and compare directly |

No workflow manager. Snakemake/Nextflow were rejected earlier for good reasons, and
they still hold: the gates are interactive, and dsub already works.

## Layout

```
aoucov/                 run bookkeeping (exists): config, store, manifest, provenance, runner
  report.py             NEW  figures/tables/summary with count suppression
src/
  grm_class_tool.cpp, he_estimators.py, covplot.py, make_bins.py   (unchanged)
resources/              NEW  committed, public, non-participant files
  high_ld_regions_hg38.txt   was inlined in 02_round2_gate
runs/
  1kg_eur.json          exists, describes the finished CEU/GBR run
  eur.json  afr.json  eas.json  sas.json      NEW, one per target population
aucov_across_pops/      NEW  stage notebooks, named by stage id, parametrised by RUN_ID
scripts/sim/, scripts/dev/   simulator, run-all, local setup
notebooks/              the 1kg_eur record; not edited
```

A stage notebook opens with:

```python
RUN_ID = "afr"                       # papermill parameter cell
cfg = load_run(RUN_ID)
rep = Report(cfg, "ancestry_round2")
```

Run each notebook interactively the first time. Decisions you make by eye (K, radius,
which PCs are LD artefacts) go into the run config, not into a cell. After that,
`papermill aucov_across_pops/X.ipynb - -p RUN_ID eas` re-executes a stage for another
population, and the executed copy, with outputs, lands in the bucket next to the
report. The repo copy stays output-free.

## Run config: what varies per population

```json
{
  "run_id": "afr",
  "root": "gs://…/phenotypic_covariance_v9/afr",
  "params": {
    "cdr": "v9",
    "r1_ref_pops": ["YRI", "LWK", "GWD", "MSL", "ESN"],
    "r1_n_pcs": 5,
    "r1_width": 0.99,
    "r2_train_pops": ["YRI", "LWK", "GWD", "MSL", "ESN"],
    "r2_anchor_pops": ["YRI", "ESN"],
    "r2_k": null,
    "r2_width": null,
    "grm_maf": 0.01,
    "bins": {"fine": 0.001, "mid": 0.005, "coarse": 0.01},
    "nblocks": 50
  },
  "stages": { "...": "same stage list for every run, see below" }
}
```

`null` means not decided yet. A stage that needs it refuses to run.

## Stages

Every stage reads only what `needs` names and writes only under its own `dir`.

### 1. `ancestry_round1`: premade PCs, ellipsoid around a 1000G superpopulation
- **In:** `ancestry_preds.tsv` (16 PCs per participant), `training_pca.tsv` (reference panel in the same space), 1000G/HGDP labels.
- **Method:** `01_round1_gate`'s code, copied into the notebook. MVN fit to the founders of `r1_ref_pops` on the first `r1_n_pcs` PCs; radius = the `r1_width` quantile of the founders' own distances.
- **Out:** keep list. **Report:** kept-N vs coverage curve (a plateau means a cluster, steady growth means an admixture continuum), gate in PC space, kinship-class counts inside the gate.
- Same code for every superpopulation, so `1kg_eur` round 1 and the new runs can't drift apart.

### 2. `ancestry_round2`: PCA on the 1000G superpopulation, AoU projected in
This is where Privé et al. 2020 (doi:10.1093/bioinformatics/btaa520) applies most directly.

- **Variants:** HM3 SNPs, QC'd in the round-1 set (MAF ≥ 1%, HWE, missingness), matched to 1000G on ID+REF+ALT, **A/T and C/G removed**, MAF ≥ 1% in the training subset too.
- **Training set:** 1000G founders of `r2_train_pops`, KING-unrelated (< 0.0884), with LOF outliers removed (`bigutilsr::prob_dist`). That's about 500 people for a superpopulation.
- **PCA:** `bigsnpr::bed_autoSVD`: clumping at r² = 0.2, then iterative removal of long-range-LD outlier variants. Our hardcoded hg38 LD list is excluded as well.
- **Projection:** `bigsnpr::bed_projectPCA` with **OADP**, not `plink2 --score`. With ~500 training samples, simple projection shrinks AoU scores substantially from about PC4 onward. The paper reports shrinkage factors of 1.36 at PC5 and 2.84 at PC10 even with all 2,500 1000G samples. Because the gate compares in-sample 1000G scores against projected AoU scores, shrinkage would put the gate in the wrong place.
- **K:** pick PCs that show population structure in the scores and have no single-region peak in the loadings. Probably 3–5 within a superpopulation. Set `r2_k` in the config.
- **Gate:** robust Mahalanobis distance (OGK in bigutilsr, or MCD in Python) around `r2_anchor_pops` in the first K PCs.
- **Runs on:** one Batch job (R + bigsnpr), about 300K people × ~100–200K variants. The paper projected 488K UKBB in 12 minutes on 20 cores. The gate itself runs in the notebook.
- **Report:** simple vs OADP scores (shows the shrinkage), loadings by position, scores coloured by 1000G population, kept-N vs radius, overlap with `1kg_eur`'s keep list (EUR run only).

### 3. `covariate_pcs`
HM3 common SNPs, r² = 0.1, 20 PCs, fit within the final sample. This is the 02a decision, kept for consistency. See "Relatives" below for which people are used to fit them.

### 4. `grm_panel`: about 10⁶ variants
- HM3 biallelic SNPs from the ACAF callset, QC'd **within the target population**: MAF ≥ 1%, HWE p > 1e-6, missingness ≤ 5%, no LD pruning (standard GCTA practice). In a single ancestry this comes out at about 1.0–1.1M SNPs, close to your target.
- Allele frequencies are computed once and passed to every shard with `--read-freq`, as `04_grm_panel_qc` does.
- **Differs from `1kg_eur`:** its GRM used the unified panel (HM3 ∪ 20% random thin, about 1.69M variants). See open decision 2.

### 5. `grm_shards`
Unchanged from `06_grm_shards`: plink 1.9 `--make-grm-bin --parallel k N` on Batch. The shard format is the **interface** every relatedness measure has to produce (see 9).

### 6. `relatedness`: PO vs FS
`06a_close_pair_screen` → KING on the screened subset → `deg1_classified.tsv`. The code is the same; only the config changes.

### 7. `phenotypes`
Same `residualize_lib.R` and `docs/phenotype_list.tsv`, run on this run's keep list and covariate PCs. "Well behaved" becomes explicit per-run checks, recorded in the report: N after QC ≥ a floor, enough distinct values (excludes near-binary survey items), skew after invnorm, and plausible-range drop rate. A phenotype that fails is listed and excluded, not silently dropped.

### 8. `crossproducts` → `merged` → `curves`
`grm_class_tool accumulate/merge` with `--pheno-list` batching, then `he_estimators.py` and the curve figures. These are already declared in `runs/1kg_eur.json`.

### 9. Later: other relatedness measures, same downstream
To let PC-Relate and IBD reuse 8 unchanged, each produces **pair values on the GRM's scale** (relatedness = 2φ):

- **`pcrelate_shards`:** PC-Relate kinship is a ratio of two matrix products (centred genotypes / individual-specific variances, with μᵢₛ from the covariate PCs). It can be computed in the same row-block partition and written as `.grm.bin` shards, so `grm_class_tool` reads it as is. Hail's `pc_relate` (what AoU used for its own table) has to filter pairs to stay tractable, so it can't give every pair.
- **`ibd`:** iLASH needs **phased** haplotypes. Check whether CDR v9 ships a phased callset. If not, phasing (Beagle/SHAPEIT5 on Batch) is a prerequisite and a significant cost. You declined phasing for the haplotype components earlier, so this needs a separate decision. IBD is sparse (most pairs = 0), so it would need a sparse-input mode in `grm_class_tool`. The zero bin follows by subtraction: Σᵢ<ⱼ yᵢyⱼ = ((Σy)² − Σy²)/2, per jackknife block too. IBD can't resolve the a < 0.02 region, so it informs the relatives' end of the curve only.

## Relatives: nobody is discarded

The paper removes relatives only from the set used to **fit** PCs, so that projected people aren't related to anyone in the fit. That requirement plays out differently at each PCA step:

- **Round 2:** the PCs are fit on 1000G. AoU relatives are only projected, and they're never related to 1000G samples. Every AoU participant keeps a score and faces the same gate. Nothing to do.
- **Covariate PCs:** these are fit in AoU, so there's a choice. (a) Fit on everyone, as `1kg_eur` did. (b) Fit on the unrelated subset (both members of every pair with AoU kinship > 0.1 left out of the fit), then project everyone else with `bed_projectSelfPCA`. Under both options **every participant gets PCs and stays in the analysis**. I'd pick (b): it's the paper's procedure, and it's also what PC-Relate expects (PC-AiR PCs), so step 9 gets its PCs for free. It costs one projection call. It does differ from `1kg_eur`.
- **GRM, phenotypes, cross-products:** relatives stay. They're the point of the analysis.

## Things to fix in passing (not blocking)

- **Coverage quantiles above ~0.99 aren't estimable from the anchor.** `02_round2_gate` used `COVERAGE = 0.999` on CEU+GBR, about 190 samples, which is effectively the anchor's maximum distance. For the new runs, set the radius from χ²_K of the robust fit, or report it as a radius with the coverage curve beside it.

## Reporting and GitHub

- `Report` writes `{root}/report/{stage}/`: PNGs, aggregate TSVs, `summary.json` (N in/out, parameters, chosen K and radius). `rep.table(df, counts=[...])` raises if a count column holds 1–20 unmasked. This makes the existing `show()`/`≤20` convention impossible to forget.
- `python -m aoucov --run X report` stitches the summaries into `{root}/report/RESULTS.md`.
- The VM never pushes (no write credentials). To put results on GitHub, you download a `report/` folder, check it against AoU's dissemination policy, and commit it to `results/<run_id>/`. That check stays a manual step.

## Open decisions

1. **Round-2 anchors** per population: CEU+GBR for EUR as before? YRI(+ESN) for AFR, CHB+JPT for EAS, GIH/PJL for SAS? Most AoU African-ancestry participants sit on an admixture cline, so the AFR round-1 gate may hold few people. The curve in stage 1 will show it.
2. **GRM variant set:** HM3-only (~1.0M) for the new runs, against `1kg_eur`'s ~1.69M. Keep the difference, or rebuild `1kg_eur`'s GRM on HM3 so all runs match?
3. **Covariate PCA fit set:** (a) or (b) above.
4. **IBD:** is there a phased callset, and is phasing in scope if not?

Bins are a run parameter (`params.bins`: fine/mid/coarse widths), generated by `src/make_bins.py` at the cross-product stage.

## Starting point: PLINK, not Hail

The ACAF callset ships as Hail MT, VCF and per-chromosome PLINK. Every HapMap3 SNP is in ACAF, so the pipeline starts from the PLINK files with plink2 on Batch (as `02_build_unified_panel` already does). Hail would need a Spark cluster for no gain here. If a later step needs Hail (rare variants from the VDS), keep it to one stage that ends in `hl.export_plink`, so everything downstream is unchanged.

## Robust fit: MCD, not OGK

The round-2 radius comes from chi-square, which needs a calibrated scatter. `bigutilsr::covrob_ogk` isn't consistency-corrected (variances ≈ 0.61 on N(0,1) data), so the gate uses reweighted MCD (`sklearn.covariance.MinCovDet`). Both are high-breakdown estimators.

## Local runs on simulated data

`docs/local_walkthrough.md`. `runs/sim.json` and `runs/eur.json` share stages and parameters; only `inputs`, `tools`, `backend` and sizes differ.
