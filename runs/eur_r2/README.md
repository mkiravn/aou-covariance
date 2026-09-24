# eur_r2

European-ancestry sample set anchored on 1000G. Round 1 is the 1000G European
Mahalanobis gate from the ancestry survey; round 2 refits PCA on that set and
keeps those closest to the projected CEU + GBR centroid, the anchor Kemper et al.
used. Everything downstream is rerun from scratch on the final set.

Target size: eur_D2's.

```
R = gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9/eur_r2
```

## Steps

| step | cell set | status | writes |
|---|---|---|---|
| round 1 | ancestry survey (prerequisite) | external | `R/01_ancestry/round1/` |
| gate | `gate.ipynb` | ready | `R/01_ancestry/round2/` |
| covariate PCs | `covariate_pcs.ipynb` | ready | `R/01_ancestry/covariate_pca/` |
| GRM panel QC | to port from `EUR_D2_GRM_QC.ipynb` | to write | `R/03_grm/grm_input/` |
| residualize | to port from `EUR_D2_resid.ipynb` | to write | `R/02_phenotype/` |
| GRM shards | to port from `eur_D2_grm_batch_cells.md` | to write | `R/03_grm/shards/` |
| relatedness screen | to port from `00_relatedness_screen_cells.md` | to write | `R/03_grm/relatedness/` |
| accumulate | to port from `eur_D2_accumulate_cells.md` | to write | `R/04_crossproducts/` |
| merge + estimators | to port from `eur_D2_meta_analysis_cells.md` | to write | `R/04_crossproducts/` |
| model ladder | `analyses/model_ladder/ladder_cells.md`, `SAMPLE_SET="eur_r2"` | ready | `B/analyses/model_ladder/eur_r2/` |

## Design

**The gate** is a multivariate normal fitted to the projected CEU + GBR samples:
Mahalanobis distance under it, radius set by a coverage quantile of the anchor's
own distances. The anchor fixes the centre and the relative PC weighting to a
population defined independently of this cohort. Projection shrinks the anchor's
scores toward the origin, but that only rescales distances, which the quantile
absorbs.

Two comparison gates are built alongside it, to bound how much the anchor choice
matters: the participants' own centroid (what eur_D2 did), and the Kemper
direction (PCA fit on the 1000G Europeans, participants projected in).

**Pruning differs by purpose.** The gate prunes at r²=0.05 — it wants clean axes.
The covariate PCA prunes at r²=0.1 — it wants coverage, and at 0.05 the HM3 common
set loses 77% of its variants.

## Not carried over from `analyses/eur_pipeline/`

Both remain in git history at `2d5c9c0`:

- `01b_variant_set_exploration_cells.md` — settled the two variant-set decisions
  above; those conclusions are recorded in `covariate_pcs.ipynb`, so the exploration
  itself is evidence, not pipeline.
- `01c_rare_variant_pca_cells.md` — a side question about rare-variant structure,
  not a step anything downstream depends on.

The per-notebook LD-peak, batch-effect and loadings checks were duplicated 3–4
times; they are now `aoucov.checks` / `aoucov.plots`, called where they matter.
