# 1kg_eur pipeline

European-ancestry sample set anchored on external references. Round 1 (the
1000G European Mahalanobis gate) is produced by the ancestry survey notebook
and consumed here as a prerequisite. Round 2 refits PCA on that set and keeps
those closest to the projected **CEU + GBR** centroid, the anchor Kemper et al.
used. Everything downstream is rerun from scratch on the final set.

Target for the final set: eur_D2's size.

```
B = gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9
R = B/eur_r2
```

## Steps

| # | cell set | status | writes |
|---|---|---|---|
| 0 | ancestry survey notebook (prerequisite) | external | `R/01_ancestry/round1/` |
| 1 | `runs/eur_r2/gate.ipynb` | ready | `R/01_ancestry/round2/` |
| 2a | `runs/eur_r2/covariate_pcs.ipynb` | ready | `R/01_ancestry/covariate_pca/` |
| 3 | GRM panel QC | to write, from `EUR_D2_GRM_QC.ipynb` | `R/03_grm/grm_input/` |
| 4 | residualize phenotypes | to write, from `EUR_D2_resid.ipynb` | `R/02_phenotype/` |
| 5 | sharded GRM on Batch | to write, from `eur_D2_grm_batch_cells.md` | `R/03_grm/shards/` |
| 6 | relatedness screen, PO/FS | to write, from `00_relatedness_screen_cells.md` | `R/03_grm/relatedness/` |
| 7 | classed accumulation | to write, from `eur_D2_accumulate_cells.md` | VM, then `R/04_crossproducts/` |
| 8 | merge, estimators | to write, from `eur_D2_monitor_merge_cells.md`, `eur_D2_meta_analysis_cells.md` | `R/04_crossproducts/` |
| 9 | model ladder | `analyses/model_ladder/ladder_cells.md` with `SAMPLE_SET = "1kg_eur"` | `B/analyses/model_ladder/1kg_eur/` |

Steps 3 and 4 exist for eur_D2 only as notebooks in the bucket; each needs its
eur_D2 original to adapt.

## Step 1 — round 2

- **gate** Mahalanobis distance under a multivariate normal fitted to the
  projected CEU + GBR samples, radius tuned to the target size
- two comparison gates: the participants' own centroid (as eur_D2 did) and the
  Kemper direction (PCA fit on the 1000G Europeans, participants projected)
