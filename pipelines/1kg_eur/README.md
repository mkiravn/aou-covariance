# 1kg_eur pipeline

European-ancestry sample set anchored on external references. Round 1 keeps
AoU participants near the centroid of the five **1000G European** populations
in AoU's ancestry PC space. Round 2 refits PCA on that set and keeps those
closest to the projected **CEU + GBR** centroid, the anchor Kemper et al. used.
Everything downstream is rerun from scratch on the final set.

Target for the final set: eur_D2's size.

```
B = gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9
R = B/1kg_eur
```

## Steps

| # | cell set | status | writes |
|---|---|---|---|
| 1 | `01_round1_eur_gate_cells.md` | ready | `R/01_ancestry/round1/` |
| 2 | `02_round2_ceu_gbr_gate_cells.md` | ready | `R/01_ancestry/round2/` |
| 2a | covariate PCs, from `04_final_pca.ipynb` on the round-2 keep list | to write | `R/01_ancestry/final_pca/` |
| 3 | GRM panel QC | to write, from `EUR_D2_GRM_QC.ipynb` | `R/03_grm/grm_input/` |
| 4 | residualize phenotypes | to write, from `EUR_D2_resid.ipynb` | `R/02_phenotype/` |
| 5 | sharded GRM on Batch | to write, from `eur_D2_grm_batch_cells.md` | `R/03_grm/shards/` |
| 6 | relatedness screen, PO/FS | to write, from `00_relatedness_screen_cells.md` | `R/03_grm/relatedness/` |
| 7 | classed accumulation | to write, from `eur_D2_accumulate_cells.md` | VM, then `R/04_crossproducts/` |
| 8 | merge, estimators | to write, from `eur_D2_monitor_merge_cells.md`, `eur_D2_meta_analysis_cells.md` | `R/04_crossproducts/` |
| 9 | model ladder | `analyses/model_ladder/ladder_cells.md` with `SAMPLE_SET = "1kg_eur"` | `B/analyses/model_ladder/1kg_eur/` |

Steps 3 and 4 exist for eur_D2 only as notebooks in the bucket; each needs its
eur_D2 original to adapt.

## Step 1 — round 1

- **inputs** `ancestry_preds.tsv` (AoU participants, 16 PCs) and
  `training_pca.tsv` (HGDP + 1000G reference samples in the same space)
- **gate** Mahalanobis distance over PCs 1–5 under a multivariate normal fitted
  to the 1000G European founders (CEU, GBR, FIN, TSI, IBS): their mean and
  covariance, so the gate follows the cloud's shape; width 1× the distance
  holding 99% of those founders
- same gate as the ancestry survey's EUR group, written here for the pipeline
- **outputs** `1kg_eur_round1_keep_ids.txt` (one `person_id` per line),
  `round1_provenance.txt`, plots
- AoU labels, `eur_loose` and eur_D2 are compared for information only; they
  play no part in selection.

## Step 2 — round 2

- **gate** Mahalanobis distance under a multivariate normal fitted to the
  projected CEU + GBR samples, radius tuned to the target size
- two comparison gates: the participants' own centroid (as eur_D2 did) and the
  Kemper direction (PCA fit on the 1000G Europeans, participants projected)
