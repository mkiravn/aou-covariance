# Analyses

Two kinds of work, kept apart:

- **`pipelines/<sample_set>/`**: one full run per sample set, from ancestry
  selection to merged cross-products. A finished run is frozen; changes go into
  a new sample set or a new analysis, not into its notebooks.
- **`analyses/<name>/`**: work that reads the outputs of one or more runs,
  or AoU's own reference data, and does not define a sample set.

| directory | reads | question |
|---|---|---|
| `pipelines/eur_D2/` | | reference-free European set (done) |
| `pipelines/1kg_eur/` | | European set anchored on 1000G CEU (in progress) |
| `analyses/ancestry_survey/` | AoU ancestry PCs, original 1000G and HGDP labels, AoU relatedness | how large and how cluster-like is each reference group in AoU, and how many first-degree pairs does each have |
| `analyses/model_ladder/` | a run's merged cross-products | does the relatedness slope differ between unrelated and related pairs, before and after degree offsets |

Shared code stays in `04_process_shards/src/` (`grm_class_tool`,
`he_estimators.py`, `make_bins.py`).

## Bucket layout

```
B = gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9

B/<sample_set>/                 new runs keep everything under one root
    01_ancestry/  02_phenotype/  03_grm/  04_crossproducts/  notebooks/
B/analyses/<name>/[<sample_set>/]
```

eur_D2 predates this and keeps its data in the stage folders; see
`pipelines/eur_D2/README.md`.

Every notebook is a markdown cell set; paste the cells into a notebook in the
matching bucket `notebooks/` folder. Clear outputs before anything goes to git.
