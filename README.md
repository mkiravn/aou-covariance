# aou-covariance

Phenotypic covariance / heritability pipeline on All of Us, following the
Kemper et al. 2021 design (*"Phenotypic covariance across the entire
spectrum of relatedness for 86 billion pairs of individuals"*) adapted to
the AoU Researcher Workbench. Uses [GRM-pairs](https://github.com/mkiravn/GRM-pairs)
(vendored as a submodule) for the sparse/binned GRM tooling and the
individual-block jackknife.

**Status: scaffolding only.** Directory structure reflects the intended
pipeline; nothing inside is a finalized decision until it's been tried
against the real workbench. See each stage's own `README.md`.

## Pipeline stages

| Stage | Purpose |
|---|---|
| `01_ancestry_filtering/` | Restrict to a genetically homogeneous reference cluster via a PCA ellipsoid filter, validated against 1000G. |
| `02_phenotype/` | Pull phenotypes + covariates from AoU, residualize, trim outliers, standardize within sex. |
| `03_grm_shards/` | Compute the GRM in row-chunked shards (`plink --make-grm-bin --parallel k n`) as AoU batch jobs. |
| `04_process_shards/` | Bin each shard by relatedness and accumulate phenotype cross-products (`GRM-pairs/grm_bin_sharded`), merge. |

The final analysis ID list will be the intersection of: passes ancestry
filtering (01), has a usable phenotype after normalization (02), and has
genotype data included in the GRM (03/04).

## Data handling — read before adding anything here

This repo is code only. **Never commit, in any form:**

- Participant-level data of any kind — phenotypes, survey responses,
  covariates, demographics, genotypes.
- Participant identifiers — `person_id`/`research_id`, sample IDs, or any
  ID list derived from an AoU cohort (a list of *public* 1000 Genomes or
  other public reference-panel sample IDs is fine — that's not
  participant data).
- Jupyter notebook **cell outputs** if the notebook was ever run against
  real AoU data — outputs can silently embed real values, counts, or
  identifiers. Clear all outputs before committing
  (`jupyter nbconvert --clear-output --inplace <notebook>.ipynb`), or set
  up [`nbstripout`](https://github.com/kynan/nbstripout) so it happens
  automatically: `pip install nbstripout && nbstripout --install`.
- Any extracted results table, even aggregated, that hasn't gone through
  the AoU publication/egress review — that review, not this repo, is the
  only sanctioned path for anything to leave the workbench.

If you're ever unsure whether something is safe to commit, don't — ask
first. `.gitignore` blocks the obvious file types (genomic formats,
tabular data extensions, notebook checkpoints) as a backstop, not a
substitute for checking before `git add`.

## Layout convention

Each numbered stage gets a `docs/` (parameters, validation results),
`scripts/` (`local/` for anything prototyped outside the workbench, `aou/`
for what actually runs there), and `notebooks/` (`local/` vs `remote/`,
same split) once there's something real to put in them.

### Where data actually lives

This repo is code only (see Data handling above) — everything the code
reads or writes lives in the workspace bucket, one folder per stage,
mirroring the repo layout 1:1. Everything on the workbench is already
access-controlled, so there's no split by sensitivity — egress review
happens once, at the point something actually leaves the workbench, not
at the folder level.

```
~/workspace/<bucket resource>/data/
├── 01_ancestry_filtering/
│   ├── 1000g_hm3/      # merged whole-genome bfile + HM3 snplist
│   └── round1_pca/     # PCA loadings/scores
├── 02_phenotype/
│   └── modeling_tables/  # one neat TSV per phenotype (prepare_modeling_tables() output) --
│                         # person_id, phenotype, phenotype__invnorm, age, sex_at_birth + covariates,
│                         # so retuning the residualization procedure doesn't need re-pulling from BigQuery
├── 03_grm_shards/
└── 04_process_shards/
```

Bucket paths under `data/<stage>/` are the source of truth. On Verily
Workbench, the workspace bucket is auto-mounted at `~/workspace/<resource
name>/` — tools read and write there directly as an ordinary local path, no
manual mount step needed (`ls ~/workspace/` to find the resource name; `wb
resource list` / `env | grep WORKBENCH` also work). Fall back to plain
local disk + explicit `gsutil cp`/`rsync` only for tools that genuinely need
low-level filesystem access the mount doesn't support well (heavy random
I/O); local disk is ephemeral either way and disappears if the environment
is deleted.

## GRM-pairs submodule

```bash
git submodule update --init --recursive
cd GRM-pairs/grm_bin_sharded && make   # build grm_shard_tool for whatever architecture you're on
```

Pinned to whatever commit `git submodule status` shows — bump with
`git submodule update --remote GRM-pairs` and commit the pointer change.
