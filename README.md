# aou-covariance

Phenotypic covariance / heritability pipeline on All of Us, following the
Kemper et al. 2021 design (*"Phenotypic covariance across the entire
spectrum of relatedness for 86 billion pairs of individuals"*) adapted to
the AoU Researcher Workbench. Uses [GRM-pairs](https://github.com/mkiravn/GRM-pairs)
(vendored as a submodule) for the sparse/binned GRM tooling and the
individual-block jackknife.

## Pipeline stages

| Stage | Purpose |
|---|---|
| `01_ancestry_filtering/` | Restrict to a genetically homogeneous reference cluster via a PCA ellipsoid filter, validated against 1000G (see its own docs). |
| `02_phenotype/` | Pull phenotypes + covariates from AoU, residualize, trim outliers, standardize within sex — the Kemper protocol. |
| `03_grm_shards/` | Compute the GRM in row-chunked shards (`plink --make-grm-bin --parallel k n`) as AoU batch jobs, store shards in the workspace bucket. |
| `04_process_shards/` | Bin each shard by relatedness and accumulate phenotype cross-products (`GRM-pairs/grm_bin_sharded`), merge, get bin means + jackknife SEs. |

The final analysis ID list is the intersection of: passes ancestry
filtering (01), has a usable phenotype after normalization (02), and has
genotype data included in the GRM (03/04).

## Layout convention

Each numbered stage has its own `docs/` (parameters, validation results,
decisions) and `scripts/` (`local/` for anything prototyped outside the
workbench, `aou/` for what actually runs there — dsub job scripts, workbench
notebooks). Keep genomic data files, GRM binaries, and large intermediate
results out of git — see `.gitignore`.

## GRM-pairs submodule

```bash
git submodule update --init --recursive
cd GRM-pairs/grm_bin_sharded && make   # build grm_shard_tool for whatever architecture you're on
```

Pinned to whatever commit `git submodule status` shows — bump with
`git submodule update --remote GRM-pairs` and commit the pointer change.
[GRM-pairs#16](https://github.com/mkiravn/GRM-pairs/pull/16) and
[#17](https://github.com/mkiravn/GRM-pairs/pull/17) (untracking every
accidentally-committed compiled binary in that repo) are both merged, so
the pin here is clean.

## Covariate decisions so far

- **Geographic factor**: AoU's `zip3_ses_map` table (per-zip3 SES: median
  income, poverty fraction, deprivation index, etc.) rather than building a
  Kemper-style geographic cluster from scratch — see `02_phenotype/docs/`.
- **Genotyping/sequencing batch**: not used. AoU's harmonized single-platform
  processing appears to make this unnecessary in practice (checked against
  a real published AoU pipeline — no batch covariate anywhere in it).
- **Ancestry PCs**: computed in a separate round, not detailed here yet.
