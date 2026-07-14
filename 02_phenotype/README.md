# 02_phenotype

Download phenotypes/covariates from AoU, normalize (Kemper et al. 2021
style: residualize → trim outliers → standardize within sex — see main
README for background).

`notebooks/remote/residualize_phenotypes.ipynb`: takes round 2's
ancestry-filtered keep-list and a phenotype list TSV, and for every
phenotype exports one `FID IID Y` file (matching
`GRM-pairs/full_grm_bin/prep_pheno.R`'s expected format) per combination of:

- raw vs. rank-inverse-normal-transformed
- covariate-set (base = sex-at-birth + age; PCs / 3-digit zip factor / SES
  vars each independently toggled on top)

The residualization/transform/export/diagnostics logic is real and runnable.
**Not yet filled in**, needs the real workbench to pin down:
- the actual BigQuery pulls per phenotype source (`measurement` / `survey`
  / `condition`)
- `zip3_ses_map`'s join path (likely via `observation`, not directly on
  `person`)
- age definition — currently written as age at each phenotype's own
  measurement/response time, since AoU doesn't have a single UKB-style
  assessment visit
