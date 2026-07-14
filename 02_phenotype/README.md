# 02_phenotype

Download phenotypes/covariates from AoU, normalize (Kemper et al. 2021
style: residualize → trim outliers → standardize within sex — see main
README for background).

`notebooks/remote/residualize_phenotypes.ipynb` (IRkernel) /
`residualize_phenotypes.Rmd` (R Markdown, identical content, pick whichever
your environment prefers): takes round 2's ancestry-filtered keep-list and
a phenotype list TSV, and for every phenotype exports one `FID IID Y` file
(matching `GRM-pairs/full_grm_bin/prep_pheno.R`'s expected format) per
combination of:

- raw vs. rank-inverse-normal-transformed
- covariate-set (base = sex-at-birth + age; PCs / 3-digit zip factor / SES
  vars each independently toggled on top)

The residualization/transform/export/diagnostics logic lives in
`scripts/local/residualize_lib.R`, shared between the real notebook and
`notebooks/local/test_residualize_fake_data.ipynb` — a synthetic-data smoke
test (fake IDs, fake phenotypes, fake PCs, no AoU access needed) that
exercises the full pipeline end to end. Run that first if changing
anything in the lib.

`docs/phenotype_list.tsv`: a starter real phenotype list, anthropometric +
metabolic panel (height, weight, BMI, systolic/diastolic BP, glucose,
HbA1c, HDL/LDL cholesterol, triglycerides), OMOP concept_ids looked up
against public AoU/OHDSI documentation. 7 of 10 concept_ids are confirmed;
the 3 lipid values have a real LOINC code but an unverified concept_id
(marked `UNCONFIRMED` in the `concept_id` column) — check the AoU Data
Browser for the exact concept_id before running those. These are all
public, standard vocabulary identifiers describing *which* concepts to
pull — not participant data, fine to have in git.

**Not yet filled in**, needs the real workbench to pin down:
- the actual BigQuery pulls per phenotype source (`measurement` / `survey`
  / `condition`)
- `zip3_ses_map`'s join path (likely via `observation`, not directly on
  `person`)
- age definition — currently written as age at each phenotype's own
  measurement/response time, since AoU doesn't have a single UKB-style
  assessment visit
