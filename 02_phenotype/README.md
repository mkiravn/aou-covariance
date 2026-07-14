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

`pull_phenotype()` is implemented for `source == "measurement"` via the
`allofus` package (`aou_connect()` + `aou_sql()`), most recent value per
person, joined to age from `person` and sex from `person.
sex_at_birth_concept_id` — a direct AoU-specific column, distinct from
`person.gender_concept_id` (gender identity). Age is computed from a
single fixed `REFERENCE_DATE` (set per run) rather than at each
phenotype's own measurement time — simpler, and matches validated
real-world usage.

`pull_covariates()`'s zip3/SES join is implemented against
`zip3_ses_map`: 3-digit zip comes from a masked (privacy-protected)
`observation` row (marked with a `*` in `value_as_string`), joined on
zip3 to get `median_income`, `poverty`, `deprivation_index`.

Raw `pull_phenotype()` output is cached per phenotype under
`RAW_PHENO_CACHE_DIR` (one TSV per phenotype), so re-running the notebook
while iterating on `residualize_lib.R` doesn't re-hit BigQuery — delete a
phenotype's cache file to force a refresh. Any covariate-set combo whose
covariates are entirely `NA` for a given phenotype is skipped rather than
crashing the whole run; check `combo_summary_table$status` for which
combos actually ran.

**Not yet filled in**, needs the real workbench to pin down:
- `survey` / `condition` phenotype sources (only `measurement` is wired up)
