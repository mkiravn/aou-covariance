# Phenotype normalization protocol

Following Kemper et al. 2021 (*"Phenotypic covariance across the entire
spectrum of relatedness for 86 billion pairs of individuals"*, Nature
Communications, PMC7886899), adapted for AoU's available covariates.

## Order of operations (this is not interchangeable — see paper's Methods)

1. **Residualize**: fit `phenotype ~ covariates`, keep residuals.
2. **Trim outliers**: exclude residuals more than 5 SD from the mean.
3. **Standardize within sex**: normalize the surviving residuals to mean 0,
   variance 1, separately within each sex. This is the final phenotype.

## Covariates (UKB original -> AoU adaptation)

| Kemper (UKB) | AoU adaptation | Status |
|---|---|---|
| sex | sex | direct |
| age at assessment | age | direct |
| year of birth | — | dropped; age already captures this for AoU's single assessment wave |
| genotyping batch (106 levels) | — | dropped; see note below |
| PC1-PC25 (array, unrelated individuals) | ancestry/PCA PCs | being computed in a separate round — not wired in yet |
| birth contemporary group (378 UK local authority areas) | `zip3_ses_map` (median income, poverty fraction, deprivation index, etc., per zip3) | see below |

### Why no batch covariate

Checked against a real published AoU pipeline
([asgarilab/SDoH-GeneticRisk-Biobank-MCA](https://github.com/asgarilab/SDoH-GeneticRisk-Biobank-MCA)) --
zero uses of any batch/center/platform covariate anywhere in their QC or
models. AoU's harmonized single-platform genotyping (and joint-called,
harmonized WGS processing) appears to make this a non-issue in practice,
unlike UK Biobank's two-array system.

### Geographic factor: `zip3_ses_map`

Rather than building a geographic cluster from scratch (Kemper's approach,
not really replicable given AoU's zip3-only geographic granularity), use
AoU's own `zip3_ses_map` BigQuery table directly as covariates:
`median_income`, `poverty` (fraction), `deprivation_index`,
`high_school_education`, `no_health_insurance`, `vacant_housing`,
`assisted_income`. Confirmed in use by the same real AoU pipeline above
(`g_extract_surveys_and_get_mca.ipynb`). Exact join (zip3 -> person_id) not
yet worked out here — TODO.

### Zip3 caveats

- Controlled Tier gives 3-digit zip only; Registered Tier gives state only.
- Any zip3 with <20,000 residents (Census) gets randomly reassigned to a
  neighboring zip3.
- Zip3 is fully suppressed (`000*`) for participants who self-identify as
  American Indian/Alaska Native (CDR v7/v8) — the geographic covariate is
  unusable for this subgroup specifically.

## Open items

- [ ] Wire in ancestry PCs once the separate PC-computation round is done.
- [ ] Confirm the exact `zip3_ses_map` join path (likely via `observation`,
      not directly on `person`).
- [ ] Decide what counts as "age at assessment" given AoU doesn't have a
      single UKB-style assessment visit.
