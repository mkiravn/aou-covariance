# 02_phenotype

Download phenotypes/covariates from AoU, normalize (trim implausible
values → Kemper et al. 2021 style: residualize → trim outliers →
standardize within sex — see main README for background).

`notebooks/remote/query_filter_check.ipynb`: a smaller smoke test to run
first, against the real CDR, before trusting the full pipeline below — pulls
just the 3 fully-confirmed phenotypes from `docs/phenotype_list.tsv` (no
`UNCONFIRMED` lipid rows) and checks concept_id validity, pull row/person
counts, the keep-list filtering funnel, sex_at_birth breakdown, and value
ranges. Every cell prints aggregate counts/summary stats only, never a
person-level row. Connects via `allofus::aou_connect()` / `aou_sql()`, same
as the main pipeline below — confirmed working on Workbench 2.0 (Verily) in
practice. Reads both round 2 (`round2_filter.ipynb`, 1000G-fit ellipsoid,
the default) and round 2b (`reverse_pca_aou.ipynb`, AoU-fit ellipsoid,
provisional) keep-lists, so either can be checked against — flip
`KEEP_LIST_PATH` to compare. Step 5 runs `residualize_lib.R`'s real
`run_residualization()` — the same function `residualize_phenotypes.ipynb`
calls, wrapped around the data already pulled in earlier steps instead of
issuing new BigQuery calls — across all 4 covariate-set combos (round 2b's
`PC1..PC5` for PCs; round 2b's file has 20, only the top 5 are used, since
beyond that isn't considered informative for this cohort) crossed with
`{raw, invnorm}`, exercising the pipeline's actual statistical step against
real AoU values, not just synthetic data like
`test_residualize_fake_data.ipynb`. `.pheno`-shaped files land in a
throwaway `/tmp` directory. Step 6 goes deeper: pulls SES data (same
`zip3_ses_map` join as `pull_covariates()`, zip3 itself included this
time), fits `base_pcs_zip3_ses` (the richest combo) directly with `lm()` to
inspect the actual coefficients (how each phenotype loads onto
age/PCs/zip3/SES), checks covariate encoding (numeric types
post-integer64-fix, PC mean/SD from `variance-standardize` scoring,
`sex_at_birth` category counts, zip3 as a one-hot-encoded factor —
summarized as dummy-term/significant counts rather than printing all ~800
possible zip3 coefficients), and reports phenotype distributions
(histograms, skewness before/after `inverse_normal_transform()`, quantiles,
by-sex boxplots) — still model-level/aggregate output only, never a
person-level row.

`notebooks/remote/residualize_phenotypes.ipynb` (IRkernel) /
`residualize_phenotypes.Rmd` (R Markdown, identical content, pick whichever
your environment prefers): takes round 2's ancestry-filtered keep-list and
a phenotype list TSV, and for every phenotype exports one `FID IID Y` file
(matching `GRM-pairs/full_grm_bin/prep_pheno.R`'s expected format) per
combination of:

- raw vs. rank-inverse-normal-transformed
- covariate-set: `build_covariate_sets()` is a nested staircase, not
  independently-toggled combinations — `base` (sex-at-birth + age) →
  `base_pcs` (+ PCs) → `base_pcs_zip3` (+ 3-digit zip factor) →
  `base_pcs_zip3_ses` (+ SES vars)

The residualization/transform/export/diagnostics logic lives in
`scripts/local/residualize_lib.R`, shared between the real notebook and
`notebooks/local/test_residualize_fake_data.ipynb` — a synthetic-data smoke
test (fake IDs, fake phenotypes, fake PCs, no AoU access needed) that
exercises the full pipeline end to end. Run that first if changing
anything in the lib.

Split into two stages, each its own cell in the real notebook:
`prepare_modeling_tables()` (pulls, applies the plausible-range filter,
joins PCs/zip3/SES, adds the invnorm variant, writes one neat TSV per
phenotype — `<phenotype_name>.tsv`: `person_id, phenotype,
phenotype__invnorm, age, sex_at_birth`, plus every covariate column — to
`MODELING_TABLE_DIR`) is the only step that touches BigQuery.
`run_residualization_from_tables()` reads those TSVs back in and runs the
`{raw, invnorm} x covariate-set` cross product. `MODELING_TABLE_DIR` should
live in the workspace bucket (`data/02_phenotype/modeling_tables/`, per
root README's bucket layout), not local disk — the whole point is that
retuning the residualization procedure itself (`covariate_sets`,
`outlier_sd`, which phenotypes) only needs re-running the second stage, not
re-pulling from BigQuery. `run_residualization()` still exists as a
convenience wrapper that calls both stages back-to-back (writing tables to
an ephemeral `tempfile()` dir) — what `test_residualize_fake_data.ipynb`
and `query_filter_check.ipynb`'s smoke tests use, since they don't need the
tables to persist.

`docs/phenotype_list.tsv`: the real phenotype list — anthropometric +
metabolic panel (height, weight, BMI, systolic/diastolic BP, glucose,
HbA1c, HDL/LDL/total cholesterol, triglycerides, waist/hip circumference,
heart rate, creatinine, hemoglobin, WBC count, platelet count, ALT,
waist_hip_ratio). Concept_ids for the rows promoted from
`docs/candidate_phenotypes.tsv` were verified directly against the real
CDR (the standard `concept`/`concept_relationship` "Maps to" join from
each LOINC code, plus an `n_persons` sanity check against
`{CDR}.measurement`) — `hip_circumference` is still `UNCONFIRMED`, since
its LOINC (8286-7) turned out to map to *thigh* circumference, not hip
(`n_rows = 0`), caught exactly by this check. `waist_circumference`'s
concept_id resolved correctly but has only `n_persons = 22` in the whole
CDR — confirmed correct, essentially unusable — so it, `hip_circumference`,
and `waist_hip_ratio` are flagged in their `notes` as deprioritized rather
than actively pursued further. These are all public, standard vocabulary
identifiers describing *which* concepts to pull — not participant data,
fine to have in git.

`residualize_phenotypes.ipynb`/`.Rmd` filter `UNCONFIRMED` rows out of
`pheno_list` right after reading it (with a `message()` naming which ones
got skipped) — `pull_phenotype()` also refuses to run on a non-numeric
`concept_id` as a second line of defense, so a bad row fails with a clear
message instead of building `WHERE measurement_concept_id IN (UNCONFIRMED)`
and failing deep inside BigQuery with an opaque `aou_sql()` error. Beyond
that, `prepare_modeling_tables()`/`run_residualization_from_tables()` (in
`residualize_lib.R`) catch a bad row at every remaining failure point —
`pull_phenotype()` throwing (wrong/unimplemented `source`, a `concept_id`
that looks numeric but doesn't actually exist in this CDR version, a
transient BigQuery error) or a missing modeling table — and skip just that
phenotype with a `message()` and a `status` column recording why, rather
than taking the whole run down. `pull_covariates()` is deliberately *not*
given this treatment: it's the same call for every phenotype, so a failure
there is systemic, not phenotype-specific, and should stop the run loudly.
`test_residualize_fake_data.ipynb` exercises this with a `fake_broken`
phenotype whose pull always throws, confirming the other fake phenotypes
still complete.

`docs/candidate_phenotypes.tsv`: the same curated wishlist this panel grew
from — same schema plus `expected_completeness` (qualitative) and
`approx_snp_h2` (rough published SNP-heritability, general-population
GWAS/twin literature, *not* AoU-specific — for prioritization, not a claim
about what this pipeline will actually find), kept around for the rows not
yet promoted: 8 Fitbit-derived phenotypes (mean daily steps, sedentary
minutes, very-active minutes, calories out, sleep minutes, minutes awake,
minutes in bed, mean heart rate) flagged `source == "fitbit"` — genuinely
distinctive to AoU (longitudinal wearable data linked to WGS + EHR, not
available at this scale in most public biobanks) but a real tradeoff
against completeness, since only participants who shared Fitbit data are
represented. Schema confirmed directly against the real CDR (via
`INFORMATION_SCHEMA.TABLES`/`.COLUMNS`, not the OMOP `measurement` table):
`activity_summary` (daily steps/calories/active-minutes breakdown),
`sleep_daily_summary` (daily sleep-stage minutes), `heart_rate_minute_level`
(raw per-minute heart rate — no clean "resting heart rate" field exists
anywhere in the schema, so `fitbit_mean_heart_rate` averages every
recorded minute regardless of activity level, not a true resting rate).
Each of these tables exists in **two variants** in this CDR — lowercase
with `row_id`/`src_id` columns, and an `UPPERCASE` variant without them —
not yet confirmed which is the live/current one to actually query.

None of this fits `pull_phenotype()`'s "most recent value" shape even once
`source == "fitbit"` is implemented: `activity_summary`/
`sleep_daily_summary` are one row per person *per day*,
`heart_rate_minute_level` one row per person *per minute*, so every
Fitbit candidate needs averaging across a person's available days first —
a real methodological choice, not just an implementation detail. Each
row's notes carry the same open questions: minimum valid days required
per person (wearable studies, e.g. UK Biobank's accelerometer data,
typically enforce a wear-time/day-count threshold so 1-2 noisy days don't
dominate someone's value), mean vs. median across days, and whether to use
a date window near `REFERENCE_DATE` or all available history. Promoting
one means resolving those questions, picking the lowercase-vs-UPPERCASE
table variant, and implementing the `fitbit` source, not just copying a
row into `phenotype_list.tsv`.

2 Lifestyle-survey phenotypes: `alcohol_audit_c_score` (`source ==
"survey_composite"`, AUDIT-C — the standard validated 3-item alcohol-use
score, `concept_id` holding all 3 underlying question IDs) and
`cigarettes_per_day` (`source == "survey"`, from the Lifestyle survey's
TUS-CPS-derived tobacco items). AoU's core survey battery has no dedicated
diet/nutrition module (checked the [Survey
Explorer](https://researchallofus.org/survey-explorer) directly — none of
the 8 modules cover food intake or a food-frequency questionnaire, unlike
e.g. UK Biobank), so these are the closest fit to "lifestyle phenotypes"
actually available. Neither `source` is wired up (`survey_composite` would
need a new code path entirely — pulling multiple survey answers from
`{CDR}.observation`, not `{CDR}.measurement`, and summing them, similar in
spirit to `waist_hip_ratio`'s `derived_ratio` but summing instead of
dividing). `cigarettes_per_day` is only answered by participants who
report current/former smoking, so its completeness is meaningfully lower
than the survey's own response rate — flagged as cutting against the "few
missing" goal more than any EHR row above.

`waist_hip_ratio` (`source == "derived_ratio"`, already promoted) is a
different shape from every other row: `concept_id` holds both underlying
concept_ids (`waist,hip`), and it's not a single pulled value but
`waist_circumference / hip_circumference` for the same person.
`pull_phenotype()` still throws on it (unimplemented source, caught and
skipped like any other bad row) — it's built instead as its own cell in
`residualize_phenotypes.ipynb`/`.Rmd`, "Derived phenotype: waist/hip
ratio", which combines the already-prepared `waist_circumference.tsv`/
`hip_circumference.tsv` modeling tables (no new BigQuery calls) rather
than pulling fresh. **Caveat carried into that cell's markdown too:**
`waist`/`hip` are each independently "most recent value as of
`REFERENCE_DATE`", joined on `person_id` — not matched to the same
visit/date. AoU's Physical Measurements module tends to collect both
together, so this is usually fine in practice, but it's worth checking
first if the ratio's distribution looks off.

`plausible_min`/`plausible_max` columns (original units, e.g. cm for
height, mg/dL for glucose): generous physiological-plausibility bounds,
applied by `filter_plausible_range()` in `residualize_lib.R` before any
modeling. This is deliberately not the same thing as Kemper et al. 2021's
existing post-residual 5-SD trim (still in `residualize_phenotype()`,
unchanged) — that trim only catches values extreme *relative to the fitted
model*, computed separately per covariate-set, so a data-entry error that
happens to land near the model's mean can slip through it. EHR-sourced
measurements (unlike a curated cohort assessment like UK Biobank's) carry a
non-trivial rate of exactly this kind of error — AoU's own analysis of
height/weight in the EHR found ~1.4-1.5% flagged as erroneous
([Zhou et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11973958/)) —
so a wide absolute-range filter is a cheap, standard complement. Bounds
here are intentionally generous (e.g. height 100-250cm, not a "normal"
range) — the goal is dropping impossible values, not disease-range ones,
which are real biological variation GWAS wants to keep.

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

Every numeric covariate pulled from BigQuery (`age`, `median_income`,
`poverty`, `deprivation_index`) is explicitly `CAST(... AS FLOAT64)` in the
SQL. Bare `INT64` columns collect via `bigrquery`/`allofus` as
`bit64::integer64`, not a plain double, and `lm()` silently mis-coerces
that into a degenerate fit (`r_squared = NaN`, everything trimmed as an
outlier) rather than erroring — caught by `query_filter_check.ipynb`'s
mock residualization step on `age`, fixed the same way everywhere else a
raw `INT64` covariate is pulled in.

Raw `pull_phenotype()` output is cached per phenotype under
`RAW_PHENO_CACHE_DIR` (one TSV per phenotype) — separate from and upstream
of `MODELING_TABLE_DIR`'s prepared tables, this only saves re-hitting
BigQuery within `prepare_modeling_tables()` itself; delete a phenotype's
cache file to force a refresh. Any covariate-set combo whose covariates
are entirely `NA` for a given phenotype is skipped rather than crashing
the whole run; check `combo_summary_table$status` for which combos
actually ran. `prepare_modeling_tables()`'s returned `range_summary_table`
reports how many values `filter_plausible_range()` dropped per phenotype,
before any of that.

**Not yet filled in**, needs the real workbench to pin down:
- `survey` / `condition` phenotype sources (only `measurement` is wired up)
