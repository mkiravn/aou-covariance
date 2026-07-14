library(dplyr)

# Rank-based inverse-normal transform. Chosen over log/Box-Cox since it
# doesn't assume a particular skew direction or require positive values.
inverse_normal_transform <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  qnorm((r - 0.5) / n)
}

add_transformed_variant <- function(df, pheno_col) {
  df[[paste0(pheno_col, "__invnorm")]] <- inverse_normal_transform(df[[pheno_col]])
  df
}

skew_summary <- function(x, label) {
  tibble(variant = label, n = sum(!is.na(x)), skewness = e1071::skewness(x, na.rm = TRUE))
}

# Named list of covariate-set formula RHS vectors. sex_at_birth is handled
# separately (the stratification variable for step 3, not a residualization
# covariate -- see residualize_phenotype()), so it isn't in these formulas.
build_covariate_sets <- function(pc_cols) {
  list(
    base              = c("age"),
    base_pcs          = c("age", pc_cols),
    base_pcs_zip3     = c("age", pc_cols, "zip3"),
    base_pcs_ses      = c("age", pc_cols, "median_income", "poverty", "deprivation_index"),
    base_pcs_zip3_ses = c("age", pc_cols, "zip3", "median_income", "poverty", "deprivation_index")
  )
}

# Protocol (order matters, see Kemper et al. 2021 Methods):
#   1. residualize phenotype ~ covariates
#   2. trim residuals > outlier_sd SD from the mean
#   3. standardize surviving residuals to mean 0, variance 1, within each sex
residualize_phenotype <- function(df, pheno_col, sex_col, covariate_cols, outlier_sd = 5) {
  formula <- as.formula(paste(pheno_col, "~", paste(covariate_cols, collapse = " + ")))
  fit <- lm(formula, data = df, na.action = na.exclude)
  df$.residual <- residuals(fit)

  keep <- !is.na(df$.residual)
  m <- mean(df$.residual[keep])
  s <- sd(df$.residual[keep])
  not_outlier <- keep & abs(df$.residual - m) <= outlier_sd * s

  df$phenotype_norm <- NA_real_
  for (s_level in unique(df[[sex_col]][not_outlier])) {
    idx <- not_outlier & df[[sex_col]] == s_level
    r <- df$.residual[idx]
    df$phenotype_norm[idx] <- (r - mean(r)) / sd(r)
  }

  list(
    data = df %>% select(-.residual),
    n_input = nrow(df),
    n_retained = sum(!is.na(df$phenotype_norm)),
    r_squared = summary(fit)$r.squared
  )
}

# Matches GRM-pairs/full_grm_bin/prep_pheno.R's expected input format.
write_grm_pheno <- function(df, file) {
  out <- df %>%
    filter(!is.na(phenotype_norm)) %>%
    transmute(FID = person_id, IID = person_id, Y = phenotype_norm)
  write.table(out, file = file, quote = FALSE, sep = " ",
              row.names = FALSE, col.names = TRUE, na = "NA")
}

# Runs the full phenotype x {raw, invnorm} x covariate-set cross product,
# writing one .pheno file per combination and returning diagnostics.
# `pull_phenotype(row, keep_ids)` and `pull_covariates(keep_ids)` are
# supplied by the caller -- real AoU pulls in the remote notebook, synthetic
# generators in the local fake-data test.
run_residualization <- function(pheno_list, keep_ids, pull_phenotype, pull_covariates,
                                 covariate_sets, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  skew_diagnostics <- list()
  combo_diagnostics <- list()

  for (i in seq_len(nrow(pheno_list))) {
    row <- pheno_list[i, ]
    name <- row$phenotype_name

    pheno_df <- pull_phenotype(row, keep_ids)
    covars <- pull_covariates(keep_ids)

    df <- pheno_df %>%
      left_join(covars$pcs, by = "person_id") %>%
      left_join(covars$zip3, by = "person_id") %>%
      left_join(covars$ses, by = "person_id") %>%
      add_transformed_variant("phenotype")

    skew_diagnostics[[paste0(name, "__raw")]] <-
      skew_summary(df$phenotype, "raw") %>% mutate(phenotype = name, .before = 1)
    skew_diagnostics[[paste0(name, "__invnorm")]] <-
      skew_summary(df[["phenotype__invnorm"]], "invnorm") %>% mutate(phenotype = name, .before = 1)

    for (variant_col in c("phenotype", "phenotype__invnorm")) {
      variant_label <- ifelse(variant_col == "phenotype", "raw", "invnorm")

      for (covset_name in names(covariate_sets)) {
        result <- residualize_phenotype(df, variant_col, "sex_at_birth", covariate_sets[[covset_name]])

        out_name <- sprintf("%s__%s__%s", name, variant_label, covset_name)
        write_grm_pheno(result$data, file.path(out_dir, paste0(out_name, ".pheno")))

        combo_diagnostics[[out_name]] <- tibble(
          combo = out_name, phenotype = name, variant = variant_label,
          covariate_set = covset_name,
          n_input = result$n_input, n_retained = result$n_retained,
          r_squared = result$r_squared
        )
      }
    }
  }

  list(
    skew_summary_table = bind_rows(skew_diagnostics),
    combo_summary_table = bind_rows(combo_diagnostics)
  )
}
