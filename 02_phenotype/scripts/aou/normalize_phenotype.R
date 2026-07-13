#!/usr/bin/env Rscript
# Kemper et al. 2021 phenotype normalization protocol:
#   1. residualize phenotype on covariates
#   2. drop residuals more than `outlier_sd` SD from the mean
#   3. standardize surviving residuals to mean 0, variance 1, within each sex
#
# This is the generic step -- independent of *how* phenotype/covariates were
# pulled from AoU (see docs/protocol.md for the covariate mapping decisions).

normalize_phenotype <- function(df, pheno_col, sex_col, covariate_cols, outlier_sd = 5) {
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

  df$.residual <- NULL
  df
}

if (sys.nframe() == 0) {
  # quick self-test with synthetic data, run via: Rscript normalize_phenotype.R
  set.seed(1)
  n <- 2000
  df <- data.frame(
    sex = sample(c("Female", "Male"), n, replace = TRUE),
    age = round(runif(n, 18, 80)),
    median_income = rnorm(n, 60000, 15000)
  )
  # sex-specific mean/variance in the raw phenotype, plus age/income effects
  sex_effect <- ifelse(df$sex == "Female", 2, -2)
  df$bmi <- 25 + sex_effect + 0.05 * df$age - 0.00002 * df$median_income + rnorm(n, sd = 4)
  df$bmi[sample(n, 5)] <- df$bmi[sample(n, 5)] + 50  # inject outliers

  out <- normalize_phenotype(df, "bmi", "sex", c("age", "median_income"))

  stopifnot(sum(!is.na(out$phenotype_norm)) < n)  # outliers were dropped
  for (s_level in unique(out$sex)) {
    x <- out$phenotype_norm[out$sex == s_level]
    x <- x[!is.na(x)]
    stopifnot(abs(mean(x)) < 1e-8, abs(sd(x) - 1) < 1e-8)
  }
  cat("normalize_phenotype: self-test passed --",
      sum(is.na(out$phenotype_norm)), "of", n, "dropped as outliers\n")
}
