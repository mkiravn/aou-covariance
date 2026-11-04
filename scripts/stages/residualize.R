# Phenotype cleaning with residualize_lib.R (protocol as in 05_phenotype_residualize).
# Writes FID = 0 to match plink's .grm.id.
# env: LIB TABLE_DIR PHENO_LIST PHENOS TRANSFORMS COVSETS KEEP PCS N_PCS OUTLIER_SD OUT_DIR
suppressMessages({
  library(data.table)
  library(dplyr)
})
env <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("missing env ", k)
  v
}
source(env("LIB"))
split <- function(k) strsplit(env(k), ",")[[1]]
out_dir <- env("OUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- fread(env("KEEP"), header = FALSE, colClasses = "character")$V1
pcs <- fread(env("PCS"), colClasses = list(character = "person_id"))
pc_cols <- paste0("PC", seq_len(as.integer(env("N_PCS"))))
covsets <- list(base = "age", base_pcs = c("age", pc_cols))
plist <- fread(env("PHENO_LIST"), sep = "\t")
outlier_sd <- as.numeric(env("OUTLIER_SD"))

qc <- list()
for (ph in split("PHENOS")) {
  row <- plist[phenotype_name == ph]
  stopifnot(nrow(row) == 1)
  tab <- fread(file.path(env("TABLE_DIR"), paste0(ph, ".tsv")),
               colClasses = list(character = "person_id"))
  n_table <- nrow(tab)
  tab <- tab[person_id %in% keep]
  n_keep <- nrow(tab)
  pr <- filter_plausible_range(as.data.frame(tab), "phenotype", row$plausible_min, row$plausible_max)
  df <- merge(as.data.table(pr$data), pcs, by = "person_id")
  df <- as.data.frame(add_transformed_variant(as.data.frame(df), "phenotype"))
  for (tr in split("TRANSFORMS")) {
    col <- if (tr == "raw") "phenotype" else "phenotype__invnorm"
    for (cs in split("COVSETS")) {
      res <- residualize_phenotype(df, col, "sex_at_birth", covsets[[cs]], outlier_sd)
      out <- res$data %>% filter(!is.na(phenotype_norm)) %>%
        transmute(FID = 0, IID = person_id, Y = phenotype_norm)
      write.table(out, file.path(out_dir, sprintf("%s__%s__%s.pheno", ph, tr, cs)),
                  quote = FALSE, sep = " ", row.names = FALSE)
      qc[[length(qc) + 1]] <- data.table(
        phenotype = ph, transform = tr, covset = cs,
        n_table = n_table, n_in_set = n_keep, n_implausible = pr$n_excluded,
        n_with_pcs = nrow(df), n_retained = res$n_retained,
        n_distinct = length(unique(df[[col]])), r_squared = res$r_squared,
        skew_input = e1071::skewness(df[[col]], na.rm = TRUE))
    }
  }
}
fwrite(rbindlist(qc), file.path(out_dir, "qc_summary.tsv"), sep = "\t")
print(rbindlist(qc))
