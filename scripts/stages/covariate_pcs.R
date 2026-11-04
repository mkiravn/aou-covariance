# Covariate PCs: fit on the unrelated part of the final set, everyone else
# projected by OADP. Nobody is dropped.
# env: BED KEEP RELATEDNESS KIN_CUT K THR_R2 LD_REGIONS OUT NCORES
suppressMessages({
  library(bigsnpr)
  library(data.table)
})
env <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("missing env ", k)
  v
}
K <- as.integer(env("K"))
dir.create(dirname(env("OUT")), recursive = TRUE, showWarnings = FALSE)
ncores <- as.integer(Sys.getenv("NCORES", "1"))
obj <- bed(env("BED"))
ids <- as.character(obj$fam$sample.ID)
final <- fread(env("KEEP"), header = FALSE, colClasses = "character")$V1
rows <- which(ids %in% final)

rel <- fread(env("RELATEDNESS"), colClasses = list(character = c("i.s", "j.s")))
rel <- rel[kin > as.numeric(env("KIN_CUT")) & i.s %in% final & j.s %in% final]
related <- unique(c(rel$i.s, rel$j.s))
fit <- rows[!(ids[rows] %in% related)]
rest <- setdiff(rows, fit)
cat(sprintf("%d in the final set: fitting on %d, projecting %d with a relative in the set\n",
            length(rows), length(fit), length(rest)))

ld <- fread(env("LD_REGIONS"), header = FALSE)
ld[, V1 := sub("^chr", "", V1)]
map <- obj$map
in_ld <- rep(FALSE, nrow(map))
for (r in seq_len(nrow(ld))) {
  in_ld <- in_ld | (as.character(map$chromosome) == ld$V1[r] &
                      map$physical.pos >= ld$V2[r] & map$physical.pos <= ld$V3[r])
}

svd <- bed_autoSVD(obj, ind.row = fit, ind.col = which(!in_ld), k = K,
                   thr.r2 = as.numeric(env("THR_R2")), ncores = ncores)
pcs <- paste0("PC", seq_len(K))
pc_fit <- predict(svd)
pc_rest <- if (length(rest)) {
  bed_projectSelfPCA(svd, obj, ind.row = rest, ncores = ncores)$OADP_proj
} else matrix(numeric(), 0, K)
colnames(pc_fit) <- colnames(pc_rest) <- pcs
out <- rbind(data.table(person_id = ids[fit], in_fit = TRUE, pc_fit),
             data.table(person_id = ids[rest], in_fit = FALSE, pc_rest))
lrldr <- attr(svd, "lrldr")
jsonlite::write_json(list(
  n_final = length(rows), n_fit = length(fit), n_projected = length(rest),
  n_variants_used = length(attr(svd, "subset")),
  n_lrldr_regions = if (is.null(lrldr)) 0L else nrow(lrldr),
  singular_values = svd$d), sub("\\.tsv$", "_meta.json", env("OUT")),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(out, env("OUT"), sep = "\t")
