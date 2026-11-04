# Round-2 PCA (Privé et al. 2020): bed_autoSVD on 1000G founders of TRAIN_POPS,
# round-1 participants projected by OADP.
# env: REF_BFILE REF_LABELS TRAIN_POPS NEW_BED K LD_REGIONS OUT_DIR NCORES
suppressMessages({
  library(bigsnpr)
  library(data.table)
})
env <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("missing env ", k)
  v
}
out_dir <- env("OUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
K <- as.integer(env("K"))
ncores <- as.integer(Sys.getenv("NCORES", "1"))

obj.ref <- bed(paste0(env("REF_BFILE"), ".bed"))
obj.new <- bed(env("NEW_BED"))

lab <- fread(env("REF_LABELS"), colClasses = "character")
lab <- lab[match(obj.ref$fam$sample.ID, lab$SampleID)]
stopifnot(!anyNA(lab$SampleID))
train_pops <- strsplit(env("TRAIN_POPS"), ",")[[1]]
founder <- lab$FatherID == "0" & lab$MotherID == "0"
ind.train <- which(lab$Population %in% train_pops & founder)
cat(sprintf("training on %d founders of %s (%d non-founders left out)\n",
            length(ind.train), paste(train_pops, collapse = ","),
            sum(lab$Population %in% train_pops & !founder)))

# listed long-range LD regions, excluded before autoSVD looks for more
ld <- fread(env("LD_REGIONS"), header = FALSE)
ld[, V1 := sub("^chr", "", V1)]
map <- obj.ref$map
in_ld <- rep(FALSE, nrow(map))
for (r in seq_len(nrow(ld))) {
  in_ld <- in_ld | (as.character(map$chromosome) == ld$V1[r] &
                      map$physical.pos >= ld$V2[r] & map$physical.pos <= ld$V3[r])
}

proj <- bed_projectPCA(obj.ref, obj.new, k = K,
                       ind.row.ref = ind.train, ind.col.ref = which(!in_ld),
                       strand_flip = TRUE, join_by_pos = TRUE,
                       ncores = ncores, verbose = TRUE)
svd <- proj$obj.svd.ref
used <- attr(svd, "subset")
lrldr <- attr(svd, "lrldr")
pcs <- paste0("PC", seq_len(K))

ref_scores <- predict(svd)
colnames(ref_scores) <- pcs
oadp <- proj$OADP_proj
simple <- proj$simple_proj
colnames(oadp) <- colnames(simple) <- pcs

scores <- rbind(
  data.table(id = lab$SampleID[ind.train], set = "reference",
             pop = lab$Population[ind.train], ref_scores),
  data.table(id = obj.new$fam$sample.ID, set = "cohort", pop = NA_character_, oadp))
fwrite(data.table(id = obj.new$fam$sample.ID, simple),
       file.path(out_dir, "simple_projection.tsv"), sep = "\t")
fwrite(data.table(chr = map$chromosome[used], pos = map$physical.pos[used], svd$v),
       file.path(out_dir, "loadings.tsv"), sep = "\t")
fwrite(data.table(pc = seq_along(svd$d), d = svd$d), file.path(out_dir, "singular_values.tsv"), sep = "\t")
fwrite(if (is.null(lrldr) || !nrow(lrldr)) data.table(Chr = integer(), Start = integer(), Stop = integer())
       else as.data.table(lrldr), file.path(out_dir, "lrldr.tsv"), sep = "\t")
jsonlite::write_json(list(
  n_train = length(ind.train), n_cohort = nrow(obj.new$fam),
  n_ref_variants = ncol(obj.ref), n_listed_ld_excluded = sum(in_ld),
  n_variants_used = length(used), n_lrldr_regions = if (is.null(lrldr)) 0L else nrow(lrldr),
  k = K), file.path(out_dir, "round2_meta.json"), auto_unbox = TRUE, pretty = TRUE)
# the stage target last: its presence means everything above was written
fwrite(scores, file.path(out_dir, "scores.tsv"), sep = "\t")
cat(sprintf("%d variants used after autoSVD; %d long-range LD region(s) detected\n",
            length(used), if (is.null(lrldr)) 0L else nrow(lrldr)))
