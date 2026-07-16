# 03_grm_shards

Compute the GRM in row-chunked shards (`plink --make-grm-bin --parallel k n`)
as AoU batch jobs, store shards in the workspace bucket.

`notebooks/remote/chr22_qc_thinning_timing.ipynb`: first look before committing
to a genome-wide run — QC (min MAF 1%, HWE 1e-6, missingness < 5%, biallelic
only) on ACAF's chr22, restricted to round 2b's ancestry-filtered keep-list,
then random thinning (`plink2 --thin`) tuned toward a ~1M-variant genome-wide
target.

ACAF's pvar ships with the variant ID column unset (`.` for every variant), so
`--rm-dup` has nothing to compare IDs against unless `--set-all-var-ids
'@:#:$r:$a'` runs first (same convention as `submit_pca_r1.ipynb`'s ACAF
handling) — without it, `--rm-dup` silently removes 0 variants regardless of
how many duplicates actually exist.

QC and thinning are separate plink2 calls specifically so re-tuning the
thinning probability doesn't require re-running the (expensive) QC pass.
Every plink2 call is timed (`time`) — chr22 is ~1.6% of the autosomal genome
by length, so this is the seed estimate for genome-wide QC wall-clock time.
QC'd and thinned pgen sets (plus their `.log`s, which carry the exact filter
counts) get copied from local scratch to `data/03_grm_shards/` in the bucket
at the end — local scratch isn't guaranteed to survive a session restart.

`notebooks/remote/genome_wide_qc_thinning.ipynb`: all 22 autosomes, same QC +
`--thin 0.2` (chr22's calibrated value, fixed rather than recomputed per
chromosome), run interactively rather than as submitted jobs — a fully serial
22-chromosome run extrapolates from chr22's ~10 min QC time to ~9.4 hours, so
chromosomes run concurrently instead (`ThreadPoolExecutor`-managed
`plink2` subprocesses, capped at `N_CONCURRENT`, biggest chromosomes first so
no batch is left waiting on a single straggler). Reports per-chromosome
timing/counts and the effective speedup vs. a naive serial estimate — a
speedup well below `N_CONCURRENT` means the shared network-mounted ACAF reads
are contending across concurrent chromosomes, which argues for moving to job
submission (below) rather than pushing concurrency higher in a single
session.

**Not yet started**: job submission (`dsub`, one Google Batch task per
chromosome — each gets its own machine and network path, unlike the
interactive notebook's shared session; the ACAF `gs://` bucket URI and
whether Controlled Tier VPC-SC rules permit a Batch worker to read it are
both still open, to confirm before writing the `dsub` task script), GRM shard
computation, `GRM-pairs/grm_bin_sharded` parallelization tuning.

**Before running `grm_shard_tool` against `02_phenotype`'s `.pheno`
output**: format is already checked and compatible (see
`02_phenotype/README.md`'s "GRM-pairs compatibility" note) — the one open
risk is that `grm_shard_tool`'s pheno lookup keys on the full `(FID, IID)`
pair, and `write_grm_pheno()` sets `FID = IID = person_id`, while a plink
`.grm.id` commonly has `FID = "0"` by default. Check the real `.grm.id`'s
`FID` column once this stage produces one before trusting a run.
