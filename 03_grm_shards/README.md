# 03_grm_shards

Compute the GRM in row-chunked shards (`plink --make-grm-bin --parallel k n`)
as AoU batch jobs, store shards in the workspace bucket.

`notebooks/remote/chr22_qc_thinning_timing.ipynb`: first look before committing
to a genome-wide run — QC (min MAF 1%, HWE 1e-6, missingness < 5%, biallelic
only) on ACAF's chr22, restricted to round 2b's ancestry-filtered keep-list,
then random thinning (`plink2 --thin`) tuned toward a ~1M-variant genome-wide
target. QC and thinning are separate plink2 calls specifically so re-tuning
the thinning probability doesn't require re-running the (expensive) QC pass.
Every plink2 call is timed (`time`) — chr22 is ~1.6% of the autosomal genome
by length, so this is the seed estimate for genome-wide QC wall-clock time and
for sizing the shard/batch-job parallelization below. QC'd and thinned pgen
sets (plus their `.log`s, which carry the exact filter counts) get copied
from local scratch to `data/03_grm_shards/` in the bucket at the end — local
scratch isn't guaranteed to survive a session restart. Parallelization itself
(`--parallel k n`, `GRM-pairs/grm_bin_sharded`) isn't attempted yet — that's
the next phase, once this notebook's numbers say how big a genome-wide QC job
actually is.

**Not yet started**: everything past chr22 QC — genome-wide QC, GRM shard
computation, AoU batch job submission, `GRM-pairs/grm_bin_sharded`
parallelization tuning.
