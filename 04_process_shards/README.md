# 04_process_shards

Bin each GRM shard by relatedness and accumulate phenotype cross-products
using `GRM-pairs/grm_bin_sharded` (`accumulate`, then `merge`).

`notebooks/remote/grm_shard_processing.ipynb`: shells out to
`grm_shard_tool` (a C++ CLI, not a Python library) to bin the sharded GRM
by relatedness and compute the phenotype cross-product (covariance)
within each bin, producing a covariance-vs-relatedness plot per chosen
phenotype. Designed to run before `03_grm_shards`' shard construction
finishes all shards -- globs for whatever `.grm.bin.k` files are actually
present at its persist destination, uses only those, and reports how many
of the total it found. Rerunning later, as more shards finish, picks up
more coverage automatically, no code changes needed.

Rewrites each phenotype file's `FID` from the real `.grm.id` (keyed on
`IID`) before use -- confirmed necessary for this panel, whose `.grm.id`
has `FID = "0"` uniformly, while `02_phenotype`'s `write_grm_pheno()` sets
`FID = person_id`. Trusting either convention blindly would have silently
dropped every pairing; see `03_grm_shards/README.md`'s note on this.
