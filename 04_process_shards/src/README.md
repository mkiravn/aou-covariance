# grm_class_tool

Bins phenotype cross-products from a row-chunked sharded GRM, keeping
designated pair classes in their **own bin sets** rather than pooling them.

Parent-offspring and full-sib pairs have the same expected additive
relatedness (`a_ij ≈ 0.5`), so the GRM cannot separate them and they land in
the same bin. They differ in dominance (`d_ij` = 0 for PO, 0.25 for FS) and in
shared environment, so pooling them both biases the FS bin and discards a real
contrast. Given an external classification — e.g. from IBS0, via
`04_process_shards/notebooks/remote/00_relatedness_screen_cells.md` — this
routes each listed pair into a per-class copy of the bin array.

Nothing is excluded. One pass yields the mean cross-product for each class
separately, plus a synthetic `pooled` class reproducing exactly what an
unclassed run would have produced.

Standalone by design: `GRM-pairs/grm_bin_sharded/grm_shard_tool` is a git
submodule and is left untouched, remaining the reference implementation for
the unclassed case. The row-range recovery and accumulator maths here are
ported from it unchanged so classed and unclassed runs stay comparable.

## Build

```
make          # build
make test     # build, then check against a Python reference
make clean
```

## Use

```
grm_class_tool accumulate \
  --grm-id   shards/grm_shard_1_of_16.grm.id \
  --shard    shards/grm_shard_7_of_16.grm.bin.7 \
  --parallel 7 16 \
  --pheno    height__invnorm__base_pcs_aligned.pheno \
  --bins     GRM-pairs/full_grm_bin/bins.txt \
  --pair-classes deg1_classified.tsv \
  --nblocks 50 --seed 1 \
  --out      acc.7.tsv

grm_class_tool merge \
  --acc-list acc_list.txt \
  --bins     GRM-pairs/full_grm_bin/bins.txt \
  --nblocks 50 \
  --out-prefix merged
```

`--pair-classes` is optional; without it everything goes to class `other` and
the result matches `grm_shard_tool` exactly.

### Pair-class file

TSV with a header. Defaults to columns `IID1`, `IID2`, `cls`, overridable with
`--id-col1` / `--id-col2` / `--class-col`. Extra columns are ignored, so the
screen's `deg1_classified.tsv` can be passed as-is.

```
IID1        IID2        cls
<person_id> <person_id> PO
<person_id> <person_id> FS
```

Class names are taken from the file; any pair not listed is `other`.

**Keys on IID alone.** plink writes `FID = 0` for every row of a `.grm.id`
built without pedigree, so `(FID, IID)` carries no extra information and would
force every downstream file to reproduce plink's convention. IID uniqueness is
checked at load and errors if violated.

Pair orientation does not matter — pairs are normalised to `(max, min)` row
index internally. Reversed duplicates are counted and ignored.

## Output

`<prefix>.full.tsv` — one row per `(class, bin)`:

```
class  bin_index  bin_left  bin_right  bin_midpoint  full_sum  full_sum_sq
full_n  full_mean  full_sd  full_se  jk_mean  jk_var  jk_se
```

`<prefix>.jk.tsv` — per `(class, block, bin)` delete-block values.

The `pooled` class is synthesised at merge time as the sum over real classes.

## Progress

`accumulate` reports progress against the shard's entry count, which is linear
in work (row `i` holds `i` entries, so row count is not):

```
[shard 16/16]  30%  1m12s elapsed, 2m48s left, 21M entries/s
[shard 16/16] done in 4m01s (21M entries/s)
```

On a terminal it overwrites in place every 2%; when redirected to a log it
prints one line every 10%, so a 288-combo batch doesn't emit thousands of bar
fragments. The counter is updated once per row, never per entry, so the inner
loop is untouched.

The entries/s rate is the useful number for extrapolating a full run: every
shard holds the same entry count (plink balances `--parallel` on pair count),
so one shard's time times the shard count gives the per-combo total.

## Validation

`accumulate` reports to stderr, per shard:

```
[INFO] pair classes: 3412 read, 3401 mapped, 11 unmapped, 0 duplicate
[INFO] binned pairs by class: other=1538221093 PO=214 FS=602 DUP=9
```

**Check the mapped and per-class counts.** A format or ID mismatch maps zero
pairs and produces a clean run in which every classed pair silently falls into
`other` — indistinguishable from "there were no classed pairs" unless you look.
Summing the per-class hits across shards should reconcile against the pair
counts in the classification file, less pairs dropped for a missing phenotype
or for falling outside every bin.

## Cost

The per-pair check is a single integer compare against a per-row sorted list of
classed partners, with a monotonic cursor — it short-circuits immediately on
rows with no classed partners, which is nearly all of them. At N ≈ 222K there
are ~2.46e10 pairs and only a few thousand classed ones, so the overhead is not
measurable against the existing bin lookup.
