# GRM shard computation

Compute the GRM in row-chunked shards via plink's native `--parallel k n`,
one dsub job per shard (or a dsub task array — see
[asgarilab/SDoH-GeneticRisk-Biobank-MCA](https://github.com/asgarilab/SDoH-GeneticRisk-Biobank-MCA)
for the `aou_dsub` + `--tasks <param.tsv>` pattern this follows), storing
shards in the workspace bucket rather than downloading/reassembling the
full dense matrix locally.

## Why shards, not one job

See `GRM-pairs/grm_bin_sharded/README.md` for the full design rationale --
short version: `plink --make-grm-bin --parallel k n` splits the expensive
part (the O(N^2) SNP-by-SNP computation) across `n` independent jobs, each
producing a `.grm.bin.k` file covering a contiguous row range. Because that
row range balances by off-diagonal pair count, `GRM-pairs/grm_bin_sharded`'s
`accumulate` step recovers which rows a shard covers directly from plink's
own split formula (ported from its source), verified against the shard
file's actual size.

## Job pattern (sketch -- not yet run against real AoU data)

```bash
plink --bfile "${BFILE}" \
    --make-grm-bin --parallel "${K}" "${N_SHARDS}" \
    --nonfounders \
    --out "${OUT}"
```

`${OUT}.grm.bin.${K}` and `${OUT}.grm.id` (unsuffixed, identical across
shards) both go to the workspace bucket, e.g.
`${WORKSPACE_BUCKET}/data/grm_shards/${OUT}.grm.bin.${K}`.

## Open items

- [ ] Confirm the actual `--bfile` source: AoU's plink-converted array data,
      or an export from the WGS callset restricted to common variants.
- [ ] Decide `N_SHARDS` based on cohort size after ancestry filtering (01).
- [ ] Wire up `aou_dsub`/task-array submission once we're testing against
      a real workspace.
