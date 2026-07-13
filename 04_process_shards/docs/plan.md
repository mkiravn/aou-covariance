# Process GRM shards

Bin each shard by relatedness and accumulate phenotype cross-products using
`GRM-pairs/grm_bin_sharded/grm_shard_tool` (`accumulate`, then `merge`),
never assembling the full dense GRM. See the submodule's own README for the
accumulator design and the row-range recovery (ported from plink's actual
`--parallel` split formula).

## Inputs this stage needs from the earlier ones

- `${OUT_PREFIX}.grm.bin.k` shards + shared `${OUT_PREFIX}.grm.id` (03)
- Normalized phenotype, `FID IID value` (02) -- looked up by ID, no
  particular row order required
- `bins.txt` (relatedness bin edges) -- copy from
  `GRM-pairs/full_grm_bin/bins.txt` or define your own
- `--nblocks`/`--seed` for the individual-block jackknife -- must be
  identical across every shard's `accumulate` call

## Sketch

```bash
for k in $(seq 1 "${N_SHARDS}"); do
    gsutil cp "${WORKSPACE_BUCKET}/data/grm_shards/${OUT_PREFIX}.grm.bin.${k}" .
    GRM-pairs/grm_bin_sharded/grm_shard_tool accumulate \
        --grm-id "${OUT_PREFIX}.grm.id" \
        --shard "${OUT_PREFIX}.grm.bin.${k}" \
        --parallel "${k}" "${N_SHARDS}" \
        --pheno phenotype_norm.txt \
        --bins bins.txt \
        --nblocks 50 --seed 1 \
        --out "shard_${k}.acc.tsv"
done

ls shard_*.acc.tsv > acc_list.txt
GRM-pairs/grm_bin_sharded/grm_shard_tool merge \
    --acc-list acc_list.txt --bins bins.txt --nblocks 50 \
    --out-prefix merged
```

Not yet run against real shards -- this is the same shape as
`GRM-pairs/grm_bin_sharded/test/test_sharded_vs_full.sh`, which *is*
verified against real 1000 Genomes data, just pointed at AoU-produced
shards and phenotypes instead of the manufactured test fixtures.

## Open items

- [ ] Decide `--nblocks` given the actual analysis cohort size.
- [ ] Wire this into a dsub job / task array rather than a local loop, once
      shard counts get large.
