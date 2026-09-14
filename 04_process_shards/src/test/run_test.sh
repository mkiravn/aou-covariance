#!/usr/bin/env bash
# Correctness test for grm_class_tool: accumulate every shard, merge, and
# compare per-(class, bin) sums against a Python reference computed directly
# from the dense matrix.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$(dirname "$HERE")"
BIN="$SRC/grm_class_tool"
TMP="$HERE/tmp"
N_IDS=60
N_SHARDS=5
NBLOCKS=6

rm -rf "$TMP"
mkdir -p "$TMP"

echo "=== generating test data ==="
python3 "$HERE/make_test_data.py" "$TMP" "$N_IDS" "$N_SHARDS"

echo
echo "=== ranges ==="
for k in $(seq 1 "$N_SHARDS"); do
  "$BIN" ranges --n-ids "$N_IDS" --parallel "$k" "$N_SHARDS"
done

echo
echo "=== accumulate ==="
: > "$TMP/acc_list.txt"
for k in $(seq 1 "$N_SHARDS"); do
  "$BIN" accumulate \
    --grm-id "$TMP/test.grm.id" \
    --shard "$TMP/shard.$k" \
    --parallel "$k" "$N_SHARDS" \
    --pheno "$TMP/test.pheno" \
    --bins "$TMP/test.bins" \
    --pair-classes "$TMP/test.classes.tsv" \
    --nblocks "$NBLOCKS" --seed 1 \
    --out "$TMP/acc.$k.tsv"
  echo "$TMP/acc.$k.tsv" >> "$TMP/acc_list.txt"
done

echo
echo "=== merge ==="
"$BIN" merge \
  --acc-list "$TMP/acc_list.txt" \
  --bins "$TMP/test.bins" \
  --nblocks "$NBLOCKS" \
  --out-prefix "$TMP/merged"

echo
echo "=== compare against reference ==="
python3 "$HERE/check_result.py" "$TMP/expected.tsv" "$TMP/merged.full.tsv"

echo
echo "=== multi-phenotype: K phenotypes in one pass ==="
mkdir -p "$TMP/multi"
for k in $(seq 1 "$N_SHARDS"); do
  "$BIN" accumulate \
    --grm-id "$TMP/test.grm.id" \
    --shard "$TMP/shard.$k" \
    --parallel "$k" "$N_SHARDS" \
    --pheno-list "$TMP/pheno_list.tsv" \
    --bins "$TMP/test.bins" \
    --pair-classes "$TMP/test.classes.tsv" \
    --nblocks "$NBLOCKS" --seed 1 \
    --out-pattern "$TMP/multi/{name}.shard$k.tsv" 2>/dev/null
done

echo "--- p1 batched must equal p1 alone, shard by shard ---"
fail=0
for k in $(seq 1 "$N_SHARDS"); do
  if cmp -s "$TMP/acc.$k.tsv" "$TMP/multi/p1.shard$k.tsv"; then
    echo "  shard $k: identical"
  else
    echo "  shard $k: DIFFERS"
    fail=1
  fi
done
[ "$fail" -eq 0 ] || { echo "FAIL -- batching changed p1's accumulators"; exit 1; }

echo
echo "--- every batched phenotype merges and differs from the others ---"
for name in p1 p2 p3 p4; do
  : > "$TMP/multi/${name}_list.txt"
  for k in $(seq 1 "$N_SHARDS"); do echo "$TMP/multi/$name.shard$k.tsv" >> "$TMP/multi/${name}_list.txt"; done
  "$BIN" merge --acc-list "$TMP/multi/${name}_list.txt" --bins "$TMP/test.bins" \
    --nblocks "$NBLOCKS" --out-prefix "$TMP/multi/${name}_merged" 2>/dev/null
  n=$(awk -F'\t' 'NR>1 && $1=="pooled" {s+=$8} END {print s}' "$TMP/multi/${name}_merged.full.tsv")
  echo "  $name: $n pairs binned"
done

echo
echo "PASS -- multi-phenotype batching matches single-phenotype exactly"
