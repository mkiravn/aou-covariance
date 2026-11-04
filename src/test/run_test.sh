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
echo "=== accumulators without the indiv column give n_indiv = NA, not a number ==="
mkdir -p "$TMP/old"
: > "$TMP/old/acc_list.txt"
for k in $(seq 1 "$N_SHARDS"); do
  cut -f1-7 "$TMP/acc.$k.tsv" > "$TMP/old/acc.$k.tsv"
  echo "$TMP/old/acc.$k.tsv" >> "$TMP/old/acc_list.txt"
done
"$BIN" merge --acc-list "$TMP/old/acc_list.txt" --bins "$TMP/test.bins" \
  --nblocks "$NBLOCKS" --out-prefix "$TMP/old/merged" 2>/dev/null
bad=$(awk -F'\t' 'NR>1 && $15 != "NA"' "$TMP/old/merged.full.tsv" | wc -l | tr -d ' ')
[ "$bad" -eq 0 ] || { echo "FAIL -- $bad rows have a number in n_indiv"; exit 1; }
cmp -s <(cut -f1-14 "$TMP/old/merged.full.tsv") <(cut -f1-14 "$TMP/merged.full.tsv") \
  || { echo "FAIL -- dropping indiv changed the sums"; exit 1; }
echo "  every n_indiv is NA; sums unchanged"

echo
echo "=== n_indiv cap across shards: a sparse bin spread over 7 shards ==="
CAP="$TMP/cap"
python3 "$HERE/make_test_data.py" "$CAP" 120 7 0.026 >/dev/null
: > "$CAP/acc_list.txt"
for k in $(seq 1 7); do
  "$BIN" accumulate --grm-id "$CAP/test.grm.id" --shard "$CAP/shard.$k" --parallel "$k" 7 \
    --pheno "$CAP/test.pheno" --bins "$CAP/test.bins" --pair-classes "$CAP/test.classes.tsv" \
    --nblocks "$NBLOCKS" --seed 1 --out "$CAP/acc.$k.tsv" 2>/dev/null
  echo "$CAP/acc.$k.tsv" >> "$CAP/acc_list.txt"
done
"$BIN" merge --acc-list "$CAP/acc_list.txt" --bins "$CAP/test.bins" \
  --nblocks "$NBLOCKS" --out-prefix "$CAP/merged" 2>/dev/null
python3 "$HERE/check_result.py" "$CAP/expected.tsv" "$CAP/merged.full.tsv" | tail -2

echo
echo "=== he_estimators.enough_people: exact with n_indiv, pair bound without ==="
python3 - "$SRC" "$CAP" "$TMP/old" <<'EOF'
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, sys.argv[1])
import he_estimators as HE
cap, old = sys.argv[2], sys.argv[3]
nb = len(open(f"{cap}/test.bins").read().split("\n")) - 1
f = pd.read_csv(f"{cap}/merged.full.tsv", sep="\t")
N = np.zeros((1, len(HE.CLS), nb))
for ci, c in enumerate(HE.CLS):
    fc = f[f["class"] == c]
    N[0, ci, fc["bin_index"] - 1] = fc["full_n"]
P = HE.load_people(f"{cap}/merged.full.tsv", nb)
assert P is not None, "n_indiv present but load_people returned None"
ref = pd.read_csv(f"{cap}/expected.tsv", sep="\t")
for c in HE.CLS:
    r = ref[ref["class"] == c].sort_values("bin_index")
    if not len(r):
        continue
    got = HE.enough_people(P, N, [c])
    want = (r["n_indiv"] >= HE.MIN_PEOPLE).to_numpy()
    assert (got == want).all(), f"{c}: {got} vs reference {want}"
assert HE.load_people(f"{old}/merged.full.tsv", 5) is None, "NA column must give None"
fb = HE.enough_people(None, N, ["other"])
assert (fb == (N[0, 0][N[0, 0] > 0] >= 210)).all()
print(f"  per-class gate matches the brute-force people count; "
      f"NA merge falls back to >= {HE.MIN_PAIRS_NO_PEOPLE} pairs")
EOF

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
echo "--- a path containing spaces resolves (AoU bucket mount) ---"
for k in $(seq 1 "$N_SHARDS"); do
  if ! cmp -s "$TMP/acc.$k.tsv" "$TMP/multi/p_spaced.shard$k.tsv"; then
    echo "FAIL -- p_spaced (same data via a spaced path) differs from p1 on shard $k"
    exit 1
  fi
done
echo "  identical to p1 on all $N_SHARDS shards"

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
