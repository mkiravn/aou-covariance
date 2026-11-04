#!/usr/bin/env bash
# One GRM shard. plink 1.9: grm_class_tool's row ranges assume its --parallel split.
#   PLINK      plink 1.9 binary
#   BFILE      panel prefix from the grm_panel stage (has BFILE.frq beside it)
#   K N        this shard, of N
#   OUT        output prefix; plink writes OUT.grm.bin.K
#   THREADS
set -euo pipefail
: "${PLINK:?}" "${BFILE:?}" "${K:?}" "${N:?}" "${OUT:?}"
mkdir -p "$(dirname "$OUT")"
"$PLINK" --bfile "$BFILE" --read-freq "$BFILE.frq" \
  --make-grm-bin --parallel "$K" "$N" \
  --threads "${THREADS:-4}" --silent --out "$OUT"
ls -l "$OUT".grm.*
