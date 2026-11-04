#!/usr/bin/env bash
# QC'd HapMap3 SNP panel from per-chromosome ACAF files (hm3_panel and grm_panel).
#   PLINK2, PLINK          binaries (PLINK only needed when FREQ=1)
#   GENO_FLAG              --bfile or --pfile (ACAF ships as both)
#   GENO_PATTERN           per-chromosome prefix with a literal {chr}
#   NCHR                   number of autosomes to read (22)
#   KEEP                   one IID per line
#   EXTRACT                HapMap3 variant IDs
#   MAF HWE GENO           QC thresholds, applied within the KEEP set
#   REF_BFILE              optional: keep only variants matching it on ID+REF+ALT
#   DROP_AMBIGUOUS         1 (default) drops A/T and C/G; 0 for the GRM, which matches nothing
#   FREQ                   1 to write OUT.frq with plink 1.9 (for --read-freq)
#   OUT                    output prefix (writes OUT.bed/.bim/.fam, OUT.counts.tsv)
#   THREADS
set -euo pipefail
: "${PLINK2:?}" "${GENO_FLAG:?}" "${GENO_PATTERN:?}" "${NCHR:?}" "${KEEP:?}" "${EXTRACT:?}"
: "${MAF:?}" "${HWE:?}" "${GENO:?}" "${OUT:?}"
THREADS=${THREADS:-4}
mkdir -p "$(dirname "$OUT")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

count() { grep -vc '^#' "$1"; }
counts="$work/counts.tsv"
printf 'step\tvariants\n' > "$counts"

# 1. per chromosome: restrict to the sample set and HapMap3 biallelic SNPs
: > "$work/merge.txt"
n_in=0
for c in $(seq 1 "$NCHR"); do
  in=${GENO_PATTERN//\{chr\}/$c}
  "$PLINK2" "$GENO_FLAG" "$in" --keep "$KEEP" --extract "$EXTRACT" \
    --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all \
    --threads "$THREADS" --silent --make-pgen --out "$work/chr$c"
  n_in=$((n_in + $(count "$work/chr$c.pvar")))
  echo "$work/chr$c" >> "$work/merge.txt"
done
printf 'hm3_biallelic_snps\t%s\n' "$n_in" >> "$counts"
"$PLINK2" --pmerge-list "$work/merge.txt" pfile --threads "$THREADS" --silent \
  --make-pgen --out "$work/merged"
n_samples=$(grep -vc '^#' "$work/merged.psam")

# 2. QC within this sample set: missingness, MAF, HWE
"$PLINK2" --pfile "$work/merged" --geno "$GENO" --maf "$MAF" --hwe "$HWE" \
  --threads "$THREADS" --silent --write-snplist --out "$work/qc"
printf 'after_qc\t%s\n' "$(wc -l < "$work/qc.snplist" | tr -d " ")" >> "$counts"

# 3. drop strand-ambiguous SNPs (A/T, C/G)
awk 'BEGIN{OFS="\t"} !/^#/ {r=$4; a=$5;
     if (!((r=="A"&&a=="T")||(r=="T"&&a=="A")||(r=="C"&&a=="G")||(r=="G"&&a=="C")))
       print $3, $4, $5}' "$work/merged.pvar" | LC_ALL=C sort > "$work/unamb.tsv"
cut -f1 "$work/unamb.tsv" | LC_ALL=C sort > "$work/unamb.ids"
if [ "${DROP_AMBIGUOUS:-1}" = "1" ]; then
  LC_ALL=C sort "$work/qc.snplist" | LC_ALL=C comm -12 - "$work/unamb.ids" > "$work/keep.ids"
  printf 'unambiguous\t%s\n' "$(wc -l < "$work/keep.ids" | tr -d " ")" >> "$counts"
else
  LC_ALL=C sort "$work/qc.snplist" > "$work/keep.ids"
fi

# 4. optionally, keep only variants that match the reference on ID + REF + ALT
if [ -n "${REF_BFILE:-}" ]; then
  "$PLINK2" --bfile "$REF_BFILE" --extract "$work/keep.ids" --freq counts \
    --threads "$THREADS" --silent --out "$work/ref"
  awk 'NR>1 {print $2"\t"$3"\t"$4}' "$work/ref.acount" | LC_ALL=C sort > "$work/ref.tsv"
  LC_ALL=C comm -12 "$work/unamb.tsv" "$work/ref.tsv" | cut -f1 | LC_ALL=C sort > "$work/matched.ids"
  LC_ALL=C comm -12 "$work/keep.ids" "$work/matched.ids" > "$work/keep2.ids"
  mv "$work/keep2.ids" "$work/keep.ids"
  printf 'matched_to_reference\t%s\n' "$(wc -l < "$work/keep.ids" | tr -d " ")" >> "$counts"
fi

# 5. write the panel (plink 1 bed: what bigsnpr and plink 1.9 --make-grm-bin read)
"$PLINK2" --pfile "$work/merged" --extract "$work/keep.ids" \
  --threads "$THREADS" --silent --make-bed --out "$OUT"
printf 'final\t%s\n' "$(wc -l < "$OUT.bim" | tr -d " ")" >> "$counts"
printf 'samples\t%s\n' "$n_samples" >> "$counts"

if [ "${FREQ:-0}" = "1" ]; then
  : "${PLINK:?FREQ=1 needs PLINK (1.9)}"
  "$PLINK" --bfile "$OUT" --freq --threads "$THREADS" --silent --out "$OUT"
fi
cp "$counts" "$OUT.counts.tsv"
column -t "$OUT.counts.tsv"
