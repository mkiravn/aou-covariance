#!/usr/bin/env bash
# Time readers for the phased VCF on one region, one core each; extrapolates core-hours for the genome.
#   bash scripts/dev/phased_bench.sh [run] [chr] [start_mb] [length_mb]     e.g. chb_chs 20 40 1
set -euo pipefail
RUN=${1:-chb_chs}; CHR=${2:-20}; MB=${3:-40}; LEN=${4:-1}
PY=$(command -v python || command -v python3)
cd "$(dirname "$0")/../.."
eval "$($PY - "$RUN" "$CHR" <<'EOF'
import os, sys
from aoucov import load_run
c = load_run(sys.argv[1])
print(f'VCF={os.path.expanduser(c.input("phased_vcf", chr=sys.argv[2]))}')
print(f'HM3={os.path.expanduser(c.input("hm3_ids"))}')
print(f'PLINK2={os.path.expanduser(c.tool("plink2"))}')
print(f'W={os.path.expanduser(c.params.get("scratch", "~"))}/phased_bench')
EOF
)"
mkdir -p "$W"; cd "$W"
contig=$(bcftools view -h "$VCF" | awk -F'[=,>]' -v c="$CHR" '/^##contig=<ID=/ && ($3 == "chr"c || $3 == c) && !k {k = $3} END {print k}')
REG="$contig:$((MB * 1000000))-$(((MB + LEN) * 1000000 - 1))"
t() { local s=$(date +%s.%N); "$@"; echo "$(echo "$(date +%s.%N) - $s" | bc)"; }
res() { printf '%-28s %8.0f s  %6s variants  -> %6.1f core-h for 2900 Mb\n' "$1" "$2" "$3" "$(echo "$2 * 2900 / $LEN / 3600" | bc -l)"; }

if [[ ! -s slice.vcf.gz ]]; then
  echo "fetching $REG $(date +%T)"
  s=$(t bash -c "tabix -h '$VCF' '$REG' | bgzip -l 1 -@ 16 > slice.vcf.gz && tabix -p vcf slice.vcf.gz")
  echo "fetch: ${s%.*} s, $(du -h slice.vcf.gz | cut -f1), $(zcat slice.vcf.gz | grep -vc '^#') records"
fi
awk -F: -v OFS='\t' -v c="$CHR" -v k="$contig" '$1 == c {print k, $2}' "$HM3" | sort -k2,2n -u > targets.tsv
cut -f2 targets.tsv > pos.txt
P2="--vcf-half-call m --output-chr 26 --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 1000 --extract $HM3 --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --threads 1 --silent"
nvar() { grep -vc '^#' "$1" || true; }

s=$(t bash -c "bcftools view -T targets.tsv -f PASS -Ou slice.vcf.gz | bcftools norm -m- -Ou | bcftools view -m2 -M2 -v snps -Ob -o a.bcf && $PLINK2 --bcf a.bcf $P2 --make-pgen --out a")
res "A bcftools (current)" "$s" "$(nvar a.pvar)"

s=$(t bash -c "$PLINK2 --vcf slice.vcf.gz --var-filter $P2 --make-pgen --out b")
res "B plink2 --vcf" "$s" "$(nvar b.pvar)"

cat > keep.py <<'EOF'
import sys
want = set(open(sys.argv[1]).read().split())
out = sys.stdout.buffer
for line in sys.stdin.buffer:
    if line[:1] == b"#":
        out.write(line); continue
    a = line.index(b"\t"); b = line.index(b"\t", a + 1)
    if line[a + 1:b].decode() in want:
        out.write(line)
EOF
s=$(t bash -c "bgzip -dc -@ 1 slice.vcf.gz | $PY keep.py pos.txt > c.vcf && $PLINK2 --vcf c.vcf --var-filter $P2 --make-pgen --out c")
res "C skip lines, then plink2" "$s" "$(nvar c.pvar)"

if $PY -c "import hail" 2>/dev/null; then
  cat > hl_bench.py <<'EOF'
import sys, hail as hl
hl.init(master="local[1]", quiet=True, log="/dev/null")
pos = hl.literal({int(p) for p in open(sys.argv[1]).read().split()})
mt = hl.import_vcf("slice.vcf.gz", reference_genome="GRCh38", force_bgz=True, array_elements_required=False)
mt = mt.filter_rows(pos.contains(mt.locus.position) & (hl.len(mt.filters) == 0) & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (hl.len(mt.alleles) == 2))
mt.select_entries("GT").write("d.mt", overwrite=True)
print(hl.read_matrix_table("d.mt").count_rows())
EOF
  s0=$(date +%s.%N); n=$($PY hl_bench.py pos.txt 2>/dev/null | tail -1); s=$(echo "$(date +%s.%N) - $s0" | bc)
  res "D hail local[1]" "$s" "$n"
else
  echo "D hail: not installed on this VM"
fi
for x in a b c; do [[ -s $x.pgen ]] && "$PLINK2" --pfile $x --export vcf --silent --out $x.chk && grep -v '^##' $x.chk.vcf | md5sum | sed "s/ .*/  $x genotypes+phase/"; done
echo "slice kept in $W (delete when done)"
