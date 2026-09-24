

# Also get the population panel file if you haven't already
wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

awk 'NR>1 && $3=="EUR" {print $1}' \
  integrated_call_samples_v3.20130502.ALL.panel \
  > eur_samples.txt

awk 'NR>1 && ($2=="CEU" || $2=="GBR") {print $1}' \
  integrated_call_samples_v3.20130502.ALL.panel \
  > ceu_gbr_samples.txt

wc -l eur_samples.txt ceu_gbr_samples.txt
# Expect: ~503 EUR total, ~192 CEU+GBR

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz


./plink2 \
  --vcf ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
  --keep eur_samples.txt \
  --maf 0.01 \
  --max-alleles 2 \
  --rm-dup exclude-all \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 1000 \
  --make-bed \
  --out 1kg_eur_chr22

# LD pruning
./plink2 \
  --bfile 1kg_eur_chr22 \
  --indep-pairwise 200kb 1 0.1 \
  --out eur_pruned

# Check how many variants survive
wc -l eur_pruned.prune.in

./plink2 \
  --bfile 1kg_eur_chr22 \
  --extract eur_pruned.prune.in \
  --freq counts \
  --pca allele-wts 20 \
  --out eur_pca

ls -lh eur_pca.*

python ellipsoid_filter.py