# Deck text (generated from the .pptx) — edit and save, then ask for a rebuild

<!-- Slide number: 1 -->
Relatedness and phenotypic covariance across AoU ancestry groups
GRM, PC-Relate and Weir-Goudet · 8 groups · 7 traits

<!-- Slide number: 2 -->
What I'm doing
1.  ancestry filtering  — two rounds of PCA against 1000G
2.  SNP panel  — HapMap3, QC'd within each group
3.  relatedness  — GRM, PC-Relate, Weir-Goudet
4.  phenotypes  — 7 traits, residualised
5.  covariance  — y_i y_j against relatedness
Question: how well does each relatedness measure track phenotypic covariance, in unrelated pairs and in relatives?

<!-- Slide number: 3 -->
Samples
ancestry filtering, groups, SNP panel

<!-- Slide number: 4 -->
Round 1: AoU's PCs
ancestry_round1/gate_pcs.png
gate on AoU PCs
AoU's premade 16 PCs (ancestry predictions), 1000G projected on the same PCs
gate: ellipsoid on PCs 1-5 around the 1000G superpopulation (e.g. EAS)
radius: 99% quantile of the reference founders' distances
keeps the broad continental group, still structured inside

<!-- Slide number: 5 -->
Round 2: PCA on 1000G, cohort projected
round2_pca/cohort_on_reference.png
cohort on the 1000G PCs
PCA (bigsnpr bed_autoSVD) on 1000G founders of the superpopulation, 10 PCs
clumped at r² 0.2, long-range LD regions removed
round-1 set projected with OADP (bed_projectPCA)
changed from plink2 --score: with ~500 training samples, simple projection shrinks PCs toward 0 past PC4
gate: robust (MCD) Mahalanobis distance on PCs 1-4 around the target 1000G population(s), radius √χ²₄(0.999)
ancestry_round2/gate_pcs.png
kept vs dropped
Privé et al. 2018 Bioinformatics (bigsnpr); Privé et al. 2020 Bioinformatics (OADP projection); Rousseeuw & Van Driessen 1999 Technometrics (MCD)

<!-- Slide number: 6 -->
Round 2 PCA checks
round2_pca/scree.png
scree
round2_pca/simple_vs_oadp.png
simple vs OADP projection
round2_pca/reference_scores.png
1000G reference
round2_pca/loadings.png
|loadings| along the genome: no peaks

<!-- Slide number: 7 -->
How many people each gate keeps
ancestry_round1/kept_vs_coverage.png
round 1: kept vs gate coverage; plateau = cluster
pop_gates/kept_by_population.png
round 2: AoU people inside each 1000G population's gate
pop_gates/overlap_EAS.png
gate overlap (Jaccard), one per superpopulation

<!-- Slide number: 8 -->
The groups

| group | who | N | SNPs |
| --- | --- | --- | --- |
| chb_chs | East Asian, near CHB/CHS | … | … |
| yri | West African, near YRI | … | … |
| lwk | East African, near LWK | … | … |
| stu_itu_beb | South Asian, near STU/ITU/BEB | … | … |
| pjl | South Asian, near PJL | … | … |
| asw_10k | African American, admixed (10k sample) | … | … |
| uniform_10k | all ancestries, uniform sample (10k) | … | … |
| even_10k | all ancestries, even across groups (10k) | … | … |

mixed sets
uniform_10k and even_10k mix ancestries on purpose: they show what structure does to each estimator

<!-- Slide number: 9 -->
QC
samples
AoU srWGS callset (AoU's own sample QC); ancestry gates in rounds 1 and 2; relatives kept throughout
for PCA only
long-range LD regions removed; clumped (r² 0.2 on 1000G, r² 0.1 within groups); autoSVD drops SNPs with loading peaks
variants (GRM panel)
HapMap3 autosomal SNPs; missingness ≤ 5%; MAF ≥ 1%; HWE p ≥ 1e-6 in single-ancestry groups, none in mixed sets; QC'd within each group
grm/panel_qc.png
SNPs per group after QC

<!-- Slide number: 10 -->
Covariate PCs, within each group
covariate_pcs/covariate_pcs.png
scree, fit vs projected
bed_autoSVD on the group's HapMap3 common SNPs, clumped at r² 0.1, long-range LD removed
fit on unrelated people: both members of every pair with AoU kinship > 0.1 left out
relatives projected back (bed_projectSelfPCA): nobody dropped
PCs 1-10: phenotype residualisation; PCs 1-4: PC-Relate
covariate_pcs/covariate_pc_pairs.png
PC pairs, density
Privé et al. 2018 Bioinformatics (bigsnpr, autoSVD); Privé et al. 2020 Bioinformatics (projection)

<!-- Slide number: 11 -->
Genotypes against HWE
![F](eq:F)
allele_variance/variance_by_maf.png
all groups: F by MAF bin (meta page)
F > 0: fewer heterozygotes than HWE (structure, inbreeding)
expect ~0 in homogeneous groups, > 0 in mixed sets

<!-- Slide number: 12 -->
Relatedness
three estimators, one scale

<!-- Slide number: 13 -->
GRM
![grm](eq:grm)
reference
the group's allele frequencies, same for everyone
x = 0/1/2 allele count, p = allele frequency in the group
relative to the group's average: off-diagonal mean ≈ 0
structure inside the group shows up as relatedness
Yang et al. 2010 Nat Genet; Yang et al. 2011 AJHG (GCTA)

<!-- Slide number: 14 -->
PC-Relate
![pcrelate](eq:pcrelate)
reference
individual-specific allele frequencies from PCs
![pcrelate_mu](eq:pcrelate_mu)
μ: each person's own allele frequency, from PCs 1-4
removes ancestry-level sharing before measuring kinship
relatedness = 2φ
Conomos, Reiner, Weir & Thornton 2016 AJHG

<!-- Slide number: 15 -->
Weir-Goudet
![wg](eq:wg)
reference
the mean pair in the sample
![wg_m](eq:wg_m)
M: allele matching, no allele frequencies needed
relative to the average pair in the sample
relatedness = 2β (SNPRelate snpgdsIndivBeta)
Weir & Goudet 2017 Genetics; Zheng et al. 2012 Bioinformatics (SNPRelate)

<!-- Slide number: 16 -->
All on one scale
![rel](eq:rel)
GRM
group allele frequencies
PC-Relate
per-person allele frequencies from PCs
Weir-Goudet
mean pair in the sample

<!-- Slide number: 17 -->
Distributions
meta/relatedness_distributions.png
all groups: bulk and relatives
meta/relatedness_spread.png
SD, share > 0.1

<!-- Slide number: 18 -->
Estimators against each other
weir_goudet/estimator_comparison.png
chb_chs: hexes by PC distance, black = bin means

<!-- Slide number: 19 -->
Do they rank pairs the same way?
weir_goudet/rank_comparison.png
chb_chs: rank vs rank, density and PC distance
meta/rank_agreement.png
all groups: Spearman ρ

<!-- Slide number: 20 -->
Relatedness and PC distance
weir_goudet/relatedness_vs_pc_distance.png
relatedness by PC-distance bin
weir_goudet/pc_distance_by_relatedness.png
PC distance by relatedness bin
PC distance: Euclidean, round-2 PCs 1-4

<!-- Slide number: 21 -->
Against IBD (chb_chs)
ibd_vs_kinship/binned_by_ibd.png
iLASH IBD relatedness vs each estimator
IBD: shared segments from iLASH, no allele frequencies
pairs without segments at 0
Shemirani et al. 2021 Nat Commun (iLASH)

<!-- Slide number: 22 -->
Phenotypic covariance
y_i y_j against relatedness

<!-- Slide number: 23 -->
Phenotypes
height, BMI, systolic BP, HDL, LDL, platelets, hemoglobin
most recent measurement per person
residualised on age, sex and 10 PCs, then standardised
every pair contributes y_i y_j
two scopes
all pairs, and pairs with relatedness < 0.1 (the 'unrelated' slope)
![cov](eq:cov)

<!-- Slide number: 24 -->
Covariance against relatedness
kin_pheno_cov/binned_covariance_bulk.png
bulk: relatedness to 0.1
kin_pheno_cov/binned_covariance_full.png
full: relatedness to 0.6, ±SE

<!-- Slide number: 25 -->
Slopes, intercepts, r²
kin_pheno_cov/fit_coefficients.png
chb_chs: per trait and estimator
![jk](eq:jk)
delete-block jackknife: people in 50 random blocks
drop every pair touching a block
95% CI, t on 49 df
Pearson and Spearman r²

<!-- Slide number: 26 -->
Across groups
meta/fits_all_pairs.png
all pairs; ◆ inverse-variance pooled
meta/fits_unrelated.png
relatedness < 0.1

<!-- Slide number: 27 -->
Takeaways
…  ·  …  ·  …
