// Condensed deck, ~14 slides + appendix: node slides/build_condensed.js [out.pptx] [figdir]
const path = require("path");
const OUT = process.argv[2] || path.join(__dirname, "relatedness_condensed.pptx");
const D = require("./deck_lib")(process.argv[3]);
const { pres, C, GREY, GROUPS, EX, section, content, perGroup, fig, bullets, card, cite, hd } = D;
const sz = require("image-size");
const EQ = (k) => path.join(__dirname, "eq", `${k}.png`);
const MIXED = GROUPS.includes("uniform_10k") ? "uniform_10k" : GROUPS[GROUPS.length - 1];

function eqw(s, k, x, y, w, hmax) {   // equation scaled to a width, capped in height
  const d = sz(EQ(k)), h = Math.min(hmax, w * d.height / d.width);
  s.addImage({ path: EQ(k), x, y, h, w: h * d.width / d.height, altText: k, objectName: `eq:${k}` });
  return h;
}
function label(s, t, x, y, w) {
  s.addText(t, { x, y, w, h: 0.3, fontSize: 13, bold: true, color: "444444", margin: 0, isTextBox: true });
}
function forGroup(g, cb) { D.setRun(g); cb(); D.setRun(EX); }

// ---------------------------------------------------------------- main
pres.addSection({ title: "Main" }); D.setSec("Main");
let s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Main" });
s.addText("Relatedness and phenotypic covariance across AoU ancestry groups", { placeholder: "title" });
s.addText("GRM, PC-Relate and Weir-Goudet · 8 groups · 7 traits", { placeholder: "body" });

s = content("Groups");
const rows = [["chb_chs", "EAS gate, then CHB+CHS", "11,313", "5,646"], ["yri", "AFR gate, then YRI", "97,335", "553"],
              ["lwk", "AFR gate, then LWK", "97,335", "71"], ["stu_itu_beb", "SAS gate, then STU+ITU+BEB", "5,145", "1,746"],
              ["pjl", "SAS gate, then PJL", "5,145", "3,902"], ["asw_10k", "10k of the ASW gate (90,160)", "97,335", "10,000"],
              ["uniform_10k", "random 10k of all AoU", "—", "10,000"], ["even_10k", "10k spread evenly over AoU PCs 1-4", "—", "10,000"]];
s.addTable([[hd("group"), hd("how"), hd("round 1"), hd("round 2")], ...rows],
           { x: 0.6, y: 1.4, w: 7.6, colW: [1.7, 3.5, 1.2, 1.2], fontSize: 14, color: C.text1, border: { type: "solid", pt: 0.5, color: "C9D4DF" }, rowH: 0.5, valign: "middle" });
s.addText("of 535,662 AoU participants with srWGS", { x: 0.6, y: 6.0, w: 7.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });
bullets(s, ["round 1: ellipsoid on AoU's PCs 1-5 around the 1000G superpopulation (99% of its founders)",
            "round 2: PCA on that superpopulation's 1000G founders, cohort projected (OADP), robust Mahalanobis gate on PCs 1-4 around the target population(s)",
            "mixed sets (uniform_10k, even_10k) skip the gates: they show what structure does to each estimator"], 8.6, 1.4, 4.1, 5.2, 15);

s = content("Ancestry PCA");
bullets(s, ["bigsnpr bed_autoSVD on 1000G founders, 10 PCs; clumped r² 0.2, long-range LD removed",
            "cohort projected with OADP (bed_projectPCA), not plink2 --score: with ~500 training samples simple projection shrinks PCs toward 0 (yri 553 kept vs 151, lwk 71 vs 0)",
            "covariate PCs within each group: fit on unrelated (AoU kinship > 0.1 pairs out), relatives projected back; PCs 1-10 for phenotypes, 1-4 for PC-Relate"], 0.6, 1.4, 5.0, 5.3, 15);
label(s, `${EX}: cohort on the 1000G PCs (top), round-2 gate (bottom)`, 5.9, 1.3, 6.8);
fig(s, "round2_pca/cohort_on_reference.png", 5.9, 1.65, 6.8, 2.55, "cohort on the 1000G PCs");
fig(s, "ancestry_round2/gate_pcs.png", 5.9, 4.25, 6.8, 2.55, "round-2 gate: kept vs dropped");
cite(s, "Privé et al. 2018, 2020 Bioinformatics (bigsnpr, OADP); Rousseeuw & Van Driessen 1999 (MCD)");

s = content("QC and HWE");
bullets(s, ["HapMap3 autosomal SNPs, QC'd within each group",
            "missingness ≤ 5%, MAF ≥ 1%",
            "HWE p ≥ 1e-6 in single-ancestry groups; none in mixed sets",
            "same SNP panel for all three estimators",
            "F = 1 − H_obs / 2p(1−p): ~0 in homogeneous groups, > 0 in mixed sets (Wahlund)"], 0.6, 1.4, 4.6, 5.3, 16);
fig(s, "meta/f_by_maf.png", 5.4, 1.3, 7.3, 5.5, "all groups: F by MAF bin");

s = content("Three relatedness estimators");
const cols = [["GRM", "grm", null, "reference: the group's allele frequencies", "Yang et al. 2011 AJHG (GCTA)"],
              ["PC-Relate", "pcrelate", "pcrelate_mu", "reference: each person's own allele frequencies from PCs 1-4", "Conomos et al. 2016 AJHG"],
              ["Weir-Goudet", "wg", "wg_m", "reference: the average pair of distinct individuals in the run (S*); no allele frequencies", "Weir & Goudet 2026 Genetics, Box 3"]];
cols.forEach(([t, e1, e2, ref, ref2], k) => {
  const x = 0.6 + k * 4.1, w = 3.8;
  s.addText(t, { x, y: 1.35, w, h: 0.45, fontSize: 20, bold: true, color: C.text1, margin: 0, isTextBox: true });
  let y = 1.95 + eqw(s, e1, x, 1.95, w, 0.75) + 0.25;
  if (e2) y += eqw(s, e2, x, y, w, 0.42) + 0.2;
  s.addText(ref, { x, y: Math.max(y, 3.9), w, h: 1.0, fontSize: 15, color: "333333", margin: 0, valign: "top", isTextBox: true });
  s.addText(ref2, { x, y: 5.1, w, h: 0.4, fontSize: 12, color: "888888", margin: 0, isTextBox: true });
});
eqw(s, "rel", 0.6, 5.8, 6.0, 0.55);
s.addText("relatedness = 2 × kinship for all three", { x: 6.9, y: 5.85, w: 5.8, h: 0.45, fontSize: 15, color: GREY, isTextBox: true });

s = content("Distributions");
fig(s, "meta/relatedness_distributions.png", 0.6, 1.3, 7.4, 5.6, "all groups: bulk (top), relatives (bottom, log)");
fig(s, "meta/relatedness_spread.png", 8.2, 1.3, 4.5, 5.6, "SD and share > 0.1");

s = content("Estimators against each other");
for (const [k, g] of [[0, EX], [1, MIXED]]) forGroup(g, () => {
  fig(s, "weir_goudet/estimator_comparison.png", 0.6 + k * 6.15, 1.65, 5.95, 5.0, `${g}: cells by mean PC distance`);
  label(s, g, 0.6 + k * 6.15, 1.3, 5.9);
});

s = content("Do they rank pairs the same way?");
fig(s, "meta/rank_agreement.png", 0.6, 1.3, 6.2, 5.6, "all groups: Spearman ρ");
forGroup(EX, () => fig(s, "weir_goudet/rank_comparison.png", 7.0, 1.65, 5.7, 5.2, "rank vs rank"));
label(s, `${EX}: rank vs rank, density (top) and PC distance (bottom)`, 7.0, 1.3, 5.7);

s = content(`PC distance and IBD (${EX})`);
fig(s, "weir_goudet/relatedness_vs_pc_distance.png", 0.6, 1.3, 12.1, 2.7, "relatedness by PC-distance bin");
fig(s, "ibd_vs_kinship/binned_by_ibd.png", 0.6, 4.1, 7.0, 2.8, "iLASH IBD vs each estimator");
bullets(s, ["PC distance: Euclidean, round-2 PCs 1-4", "IBD: iLASH segments, no allele frequencies; pairs without segments at 0"], 7.9, 4.2, 4.8, 2.6, 15);
cite(s, "Shemirani et al. 2021 Nat Commun (iLASH)");

s = content("Phenotypes and model");
bullets(s, ["height, BMI, systolic BP, HDL, LDL, platelets, hemoglobin; most recent measurement",
            "residualised on age, sex and 10 covariate PCs, then standardised",
            "every pair contributes y_i y_j; fit on all pairs and on pairs with relatedness < 0.1",
            "delete-block jackknife: people in 50 random blocks, drop pairs touching a block; 95% CI, t on 49 df",
            "Pearson and Spearman r²"], 0.6, 1.4, 6.8, 5.3, 16);
eqw(s, "cov", 7.8, 1.7, 4.9, 0.55);
eqw(s, "jk", 7.8, 2.7, 4.9, 0.7);

s = content(`Covariance against relatedness, to 0.1 (${EX})`);
fig(s, "kin_pheno_cov/binned_covariance_bulk.png", 0.6, 1.3, 12.1, 5.6, "bulk");

s = content(`Covariance against relatedness, to 0.6 (${EX})`);
fig(s, "kin_pheno_cov/binned_covariance_full.png", 0.6, 1.3, 12.1, 5.6, "full range, ±SE");

s = content("Fits across groups: all pairs");
fig(s, "meta/fits_all_pairs.png", 0.6, 1.3, 8.0, 5.6, "slope, intercept, Pearson and Spearman r²");
bullets(s, ["…", "…"], 8.9, 1.5, 3.8, 5.0, 16);

s = content("Fits across groups: relatedness < 0.1");
fig(s, "meta/fits_unrelated.png", 0.6, 1.3, 8.0, 5.6, "slope, intercept, Pearson and Spearman r²");
bullets(s, ["…", "…"], 8.9, 1.5, 3.8, 5.0, 16);

s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Main" });
s.addText("Takeaways", { placeholder: "title" });
s.addText("…  ·  …  ·  …", { placeholder: "body" });

// ---------------------------------------------------------------- appendix: per group, then everything else
section("Appendix", "per-group figures");
for (const g of GROUPS) forGroup(g, () => {
  let a = content(`${g}: estimators against each other`);
  fig(a, "weir_goudet/estimator_comparison.png", 0.6, 1.3, 12.1, 5.6);
  a = content(`${g}: rank vs rank`);
  fig(a, "weir_goudet/rank_comparison.png", 0.6, 1.3, 12.1, 5.6);
  a = content(`${g}: covariance against relatedness, to 0.1`);
  fig(a, "kin_pheno_cov/binned_covariance_bulk.png", 0.6, 1.3, 12.1, 5.6);
  a = content(`${g}: covariance against relatedness, to 0.6`);
  fig(a, "kin_pheno_cov/binned_covariance_full.png", 0.6, 1.3, 12.1, 5.6);
  a = content(`${g}: slopes, intercepts, r²`);
  fig(a, "kin_pheno_cov/fit_coefficients.png", 0.6, 1.3, 6.0, 5.6, "top half", { t: 0, b: 0.5 });
  fig(a, "kin_pheno_cov/fit_coefficients.png", 6.7, 1.3, 6.0, 5.6, "bottom half", { t: 0.5, b: 0 });
});
D.supplementary();
D.write(OUT);
