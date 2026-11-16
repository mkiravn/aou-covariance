// Deck scaffold: node slides/build_deck.js [out.pptx] [figdir]
// Figure slots are named "fig:<file>"; with figdir, any <figdir>/<file> found is placed in its slot.
const fs = require("fs"), path = require("path");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require(process.env.PPTX_SKILL + "/scripts/apply_theme.js");

const OUT = process.argv[2] || path.join(__dirname, "relatedness_scaffold.pptx");
const FIGDIR = process.argv[3];
const EQ = (k) => path.join(__dirname, "eq", `${k}.png`);
const THEME = {
  name: "aoucov-plain", headFontFace: "Calibri", bodyFontFace: "Calibri",
  colors: { dk1: "1A1A1A", lt1: "FFFFFF", dk2: "444444", lt2: "F2F2F2", accent1: "444444", accent2: "777777",
            accent3: "999999", accent4: "BBBBBB", accent5: "DDDDDD", accent6: "EEEEEE", hlink: "1F4E79", folHlink: "444444" },
};
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";   // 13.33 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Relatedness and phenotypic covariance";
const C = pres.SchemeColor;
const GREY = "666666";

pres.defineSlideMaster({ title: "TITLE", background: { color: "FFFFFF" }, objects: [
  { placeholder: { options: { name: "title", type: "title", align: "left", x: 0.8, y: 2.3, w: 11.5, h: 1.6, fontSize: 36, bold: true, color: C.text1, valign: "bottom" }, text: "" } },
  { placeholder: { options: { name: "body", type: "body", align: "left", x: 0.8, y: 4.1, w: 11.5, h: 1.0, fontSize: 18, color: C.text2 }, text: "" } },
] });
pres.defineSlideMaster({ title: "SECTION", background: { color: "FFFFFF" }, objects: [
  { placeholder: { options: { name: "title", type: "title", align: "left", x: 0.8, y: 2.8, w: 11.5, h: 1.0, fontSize: 36, bold: true, color: C.text1 }, text: "" } },
  { placeholder: { options: { name: "body", type: "body", align: "left", x: 0.8, y: 3.8, w: 11.5, h: 0.8, fontSize: 18, color: C.text2 }, text: "" } },
] });
pres.defineSlideMaster({ title: "CONTENT", background: { color: "FFFFFF" }, slideNumber: { x: 12.4, y: 7.0, fontSize: 10, color: "999999" }, objects: [
  { placeholder: { options: { name: "title", type: "title", align: "left", x: 0.6, y: 0.35, w: 12.1, h: 0.8, fontSize: 28, bold: true, color: C.text1, valign: "top" }, text: "" } },
] });

let sec = "";
const section = (t, sub) => {
  sec = t; pres.addSection({ title: t });
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: t });
  s.addText(t, { placeholder: "title" }); if (sub) s.addText(sub, { placeholder: "body" });
};
const content = (t) => { const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: sec }); s.addText(t, { placeholder: "title" }); return s; };

function fig(s, file, x, y, w, h, note) {
  const f = FIGDIR && path.join(FIGDIR, file);
  if (f && fs.existsSync(f)) {
    const sz = require("image-size")(f), r = Math.min(w / sz.width, h / sz.height);
    const iw = sz.width * r, ih = sz.height * r;
    s.addImage({ path: f, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih, objectName: `fig:${file}`, altText: file });
    return;
  }
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: "FFFFFF" }, line: { color: "BBBBBB", width: 1, dashType: "dash" }, objectName: `fig:${file}` });
  s.addText([{ text: file, options: { fontSize: 14, color: "555555", breakLine: !!note } },
             ...(note ? [{ text: note, options: { fontSize: 12, color: "888888" } }] : [])],
            { x, y, w, h, align: "center", valign: "middle", isTextBox: true });
}
function cite(s, t) {
  s.addText(t, { x: 0.6, y: 6.85, w: 11.6, h: 0.4, fontSize: 11, color: "888888", isTextBox: true });
}
function eq(s, k, x, y, h) {
  const sz = require("image-size")(EQ(k));
  s.addImage({ path: EQ(k), x, y, h, w: h * sz.width / sz.height, altText: k, objectName: `eq:${k}` });
}
function bullets(s, items, x, y, w, h, size = 18) {
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1 } })),
            { x, y, w, h, fontSize: size, color: C.text1, paraSpaceAfter: 10, valign: "top", isTextBox: true });
}
function card(s, x, y, w, h, head, body) {   // plain: bold heading, text below
  s.addText([{ text: head, options: { bold: true, fontSize: 18, breakLine: true } }, { text: body, options: { fontSize: 16, color: "333333" } }],
            { x, y, w, h, color: C.text1, valign: "top", paraSpaceAfter: 6, isTextBox: true });
}
// ---------------------------------------------------------------- title
pres.addSection({ title: "Intro" }); sec = "Intro";
let s = pres.addSlide({ masterName: "TITLE", sectionTitle: sec });
s.addText("Relatedness and phenotypic covariance across AoU ancestry groups", { placeholder: "title" });
s.addText("GRM, PC-Relate and Weir-Goudet · 8 groups · 7 traits", { placeholder: "body" });

s = content("What I'm doing");
const steps = [["ancestry filtering", "two rounds of PCA against 1000G"], ["SNP panel", "HapMap3, QC'd within each group"],
               ["relatedness", "GRM, PC-Relate, Weir-Goudet"], ["phenotypes", "7 traits, residualised"], ["covariance", "y_i y_j against relatedness"]];
s.addText(steps.map(([h, b], k) => ({ text: `${k + 1}.  ${h}  `, options: { bold: true, breakLine: false } })).flatMap((r, k) =>
  [r, { text: "— " + steps[k][1], options: { color: GREY, breakLine: k < steps.length - 1 } }]),
  { x: 0.6, y: 1.5, w: 12.1, h: 3.4, fontSize: 20, color: C.text1, paraSpaceAfter: 12, valign: "top", isTextBox: true });
s.addText("Question: how well does each relatedness measure track phenotypic covariance, in unrelated pairs and in relatives?",
          { x: 0.6, y: 5.4, w: 12.1, h: 0.8, fontSize: 18, color: C.text1, isTextBox: true });

// ---------------------------------------------------------------- samples
section("Samples", "ancestry filtering, groups, SNP panel");
s = content("Round 1: AoU's PCs");
bullets(s, ["AoU's premade 16 PCs (ancestry predictions), 1000G projected on the same PCs",
            "gate: ellipsoid on PCs 1-5 around the 1000G superpopulation (e.g. EAS)",
            "radius: 99% quantile of the reference founders' distances",
            "keeps the broad continental group, still structured inside"], 0.6, 1.5, 5.6, 5.0, 17);
fig(s, "ancestry_round1/gate_pcs.png", 6.5, 1.4, 6.2, 5.4, "gate on AoU PCs");

s = content("Round 2: PCA on 1000G, cohort projected");
cite(s, "Privé et al. 2018 Bioinformatics (bigsnpr); Privé et al. 2020 Bioinformatics (OADP projection); Rousseeuw & Van Driessen 1999 Technometrics (MCD)");
bullets(s, ["PCA (bigsnpr bed_autoSVD) on 1000G founders of the superpopulation, 10 PCs",
            "clumped at r² 0.2, long-range LD regions removed",
            "round-1 set projected with OADP (bed_projectPCA)",
            "changed from plink2 --score: with ~500 training samples, simple projection shrinks PCs toward 0 past PC4",
            "gate: robust (MCD) Mahalanobis distance on PCs 1-4 around the target 1000G population(s), radius √χ²₄(0.999)"], 0.6, 1.5, 5.6, 5.2, 17);
fig(s, "round2_pca/cohort_on_reference.png", 6.5, 1.4, 6.2, 2.6, "cohort on the 1000G PCs");
fig(s, "ancestry_round2/gate_pcs.png", 6.5, 4.2, 6.2, 2.6, "kept vs dropped");

s = content("Round 2 PCA checks");
fig(s, "round2_pca/scree.png", 0.6, 1.4, 4.0, 2.6, "scree");
fig(s, "round2_pca/simple_vs_oadp.png", 4.75, 1.4, 4.0, 2.6, "simple vs OADP projection");
fig(s, "round2_pca/reference_scores.png", 8.9, 1.4, 3.8, 2.6, "1000G reference");
fig(s, "round2_pca/loadings.png", 0.6, 4.2, 12.1, 2.6, "|loadings| along the genome: no peaks");

s = content("How many people each gate keeps");
fig(s, "ancestry_round1/kept_vs_coverage.png", 0.6, 1.4, 4.0, 5.3, "round 1: kept vs gate coverage; plateau = cluster");
fig(s, "pop_gates/kept_by_population.png", 4.75, 1.4, 4.6, 5.3, "round 2: AoU people inside each 1000G population's gate");
fig(s, "pop_gates/overlap_EAS.png", 9.5, 1.4, 3.2, 5.3, "gate overlap (Jaccard), one per superpopulation");

s = content("The groups");
const rows = [["chb_chs", "East Asian, near CHB/CHS"], ["yri", "West African, near YRI"], ["lwk", "East African, near LWK"],
              ["stu_itu_beb", "South Asian, near STU/ITU/BEB"], ["pjl", "South Asian, near PJL"], ["asw_10k", "African American, admixed (10k sample)"],
              ["uniform_10k", "all ancestries, uniform sample (10k)"], ["even_10k", "all ancestries, even across groups (10k)"]];
const hd = (t) => ({ text: t, options: { bold: true, color: "1A1A1A", fill: { color: "F2F2F2" } } });
s.addTable([[hd("group"), hd("who"), hd("N"), hd("SNPs")], ...rows.map(([g, w]) => [g, w, "…", "…"])],
           { x: 0.6, y: 1.5, w: 8.2, colW: [2.0, 4.2, 1.0, 1.0], fontSize: 15, color: C.text1, border: { type: "solid", pt: 0.5, color: "C9D4DF" }, rowH: 0.52, valign: "middle" });
card(s, 9.2, 1.5, 3.5, 4.66, "mixed sets", "uniform_10k and even_10k mix ancestries on purpose: they show what structure does to each estimator");

s = content("QC");
card(s, 0.6, 1.4, 5.8, 2.6, "samples", "AoU srWGS callset (AoU's own sample QC); ancestry gates in rounds 1 and 2; relatives kept throughout");
card(s, 0.6, 4.1, 5.8, 2.6, "variants (GRM panel)", "HapMap3 autosomal SNPs; missingness ≤ 5%; MAF ≥ 1%; HWE p ≥ 1e-6 in single-ancestry groups, none in mixed sets; QC'd within each group");
card(s, 6.9, 1.4, 5.8, 2.6, "for PCA only", "long-range LD regions removed; clumped (r² 0.2 on 1000G, r² 0.1 within groups); autoSVD drops SNPs with loading peaks");
fig(s, "grm_panel/variants_by_step.png", 6.9, 4.1, 5.8, 2.6, "variants left after each QC step");

s = content("Covariate PCs, within each group");
cite(s, "Privé et al. 2018 Bioinformatics (bigsnpr, autoSVD); Privé et al. 2020 Bioinformatics (projection)");
bullets(s, ["bed_autoSVD on the group's HapMap3 common SNPs, clumped at r² 0.1, long-range LD removed",
            "fit on unrelated people: both members of every pair with AoU kinship > 0.1 left out",
            "relatives projected back (bed_projectSelfPCA): nobody dropped",
            "PCs 1-10: phenotype residualisation; PCs 1-4: PC-Relate"], 0.6, 1.5, 5.6, 5.2, 17);
fig(s, "covariate_pcs/covariate_pcs.png", 6.5, 1.4, 6.2, 2.6, "scree, fit vs projected");
fig(s, "covariate_pcs/covariate_pc_pairs.png", 6.5, 4.2, 6.2, 2.6, "PC pairs, density");

s = content("Genotypes against HWE");
eq(s, "F", 0.6, 1.5, 0.75);
bullets(s, ["F > 0: fewer heterozygotes than HWE (structure, inbreeding)", "expect ~0 in homogeneous groups, > 0 in mixed sets"], 0.6, 2.6, 5.4, 3.0, 16);
fig(s, "meta/f_by_maf.png", 6.4, 1.5, 6.3, 5.1, "all groups: F by MAF bin");

// ---------------------------------------------------------------- relatedness
section("Relatedness", "three estimators, one scale");
s = content("GRM");
cite(s, "Yang et al. 2010 Nat Genet; Yang et al. 2011 AJHG (GCTA)");
eq(s, "grm", 0.6, 1.6, 1.0);
bullets(s, ["x = 0/1/2 allele count, p = allele frequency in the group", "relative to the group's average: off-diagonal mean ≈ 0", "structure inside the group shows up as relatedness"], 0.6, 3.1, 6.5, 3.2);
card(s, 8.0, 1.6, 4.7, 3.0, "reference", "the group's allele frequencies, same for everyone");

s = content("PC-Relate");
cite(s, "Conomos, Reiner, Weir & Thornton 2016 AJHG");
eq(s, "pcrelate", 0.6, 1.5, 1.35);
eq(s, "pcrelate_mu", 0.6, 3.1, 0.5);
bullets(s, ["μ: each person's own allele frequency, from PCs 1-4", "removes ancestry-level sharing before measuring kinship", "relatedness = 2φ"], 0.6, 3.9, 6.8, 2.6);
card(s, 8.0, 1.6, 4.7, 3.0, "reference", "individual-specific allele frequencies from PCs");

s = content("Weir-Goudet");
cite(s, "Weir & Goudet 2017 Genetics; Zheng et al. 2012 Bioinformatics (SNPRelate)");
eq(s, "wg", 0.6, 1.5, 0.8);
eq(s, "wg_m", 0.6, 2.5, 0.6);
bullets(s, ["M: allele matching, no allele frequencies needed", "relative to the average pair in the sample", "relatedness = 2β (SNPRelate snpgdsIndivBeta)"], 0.6, 3.4, 6.8, 3.0);
card(s, 8.0, 1.6, 4.7, 3.0, "reference", "the mean pair in the sample");

s = content("All on one scale");
eq(s, "rel", 0.6, 1.5, 0.65);
card(s, 0.6, 2.6, 3.85, 3.6, "GRM", "group allele frequencies");
card(s, 4.74, 2.6, 3.85, 3.6, "PC-Relate", "per-person allele frequencies from PCs");
card(s, 8.85, 2.6, 3.85, 3.6, "Weir-Goudet", "mean pair in the sample");

s = content("Distributions");
fig(s, "meta/relatedness_distributions.png", 0.6, 1.4, 8.4, 5.4, "all groups: bulk and relatives");
fig(s, "meta/relatedness_spread.png", 9.3, 1.4, 3.4, 5.4, "SD, share > 0.1");

s = content("Estimators against each other");
fig(s, "weir_goudet/estimator_comparison.png", 0.6, 1.4, 12.1, 5.5, "chb_chs: hexes by PC distance, black = bin means");

s = content("Do they rank pairs the same way?");
fig(s, "weir_goudet/rank_comparison.png", 0.6, 1.4, 8.2, 5.5, "chb_chs: rank vs rank, density and PC distance");
fig(s, "meta/rank_agreement.png", 9.1, 1.4, 3.6, 5.5, "all groups: Spearman ρ");

s = content("Relatedness and PC distance");
fig(s, "weir_goudet/relatedness_vs_pc_distance.png", 0.6, 1.4, 6.0, 5.0, "relatedness by PC-distance bin");
s.addText("PC distance: Euclidean, round-2 PCs 1-4", { x: 0.6, y: 6.5, w: 12.1, h: 0.4, fontSize: 14, color: GREY, isTextBox: true });
fig(s, "weir_goudet/pc_distance_by_relatedness.png", 6.7, 1.4, 6.0, 5.0, "PC distance by relatedness bin");

s = content("Against IBD (chb_chs)");
cite(s, "Shemirani et al. 2021 Nat Commun (iLASH)");
fig(s, "ibd_vs_kinship/binned_by_ibd.png", 0.6, 1.4, 8.6, 5.5, "iLASH IBD relatedness vs each estimator");
bullets(s, ["IBD: shared segments from iLASH, no allele frequencies", "pairs without segments at 0"], 9.5, 1.6, 3.2, 4.0, 16);

// ---------------------------------------------------------------- covariance
section("Phenotypic covariance", "y_i y_j against relatedness");
s = content("Phenotypes");
bullets(s, ["height, BMI, systolic BP, HDL, LDL, platelets, hemoglobin", "most recent measurement per person", "residualised on age, sex and 10 PCs, then standardised", "every pair contributes y_i y_j"], 0.6, 1.5, 7.0, 3.6);
eq(s, "cov", 0.6, 5.3, 0.65);
card(s, 8.0, 1.5, 4.7, 3.4, "two scopes", "all pairs, and pairs with relatedness < 0.1 (the 'unrelated' slope)");

s = content("Covariance against relatedness");
fig(s, "kin_pheno_cov/binned_covariance_bulk.png", 0.6, 1.4, 6.0, 5.5, "bulk: relatedness to 0.1");
fig(s, "kin_pheno_cov/binned_covariance_full.png", 6.7, 1.4, 6.0, 5.5, "full: relatedness to 0.6, ±SE");

s = content("Slopes, intercepts, r²");
fig(s, "kin_pheno_cov/fit_coefficients.png", 0.6, 1.4, 8.0, 5.5, "chb_chs: per trait and estimator");
eq(s, "jk", 8.9, 1.6, 0.7);
bullets(s, ["delete-block jackknife: people in 50 random blocks", "drop every pair touching a block", "95% CI, t on 49 df", "Pearson and Spearman r²"], 8.9, 2.5, 3.8, 4.0, 15);

s = content("Across groups");
fig(s, "meta/fits_all_pairs.png", 0.6, 1.4, 6.0, 5.5, "all pairs; ◆ inverse-variance pooled");
fig(s, "meta/fits_unrelated.png", 6.7, 1.4, 6.0, 5.5, "relatedness < 0.1");

// ---------------------------------------------------------------- wrap
pres.addSection({ title: "Wrap-up" }); sec = "Wrap-up";
s = pres.addSlide({ masterName: "TITLE", sectionTitle: sec });
s.addText("Takeaways", { placeholder: "title" });
s.addText("…  ·  …  ·  …", { placeholder: "body" });

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
})();
