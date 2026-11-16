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
function tiles(title, choices, foot) {   // one cropped panel per group, 4 x 2; choices: [[file, crop], ...] first that exists
  const s = content(title);
  GROUPS.slice(0, 8).forEach((g, k) => forGroup(g, () => {
    const x = 0.6 + (k % 4) * 3.05, y = 1.3 + Math.floor(k / 4) * 2.85;
    label(s, g, x, y, 2.9);
    const [file, crop] = choices.find(([f]) => D.exists(f)) || choices[0];
    fig(s, file, x, y + 0.3, 2.9, 2.45, "", crop);
  }));
  if (foot) s.addText(foot, { x: 0.6, y: 6.95, w: 11.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });
  return s;
}

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

s = content("Round-2 gates: N per 1000G population, and overlap");
fig(s, "pop_gates/kept_by_population.png", 0.6, 1.3, 12.1, 3.0, "AoU people inside each 1000G population's gate");
["EAS", "AFR", "SAS"].forEach((sp, k) => fig(s, `pop_gates/overlap_${sp}.png`, 0.6 + k * 4.1, 4.4, 3.9, 2.5, `${sp}: gate overlap (Jaccard)`));
s.addText("same gate as round 2 (MCD on each population's founders, PCs 1-4, q = 0.999); gates overlap, so rows don't add up; bars under 21 not drawn",
          { x: 0.6, y: 6.95, w: 11.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });

tiles("Structure in each sample", [["ancestry_round2/gate_pcs.png", { r: 0.5 }], ["round2_pca/cohort_on_reference.png", { r: 0.5 }], ["covariate_pcs/covariate_pc_pairs.png", { r: 0.8 }]],
      "each sample on round-2 PCs 1-2 against 1000G (gated groups: kept vs dropped; mixed sets: whole sample; asw_10k: its own covariate PCs 1-2); cells with 21+ people");

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
{
  const F = (24 / 30) / 300, WMAX = 6.4;   // equation images are 30 pt at 300 dpi: shown at 24 pt
  const put = (k, x, y) => { const d = sz(EQ(k)), f = Math.min(F, WMAX / d.width); s.addImage({ path: EQ(k), x, y, w: d.width * f, h: d.height * f, altText: k, objectName: `eq:${k}` }); return d.height * f; };
  const rows = [["GRM", ["grm"], "the group's allele frequencies", "Yang et al. 2011 AJHG (GCTA)"],
                ["PC-Relate", ["pcrelate", "pcrelate_mu"], "each person's own allele frequencies, from PCs 1-4", "Conomos et al. 2016 AJHG"],
                ["Weir-Goudet", ["wg", "wg_m"], "the average pair of distinct individuals in the run (S*); no allele frequencies", "Weir & Goudet 2026 Genetics, Box 3"]];
  let y = 1.3;
  for (const [t, eqs, ref, cit] of rows) {
    let h = 0; for (const e of eqs) h += put(e, 2.3, y + h) + 0.1;
    h = Math.max(h, 1.0);
    s.addText(t, { x: 0.6, y, w: 1.6, h: 0.5, fontSize: 20, bold: true, color: C.text1, margin: 0, isTextBox: true });
    s.addText([{ text: "reference: ", options: { bold: true } }, { text: ref, options: { breakLine: true } }, { text: cit, options: { fontSize: 12, color: "888888" } }],
              { x: 9.0, y, w: 3.7, h, fontSize: 15, color: "333333", margin: 0, valign: "top", isTextBox: true });
    y += h + 0.2;
    s.addShape(pres.shapes.LINE, { x: 0.6, y: y - 0.1, w: 12.1, h: 0, line: { color: "DDDDDD", width: 0.75 } });
  }
  put("rel", 2.3, y + 0.05);
  s.addText("relatedness = 2 × kinship for all three", { x: 9.0, y: y + 0.05, w: 3.7, h: 0.5, fontSize: 15, color: GREY, margin: 0, isTextBox: true });
}

s = content("Distributions");
fig(s, "meta/relatedness_distributions.png", 0.6, 1.3, 7.4, 5.6, "all groups: bulk (top), relatives (bottom, log)");
fig(s, "meta/relatedness_spread.png", 8.2, 1.3, 4.5, 5.6, "SD and share > 0.1");

const PER = 2;   // the three biplots, two groups per slide, full width: bulk, then the full range
for (const [part, crop, note] of [["", { t: 0.02, b: 0.5 }, "bulk (−0.07 to 0.12)"], [", full range", { t: 0.525 }, "full range (−0.1 to 0.6)"]]) {
  for (let i = 0; i < GROUPS.length; i += PER) {
    const s = content(`Estimators against each other${part}${GROUPS.length > PER ? ` (${i / PER + 1}/${Math.ceil(GROUPS.length / PER)})` : ""}`);
    GROUPS.slice(i, i + PER).forEach((g, k) => forGroup(g, () => {
      const y = 1.25 + k * 2.85;
      label(s, g, 0.6, y, 5.9);
      fig(s, "weir_goudet/estimator_comparison.png", 0.6, y + 0.3, 12.1, 2.5, "", crop);
    }));
    s.addText(`${note}; left to right: PC-Relate vs GRM, Weir-Goudet vs GRM, Weir-Goudet vs PC-Relate; cells by mean PC distance, black = bin means`,
              { x: 0.6, y: 6.95, w: 11.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });
  }
}

s = content("Do they rank pairs the same way?");
fig(s, "meta/rank_agreement.png", 0.6, 1.3, 12.1, 5.6, "all groups: Spearman ρ");

tiles("Relatedness against PC distance", [["weir_goudet/relatedness_vs_pc_distance.png", { r: 0.5 }]],
      "mean relatedness in 40 equal-count bins of PC distance (Euclidean, round-2 PCs 1-4); GRM, PC-Relate, Weir-Goudet");

s = content("Relatedness against IBD (chb_chs)");
D.setRun("chb_chs");
fig(s, "ibd_vs_kinship/binned_by_ibd.png", 0.6, 1.3, 8.6, 5.5, "iLASH IBD vs each estimator");
D.setRun(EX);
bullets(s, ["IBD: iLASH segments, no allele frequencies", "pairs without segments at 0", "cells: pair density (21+ pairs); black: bin means"], 9.4, 1.5, 3.3, 4.0, 15);
cite(s, "Shemirani et al. 2021 Nat Commun (iLASH)");

s = content("Phenotypes and model");
bullets(s, ["height, BMI, systolic BP, HDL, LDL, platelets, hemoglobin; most recent measurement",
            "residualised on age, sex and 10 covariate PCs, then standardised",
            "every pair contributes y_i y_j; fit on all pairs and on pairs with relatedness < 0.1",
            "delete-block jackknife: people in 50 random blocks, drop pairs touching a block; 95% CI, t on 49 df",
            "Pearson and Spearman r²"], 0.6, 1.4, 6.8, 5.3, 16);
eqw(s, "cov", 7.8, 1.7, 4.9, 0.55);
eqw(s, "h2", 7.8, 2.6, 4.9, 0.5);
s.addText("x_ij: relatedness (2 × kinship) of i and j; y standardised, so the slope estimates h²", { x: 7.8, y: 3.3, w: 4.9, h: 0.8, fontSize: 14, color: GREY, margin: 0, isTextBox: true });

tiles("Covariance against relatedness: height, to 0.1", [["kin_pheno_cov/binned_covariance_bulk.png", { r: 0.76, b: 0.5 }]],
      "bin means of y_i y_j for GRM, PC-Relate and Weir-Goudet, with fits; other traits per group in the appendix");
tiles("Covariance against relatedness: height, to 0.6", [["kin_pheno_cov/binned_covariance_full.png", { r: 0.76, b: 0.5 }]],
      "full range, ±SE; other traits per group in the appendix");

s = content("Slopes and intercepts across groups");
fig(s, "meta/fits_all_pairs.png", 0.6, 1.25, 12.1, 5.5, "slope, intercept", { b: 0.56 });
s.addText("all pairs; 95% jackknife CIs, clipped at the axis; lwk and yri left out (too few pairs); colours as on the next slide", { x: 0.6, y: 6.95, w: 11.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });

s = content("r² across groups");
fig(s, "meta/fits_all_pairs.png", 0.6, 1.25, 12.1, 5.7, "Pearson and Spearman r²", { t: 0.48 });


// ---------------------------------------------------------------- appendix: per group, then everything else
section("Appendix", "per group: covariance for every trait");
for (const g of GROUPS) forGroup(g, () => {
  const a = content(`${g}: covariance against relatedness, to 0.1`);
  fig(a, "kin_pheno_cov/binned_covariance_bulk.png", 0.6, 1.3, 12.1, 5.6);
});
D.write(OUT);
