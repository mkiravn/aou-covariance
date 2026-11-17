// GRM vs Weir-Goudet, height only, 5 slides + appendix: node slides/build_grm_wg.js OUT.pptx FIGDIR
// FIGDIR is the unzipped slide bundle after `python slides/grm_wg_figs.py FIGDIR`.
const path = require("path");
const OUT = process.argv[2] || path.join(__dirname, "grm_wg.pptx");
const D = require("./deck_lib")(process.argv[3]);
const { pres, C, GREY, GROUPS, section, content, fig, bullets, cite } = D;
const sz = require("image-size");
const EQ = (k) => path.join(__dirname, "eq", `${k}.png`);
const KEY = ["chb_chs", "pjl", "asw_10k", "uniform_10k"].filter((g) => GROUPS.includes(g) || !process.argv[3]);
const label = (s, t, x, y, w) => s.addText(t, { x, y, w, h: 0.3, fontSize: 14, bold: true, color: "444444", margin: 0, isTextBox: true });
const foot = (s, t) => s.addText(t, { x: 0.6, y: 6.95, w: 11.6, h: 0.35, fontSize: 12, color: GREY, isTextBox: true });
const forGroup = (g, cb) => { D.setRun(g); cb(); D.setRun(GROUPS[0]); };
const OLD = (f) => !!process.argv[3] && !D.exists(f);   // figure not made yet: crop the existing per-group figures (they still show PC-Relate)
const PROVISIONAL = OLD("@grm_wg/height_fits.png");
function tiles(s, file, crop, groups) {
  groups.forEach((g, k) => forGroup(g, () => {
    if (!D.exists(file)) return;
    const x = 0.6 + (k % 4) * 3.05, y = 1.25 + Math.floor(k / 4) * 2.85;
    label(s, g, x, y, 2.9); fig(s, file, x, y + 0.3, 2.9, 2.45, "", crop);
  }));
}

pres.addSection({ title: "Main" }); D.setSec("Main");

let s = content("GRM and Weir-Goudet");
{
  const F = (26 / 30) / 300, WMAX = 6.6;
  const put = (k, x, y) => { const d = sz(EQ(k)), f = Math.min(F, WMAX / d.width); s.addImage({ path: EQ(k), x, y, w: d.width * f, h: d.height * f, altText: k }); return d.height * f; };
  const rows = [["GRM", ["grm"], "the group's allele frequencies p", "Yang et al. 2011 AJHG (GCTA)"],
                ["Weir-Goudet", ["wg", "wg_m"], "the average pair of distinct individuals in the sample (S*); no allele frequencies", "Weir & Goudet 2026 Genetics, Box 3"]];
  let y = 1.5;
  for (const [t, eqs, ref, cit] of rows) {
    let h = 0; for (const e of eqs) h += put(e, 2.4, y + h) + 0.12;
    h = Math.max(h, 1.2);
    s.addText(t, { x: 0.6, y, w: 1.7, h: 0.5, fontSize: 22, bold: true, color: C.text1, margin: 0, isTextBox: true });
    s.addText([{ text: "reference: ", options: { bold: true } }, { text: ref, options: { breakLine: true } }, { text: cit, options: { fontSize: 12, color: "888888" } }],
              { x: 9.2, y, w: 3.5, h, fontSize: 16, color: "333333", margin: 0, valign: "top", isTextBox: true });
    y += h + 0.35;
    s.addShape(pres.shapes.LINE, { x: 0.6, y: y - 0.18, w: 12.1, h: 0, line: { color: "DDDDDD", width: 0.75 } });
  }
  s.addText("relatedness = 2 × kinship for both; GRM ≈ Weir-Goudet with each person's mean kinship removed from every row and column (Weir & Goudet 2026, Eq. 13-14)",
            { x: 0.6, y: y + 0.05, w: 12.1, h: 0.8, fontSize: 16, color: "333333", isTextBox: true });
}

s = content("Weir-Goudet against GRM");
KEY.forEach((g, k) => forGroup(g, () => {
  const x = 0.6 + (k % 2) * 6.15, y = 1.25 + Math.floor(k / 2) * 2.85;
  label(s, g, x, y, 5.9);
  if (D.exists("weir_goudet/wg_vs_grm.png")) fig(s, "weir_goudet/wg_vs_grm.png", x, y + 0.3, 5.95, 2.5, "");
  else { fig(s, "weir_goudet/estimator_comparison.png", x, y + 0.3, 2.9, 2.5, "", { l: 0.32, r: 0.38, b: 0.5, t: 0.03 });
         fig(s, "weir_goudet/estimator_comparison.png", x + 3.0, y + 0.3, 2.9, 2.5, "", { l: 0.32, r: 0.38, t: 0.52 }); }
}));
foot(s, "left: bulk, right: full range; cells (21+ pairs) by mean PC distance; black: bin means; dotted: y = x; other groups in the appendix");

s = content("Relatedness against PC distance");
if (OLD("@grm_wg/relatedness_vs_pc_distance.png")) tiles(s, "weir_goudet/relatedness_vs_pc_distance.png", { r: 0.5 }, GROUPS); else fig(s, "@grm_wg/relatedness_vs_pc_distance.png", 0.6, 1.25, 12.1, 5.6, "");
foot(s, "mean relatedness in 40 equal-count bins of PC distance (Euclidean, round-2 PCs 1-4)");

s = content("Height: covariance against relatedness");
if (OLD("@grm_wg/covariance_bulk.png")) tiles(s, "kin_pheno_cov/binned_covariance_bulk.png", { r: 0.76, b: 0.5 }, GROUPS); else fig(s, "@grm_wg/covariance_bulk.png", 0.6, 1.25, 12.1, 5.6, "");
foot(s, "bin means of y_i y_j (residualised, standardised height) with least-squares fits over all pairs; relatedness to 0.1");

s = content(PROVISIONAL ? "Slopes by trait (height: top row)" : "Height: slope and intercept");
if (PROVISIONAL) ["chb_chs", "pjl", "asw_10k", "uniform_10k"].filter((g) => GROUPS.includes(g)).forEach((g, k) => forGroup(g, () => {
  label(s, g, 0.6 + k * 3.05, 1.25, 2.9);
  fig(s, "kin_pheno_cov/fit_coefficients.png", 0.6 + k * 3.05, 1.55, 2.9, 5.3, "", { r: 0.515, b: 0.535 });
})); else fig(s, "@grm_wg/height_fits.png", 0.6, 1.25, 12.1, 5.6, "");
foot(s, "slope = ĥ²; 95% delete-block jackknife CIs; lwk and yri left out (too few pairs)");

section("Appendix", "");
for (const g of GROUPS.filter((g) => !KEY.includes(g))) forGroup(g, () => {
  const a = content(`Weir-Goudet against GRM: ${g}`);
  if (D.exists("weir_goudet/wg_vs_grm.png")) fig(a, "weir_goudet/wg_vs_grm.png", 0.6, 1.3, 12.1, 5.6, "");
  else { fig(a, "weir_goudet/estimator_comparison.png", 0.6, 1.3, 6.0, 5.6, "", { l: 0.32, r: 0.38, b: 0.5, t: 0.03 }); fig(a, "weir_goudet/estimator_comparison.png", 6.7, 1.3, 6.0, 5.6, "", { l: 0.32, r: 0.38, t: 0.52 }); }
});
s = content("Height: covariance against relatedness, full range");
if (OLD("@grm_wg/covariance_full.png")) tiles(s, "kin_pheno_cov/binned_covariance_full.png", { r: 0.76, b: 0.5 }, GROUPS); else fig(s, "@grm_wg/covariance_full.png", 0.6, 1.25, 12.1, 5.6, "");
foot(s, "bin means ± 1.96 SE, relatedness to 0.6");

for (const f of ["relatedness_vs_pc_distance", "covariance_bulk", "height_fits"]) if (OLD(`@grm_wg/${f}.png`)) console.log(`provisional: no grm_wg/${f}.png, cropped the old figure instead`);
D.write(OUT);
