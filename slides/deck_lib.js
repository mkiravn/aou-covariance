// Shared layouts and helpers for the decks: const D = require("./deck_lib")(figdir)
module.exports = function (FIGDIR) {
  const fs = require("fs"), path = require("path");
  const pptxgen = require("pptxgenjs");
  const { applyTheme } = require(process.env.PPTX_SKILL + "/scripts/apply_theme.js");

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
  function perGroup(title, draw) {   // one slide per group
    for (const g of GROUPS) { RUN = g; draw(content(`${title}: ${g}`)); }
    RUN = EX;
  }
  const content = (t) => { const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: sec }); s.addText(t, { placeholder: "title" }); return s; };

  const USED = new Set();
  function place(s, f, file, x, y, w, h, crop) {   // fit inside the box, keep aspect; crop = { l, t, r, b } fractions trimmed (ImageMagick)
    const c = crop || CROP[file];
    if (c) {
      const d0 = require("image-size")(f), q = { l: 0, t: 0, r: 0, b: 0, ...c };
      const cw = Math.round(d0.width * (1 - q.l - q.r)), ch = Math.round(d0.height * (1 - q.t - q.b));
      const out = path.join(require("os").tmpdir(), "deck_crops", `${key(file).replace(/[\/]/g, "__")}_${q.l}_${q.t}_${q.r}_${q.b}.png`);
      fs.mkdirSync(path.dirname(out), { recursive: true });
      require("child_process").execFileSync("magick", [f, "-crop", `${cw}x${ch}+${Math.round(d0.width * q.l)}+${Math.round(d0.height * q.t)}`, "+repage", out]);
      f = out;
    }
    const sz = require("image-size")(f), k = Math.min(w / sz.width, h / sz.height), iw = sz.width * k, ih = sz.height * k;
    s.addImage({ path: f, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih, objectName: `fig:${file}`, altText: file });
  }
  const CROP = {};   // file -> { l, t, r, b } fractions to trim, e.g. "round2_pca/loadings.png": { t: 0, b: 0.5 }
  const ORDER = ["chb_chs", "yri", "lwk", "stu_itu_beb", "pjl", "asw_10k", "uniform_10k", "even_10k"];
  const GROUPS = FIGDIR && fs.existsSync(FIGDIR)
    ? fs.readdirSync(FIGDIR).filter((d) => d !== "meta" && fs.existsSync(path.join(FIGDIR, d, "weir_goudet")))   // runs used only for borrowed figures (e.g. asw) are not groups
        .sort((a, b) => (ORDER.indexOf(a) + 1 || 99) - (ORDER.indexOf(b) + 1 || 99) || a.localeCompare(b))
    : ["<group>"];
  const EX = GROUPS.includes("chb_chs") ? "chb_chs" : GROUPS[0];   // example group for one-off method slides
  let RUN = EX;
  const key = (file) => (file.startsWith("meta/") ? file : file.startsWith("@") ? file.slice(1) : `${RUN}/${file}`);   // "@run/stage/fig.png": another run's figure
  function fig(s, file, x, y, w, h, note, crop) {
    USED.add(key(file));
    const f = FIGDIR && path.join(FIGDIR, key(file));
    if (f && fs.existsSync(f)) {
      place(s, f, file, x, y, w, h, crop);
      return;
    }
    s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: "FFFFFF" }, line: { color: "BBBBBB", width: 1, dashType: "dash" }, objectName: `fig:${file}` });
    s.addText([{ text: file, options: { fontSize: 14, color: "555555", breakLine: !!note } },
               ...(note ? [{ text: note, options: { fontSize: 12, color: "888888" } }] : [])],
              { x, y, w, h, align: "center", valign: "middle", isTextBox: true });
  }
  const hd = (t) => ({ text: t, options: { bold: true, color: "1A1A1A", fill: { color: "F2F2F2" } } });
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

  function supplementary() {
    // ---------------------------------------------------------------- supplementary: every other report figure, by group
    if (FIGDIR && fs.existsSync(FIGDIR)) {
      for (const g of GROUPS) {
        const root = path.join(FIGDIR, g);
        const extra = fs.readdirSync(root).filter((d) => fs.statSync(path.join(root, d)).isDirectory()).sort()
          .flatMap((d) => fs.readdirSync(path.join(root, d)).filter((f) => f.endsWith(".png")).sort().map((f) => `${d}/${f}`))
          .filter((f) => !USED.has(`${g}/${f}`));
        if (!extra.length) continue;
        section(`Supplementary: ${g}`, "other report figures");
        for (const f of extra) {
          const s = content(`${g}: ` + f.replace(".png", "").replace("/", " / ").replace(/_/g, " "));
          place(s, path.join(root, f), f, 0.6, 1.3, 12.1, 5.6);
        }
      }
    }

  }
  async function write(OUT) {
    await pres.writeFile({ fileName: OUT });
    await applyTheme(OUT, THEME);
    console.log("wrote", OUT);
  }
  const exists = (file) => !!FIGDIR && fs.existsSync(path.join(FIGDIR, key(file)));
  return { exists, pres, C, GREY, GROUPS, EX, USED, CROP, section, content, perGroup, place, fig, eq, bullets, card, cite, hd, supplementary, write,
           setSec: (t) => { sec = t; }, getSec: () => sec, setRun: (g) => { RUN = g; } };
};
