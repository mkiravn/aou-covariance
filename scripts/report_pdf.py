"""Compile report figures and tables (only what stages already saved to report/) into one PDF.
    python scripts/report_pdf.py chb_chs uniform_10k asw_10k [--out FILE]"""
import argparse, datetime, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from aoucov import load_run

STAGES = [("weir_goudet", "Relatedness: GRM, PC-Relate, Weir-Goudet"), ("ibd_vs_kinship", "Relatedness vs iLASH IBD"),
          ("kin_pheno_cov", "Phenotypic covariance vs relatedness")]


def text_page(pdf, title, lines=(), size=11):
    fig = plt.figure(figsize=(11, 8.5)); fig.text(0.06, 0.9, title, fontsize=18, weight="bold")
    for k, l in enumerate(lines):
        fig.text(0.06, 0.84 - k * 0.032, l, fontsize=size, family="monospace" if size < 11 else None)
    pdf.savefig(fig); plt.close(fig)


def image_page(pdf, path, title, caption):
    img = plt.imread(path)
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0.03, 0.08, 0.94, 0.84]); ax.imshow(img); ax.axis("off")
    fig.text(0.03, 0.95, title, fontsize=13, weight="bold")
    fig.text(0.03, 0.03, caption, fontsize=9, wrap=True)
    pdf.savefig(fig); plt.close(fig)


def coef_table(pdf, run, path):
    t = pd.read_csv(path, sep="\t")
    t = t[(t["transform"] == "raw") & (t["covset"] == "base_pcs")] if "base_pcs" in set(t["covset"]) else t
    f = lambda r, v, fmt: f"{r[v]:{fmt}}  [{r[v + '_lo']:{fmt}}, {r[v + '_hi']:{fmt}}]"
    rows = [[r["pheno"], r["estimator"], r["pairs"], f(r, "slope", ".3f"), f(r, "intercept", ".1e"), f(r, "r2", ".1e"),
             f(r, "spearman_r2", ".1e") if "spearman_r2" in t else ""] for _, r in t.iterrows()]
    for k in range(0, len(rows), 28):
        fig = plt.figure(figsize=(11, 8.5)); ax = fig.add_axes([0.03, 0.04, 0.94, 0.86]); ax.axis("off")
        tb = ax.table(rows[k:k + 28], colLabels=["phenotype", "estimator", "pairs", "slope [95% CI]", "intercept [95% CI]", "Pearson r2 [95% CI]", "Spearman r2 [95% CI]"],
                      colWidths=[0.09, 0.09, 0.12, 0.15, 0.18, 0.185, 0.185], loc="upper center", cellLoc="left")
        tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.3)
        for (r_, c_), cell in tb.get_celld().items():
            cell.set_edgecolor("0.8")
            if r_ == 0: cell.set_text_props(weight="bold"); cell.set_facecolor("0.93")
        ax.set_title(f"{run}: fits of y_i y_j on relatedness (raw, base_pcs)", loc="left", fontsize=12, weight="bold")
        pdf.savefig(fig); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default=os.path.expanduser(f"~/aucov_runs/report_{datetime.date.today():%Y%m%d}.pdf"))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with PdfPages(a.out) as pdf:
        text_page(pdf, "Relatedness and phenotypic covariance", [f"runs: {', '.join(a.runs)}", f"compiled {datetime.datetime.now():%Y-%m-%d %H:%M}",
                  "Aggregate figures only: bins and hexagons with 21+ pairs; counts 1-20 masked.", "Relatedness = 2 x kinship (full sibs 0.5)."])
        for run in a.runs:
            cfg = load_run(run)
            for stage, title in STAGES:
                d = cfg.path("report", stage)
                if not os.path.exists(f"{d}/summary.json"):
                    print(f"{run}/{stage}: no report, skipped"); continue
                s = json.load(open(f"{d}/summary.json"))
                if s.get("simulated"):
                    print(f"{run}/{stage}: simulated data")
                for f in s["files"]:
                    if f["kind"] == "figure" and os.path.exists(f"{d}/{f['file']}"):
                        image_page(pdf, f"{d}/{f['file']}", f"{run}: {title}", f["caption"])
                if stage == "kin_pheno_cov" and os.path.exists(f"{d}/estimates.tsv"):
                    coef_table(pdf, run, f"{d}/estimates.tsv")   # the report copy: counts masked
    print("wrote", a.out)


if __name__ == "__main__":
    main()
