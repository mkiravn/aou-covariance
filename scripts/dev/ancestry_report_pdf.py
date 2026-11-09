"""PCA / gate figures from every run's report folder into one PDF (figures already mask cells under 21).
    python scripts/dev/ancestry_report_pdf.py [out.pdf] [run ...]
"""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from aoucov import load_run

STAGES = ["pop_gates", "ancestry_round1", "round2_pca", "ancestry_round2", "covariate_pcs"]   # pop_gates: report_pop_gates.ipynb
out = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/ancestry_report.pdf")
runs = sys.argv[2:] or ["chb_chs", "yri", "lwk", "stu_itu_beb", "pjl", "asw", "mxl", "uniform_10k", "even_10k"]

def text_page(pdf, title, lines):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.05, 0.94, title, fontsize=16, weight="bold")
    fig.text(0.05, 0.9, "\n".join(lines), fontsize=9, family="monospace", va="top")
    pdf.savefig(fig); plt.close(fig)

with PdfPages(out) as pdf:
    for run in runs:
        rdir = os.path.expanduser(load_run(run).path("report"))
        for stage in STAGES:
            sj = os.path.join(rdir, stage, "summary.json")
            if not os.path.exists(sj):
                continue
            s = json.load(open(sj))
            lines = [f"{k}: {v}" for k, v in (s.get("notes") or {}).items()]
            lines += ["", "checks:"] + [f"  {'PASS' if c.get('ok') else 'FAIL'}  {c.get('check')}  {c.get('detail') or ''}"
                                         for c in s.get("checks") or []]
            for f in s.get("files") or []:   # tables are written masked by Report.table
                tsv = os.path.join(rdir, stage, f["file"])
                if f["kind"] == "table" and os.path.exists(tsv):
                    lines += ["", f["file"] + ":"] + ["  " + l.replace("\t", "  ") for l in open(tsv).read().splitlines()[:40]]
            text_page(pdf, f"{run} / {stage}", lines)
            for f in s.get("files") or []:
                png = os.path.join(rdir, stage, f["file"])
                if f["kind"] != "figure" or not os.path.exists(png):
                    continue
                img = plt.imread(png)
                h, w = img.shape[:2]
                fig = plt.figure(figsize=(11, 8.5))
                ax = fig.add_axes([0.03, 0.03, 0.94, 0.86]); ax.imshow(img); ax.axis("off")
                fig.text(0.03, 0.95, f"{run} / {stage} / {f['file']}", fontsize=13, weight="bold")
                fig.text(0.03, 0.92, f.get("caption", ""), fontsize=10)
                pdf.savefig(fig); plt.close(fig)
        print(run, "done")
print("->", out)
