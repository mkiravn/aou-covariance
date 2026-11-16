"""Compile report figures and tables (only what stages already saved to report/) into one PDF.
    python scripts/report_pdf.py chb_chs uniform_10k asw_10k [--out FILE]"""
import argparse, datetime, json, os, sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["savefig.dpi"] = 400
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from aoucov import load_run

STAGES = [("weir_goudet", "Relatedness: GRM, PC-Relate, Weir-Goudet"), ("ibd_vs_kinship", "Relatedness vs iLASH IBD"), ("allele_variance", "Genotype variance vs 2p(1-p)"),
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

GCOL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#EE7733", "#000000"]
EST = ["GRM", "PC-Relate", "Weir-Goudet"]
MK = {"GRM": "o", "PC-Relate": "s", "Weir-Goudet": "^"}


def read_tables(cfgs, stage, name):
    out = []
    for run, cfg in cfgs.items():
        f = cfg.path("report", stage, name)
        if os.path.exists(f):
            out.append(pd.read_csv(f, sep="\t").assign(run=run))
    return pd.concat(out) if out else None


def save(pdf, fig, title, caption):
    fig.text(0.03, 0.975, title, fontsize=13, weight="bold", va="top")
    fig.text(0.03, 0.01, caption, fontsize=9, va="bottom")
    pdf.savefig(fig); plt.close(fig)


def meta_hist(pdf, h, col):
    h = h.assign(n=pd.to_numeric(h["n"], errors="coerce"))   # masked bins (1-20) not drawn
    fig, axes = plt.subplots(2, 3, figsize=(11, 8.5), sharey="row")
    fig.subplots_adjust(top=0.88, bottom=0.16, hspace=0.3, wspace=0.08)
    for c, e in enumerate(EST):
        for r, (lo, hi) in enumerate([(-0.07, 0.12), (0.1, 0.6)]):
            ax = axes[r, c]
            for run, g in h[h["estimator"] == e].groupby("run", sort=False):
                g = g[(g["lo"] >= lo - 1e-9) & (g["hi"] <= hi + 1e-9)]
                tot = h[(h["run"] == run) & (h["estimator"] == e)]["n"].sum()
                ax.stairs((g["n"] / tot / (g["hi"] - g["lo"])).to_numpy(), np.r_[g["lo"], g["hi"].iloc[-1:]], color=col[run], lw=1.3, label=run)
            ax.set(xlim=(lo, hi), xlabel=f"{e} relatedness", xticks=[-0.05, 0, 0.05, 0.1] if r == 0 else [0.2, 0.3, 0.4, 0.5])
            if r == 1: ax.set_yscale("log")
            if r == 0: ax.set_title(e)
    axes[0, 0].set_ylabel("density"); axes[1, 0].set_ylabel("density (log)")
    hd, lb = axes[0, 0].get_legend_handles_labels()
    fig.legend(hd, lb, loc="lower center", ncol=len(lb), fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.05))
    save(pdf, fig, "All groups: relatedness distributions", "Off-diagonal pairs; top: bulk, bottom: relatives (log). Bins with 1-20 pairs not drawn.")


def meta_spread(pdf, s, col):
    runs = list(dict.fromkeys(s["run"])); y = np.arange(len(runs))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.14, wspace=0.08, left=0.14)
    for ax, v, xl in [(axes[0], "sd", "SD of relatedness"), (axes[1], "share_gt_0.1", "share of pairs with relatedness > 0.1")]:
        for k, e in enumerate(EST):
            g = s[s["estimator"] == e].set_index("run").reindex(runs)
            ax.scatter(g[v], y + (k - 1) * 0.22, marker=MK[e], s=40, color=[col[r] for r in runs], edgecolor="k", lw=0.5, label=e, zorder=3)
        ax.set(xlabel=xl)
    axes[0].set_yticks(y, runs)
    axes[1].xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=2))
    axes[0].set_ylim(-0.6, len(runs) - 0.4)
    hd = [plt.Line2D([], [], marker=MK[e], ls="", color="0.6", mec="k", label=e) for e in EST]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.03))
    save(pdf, fig, "All groups: spread of relatedness", "Off-diagonal pairs, per estimator.")


def meta_rank(pdf, a, col):
    runs = list(dict.fromkeys(a["run"])); y = np.arange(len(runs))[::-1]
    a = a.assign(pr=a["y"] + " vs " + a["x"]); P = list(dict.fromkeys(a["pr"])); PM = ["o", "s", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.14, wspace=0.08, left=0.14)
    for ax, sc in zip(axes, ["all pairs", "relatedness < 0.1"]):
        for k, pr in enumerate(P):
            g = a[(a["pr"] == pr) & (a["pairs"] == sc)].set_index("run").reindex(runs)
            ax.scatter(g["spearman"], y + (k - 1) * 0.22, marker=PM[k % 3], s=40, color=[col[r] for r in runs], edgecolor="k", lw=0.5, zorder=3)
        ax.set(xlabel="Spearman rank correlation", title=sc)
    axes[0].set_yticks(y, runs); axes[0].set_ylim(-0.6, len(runs) - 0.4)
    hd = [plt.Line2D([], [], marker=PM[k % 3], ls="", color="0.6", mec="k", label=pr) for k, pr in enumerate(P)]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.03))
    save(pdf, fig, "All groups: rank agreement between estimators", "Spearman correlation of off-diagonal relatedness; right: pairs with relatedness < 0.1 by all three.")


def meta_variance(pdf, v, col):
    v = v.assign(snps=pd.to_numeric(v["snps"], errors="coerce")); v = v[v["snps"] >= 21]
    if "F_se" not in v: v = v.assign(F_se=v["het_ratio_se"])
    fig, ax = plt.subplots(figsize=(11, 8.5))
    fig.subplots_adjust(top=0.9, bottom=0.18, left=0.1, right=0.95)
    for run, g in v.groupby("run", sort=False):
        ax.errorbar(g["mean_maf"], g["F"], yerr=1.96 * g["F_se"], fmt="o-", ms=4, lw=1.3, capsize=2, color=col[run], label=run)
    ax.axhline(0, color="0.5", lw=0.8, ls=":")
    ax.set(xscale="log", xlabel="MAF (bin mean)", ylabel="F = 1 - het / 2p(1-p)")
    ax.set_xticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5], ["0.01", "0.02", "0.05", "0.1", "0.2", "0.5"]); ax.minorticks_off()
    ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=0.3)
    hd, lb = ax.get_legend_handles_labels()
    fig.legend(hd, lb, loc="lower center", ncol=min(len(lb), 8), fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.06))
    save(pdf, fig, "All groups: heterozygote deficit F by allele frequency", "GRM panel SNPs in MAF bins (21+ SNPs); 95% CI across SNPs. F = variance / 2p(1-p) - 1; single-ancestry panels are HWE-filtered.")

def pooled(g, v):
    g = g.dropna(subset=[v, v + "_se"]); g = g[g[v + "_se"] > 0]
    if g.empty: return np.nan, np.nan, np.nan
    w = 1 / g[v + "_se"] ** 2; m = (w * g[v]).sum() / w.sum(); se = w.sum() ** -0.5
    return m, m - 1.96 * se, m + 1.96 * se


def meta_coef(pdf, t, col, scope):
    t = t[t["pairs"] == scope]
    ph = list(dict.fromkeys(t["pheno"])); runs = list(dict.fromkeys(t["run"])); R = len(runs)
    M = [("slope", "slope"), ("intercept", "intercept"), ("r2", "Pearson r2"), ("spearman_r2", "Spearman r2")]
    fig, axes = plt.subplots(4, 3, figsize=(11, 8.5), sharex=True, sharey="row")
    fig.subplots_adjust(top=0.92, bottom=0.17, hspace=0.12, wspace=0.06, left=0.1, right=0.99)
    for c, e in enumerate(EST):
        for r, (v, lab) in enumerate(M):
            ax = axes[r, c]
            if v not in t: ax.set_visible(False); continue
            for k, run in enumerate(runs):
                g = t[(t["estimator"] == e) & (t["run"] == run)].set_index("pheno").reindex(ph)
                x = np.arange(len(ph)) - 0.08 + (k - (R - 1) / 2) * 0.55 / max(R - 1, 1)
                ax.errorbar(x, g[v], yerr=[(g[v] - g[v + "_lo"]).clip(lower=0), (g[v + "_hi"] - g[v]).clip(lower=0)],
                            fmt="o", ms=3, lw=0.8, color=col[run], label=run)
            for xi, p_ in enumerate(ph):
                m, lo, hi = pooled(t[(t["estimator"] == e) & (t["pheno"] == p_)], v)
                ax.errorbar(xi + 0.36, m, yerr=[[m - lo], [hi - m]], fmt="D", ms=5, color="k", lw=1.2, label="pooled" if xi == 0 else None)
            if v.endswith("r2"): ax.set_yscale("log")
            elif v != "intercept": ax.axhline(0, color="0.5", lw=0.6)
            if r == 0: ax.set_title(e)
            if c == 0: ax.set_ylabel(lab)
            ax.set_xticks(range(len(ph)), ph, rotation=40, ha="right", fontsize=8)
            for xi in range(1, len(ph)): ax.axvline(xi - 0.5, color="0.9", lw=0.6, zorder=0)
    hd, lb = axes[0, 0].get_legend_handles_labels()
    fig.legend(hd, lb, loc="lower center", ncol=min(len(lb), 9), fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.03))
    save(pdf, fig, f"All groups: fits of y_i y_j on relatedness ({scope})",
         "raw, base_pcs; 95% delete-block jackknife CIs; black diamond: inverse-variance-weighted mean across groups.")


def meta_pages(pdf, cfgs):
    col = {r: GCOL[k % len(GCOL)] for k, r in enumerate(cfgs)}
    h = read_tables(cfgs, "weir_goudet", "relatedness_hist.tsv")
    if h is not None: meta_hist(pdf, h, col)
    s = read_tables(cfgs, "weir_goudet", "relatedness_spread.tsv")
    if s is not None: meta_spread(pdf, s, col)
    r = read_tables(cfgs, "weir_goudet", "rank_agreement.tsv")
    if r is not None: meta_rank(pdf, r, col)
    v = read_tables(cfgs, "allele_variance", "by_maf.tsv")
    if v is not None: meta_variance(pdf, v, col)
    t = read_tables(cfgs, "kin_pheno_cov", "estimates.tsv")
    if t is not None:
        t = t[(t["transform"] == "raw") & (t["covset"] == "base_pcs")] if "base_pcs" in set(t["covset"]) else t
        for scope in dict.fromkeys(t["pairs"]):
            meta_coef(pdf, t, col, scope)
    if h is None or s is None:
        print("meta: no relatedness_hist/spread tables yet (rerun weir_goudet)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default=os.path.expanduser(f"~/aucov_runs/report_{datetime.date.today():%Y%m%d}.pdf"))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with PdfPages(a.out) as pdf:
        text_page(pdf, "Relatedness and phenotypic covariance", [f"runs: {', '.join(a.runs)}", f"compiled {datetime.datetime.now():%Y-%m-%d %H:%M}",
                  "Aggregate figures only: bins and hexagons with 21+ pairs; counts 1-20 masked.", "Relatedness = 2 x kinship."])
        cfgs = {run: load_run(run) for run in a.runs}
        meta_pages(pdf, cfgs)
        for run, cfg in cfgs.items():
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
