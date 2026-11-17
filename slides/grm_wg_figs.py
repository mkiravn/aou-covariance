"""GRM vs Weir-Goudet, height only: slide-sized figures from the report tables in a slide bundle.
    python slides/grm_wg_figs.py FIGDIR        (writes FIGDIR/grm_wg/*.png)"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = sys.argv[1]
OUT = os.path.join(FIG, "grm_wg"); os.makedirs(OUT, exist_ok=True)
ORDER = ["chb_chs", "yri", "lwk", "stu_itu_beb", "pjl", "asw_10k", "uniform_10k", "even_10k"]
SMALL = {"lwk", "yri"}   # too few pairs for fits
TRAIT = "height"
COL = {"GRM": "#E69F00", "Weir-Goudet": "#009E73"}
MK = {"GRM": "o", "Weir-Goudet": "^"}
plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.3,
                     "savefig.dpi": 300, "font.family": "sans-serif"})


def table(g, stage, name):
    f = os.path.join(FIG, g, stage, name)
    return pd.read_csv(f, sep="\t") if os.path.exists(f) else None


groups = [g for g in ORDER if os.path.isdir(os.path.join(FIG, g, "weir_goudet"))] + \
         sorted(g for g in os.listdir(FIG) if g not in ORDER and os.path.isdir(os.path.join(FIG, g, "weir_goudet")))
num = lambda s: pd.to_numeric(s, errors="coerce")


def grid(n):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.4), squeeze=False)
    for ax in axes.flat[n:]: ax.set_visible(False)
    return fig, axes.flat


def legend(fig, extra=()):
    hd = [plt.Line2D([], [], color=COL[e], marker=MK[e], lw=2, label=e) for e in COL] + list(extra)
    fig.legend(handles=hd, loc="lower center", ncol=len(hd), frameon=False, bbox_to_anchor=(0.5, -0.01))


# 1. mean relatedness against PC distance
fig, axes = grid(len(groups))
for ax, g in zip(axes, groups):
    t = table(g, "weir_goudet", "relatedness_vs_pc_distance.tsv")
    if t is None: ax.set_visible(False); continue
    for e in COL:
        d = t[t["estimator"] == e]
        ax.plot(d["pc_distance"], d["mean_rel"], marker=MK[e], ms=4, lw=1.8, color=COL[e])
    ax.axhline(0, color="k", lw=0.6); ax.set_title(g, fontweight="bold")
for ax in list(axes)[:len(groups)]:
    ax.set_xlabel("PC distance"); ax.set_ylabel("mean relatedness")
fig.tight_layout(rect=(0, 0.05, 1, 1)); legend(fig)
fig.savefig(f"{OUT}/relatedness_vs_pc_distance.png"); plt.close(fig)

# 2. height: mean y_i y_j against relatedness, to 0.1
for name, xmax, ylim in [("covariance_bulk", 0.1, (-0.04, 0.06)), ("covariance_full", 0.6, None)]:
    fig, axes = grid(len(groups))
    for ax, g in zip(axes, groups):
        b = table(g, "kin_pheno_cov", "binned_covariance.tsv"); est = table(g, "kin_pheno_cov", "estimates.tsv")
        if b is None: ax.set_visible(False); continue
        b = b[b["pheno"] == TRAIT]
        for e in COL:
            d = b[(b["estimator"] == e) & (b["hi"] <= xmax + 1e-9)]
            ax.errorbar(d["mean_rel"], d["mean_cov"], yerr=None if name.endswith("bulk") else 1.96 * d["se"], fmt=MK[e], ms=5,
                        color=COL[e], mec="white", mew=0.5, elinewidth=0.8, zorder=2)
            if est is not None:
                r = est[(est["pheno"] == TRAIT) & (est["estimator"] == e) & (est["pairs"] == "all pairs") & (est["transform"] == "raw")]
                r = r[r["covset"] == "base_pcs"] if "base_pcs" in set(r["covset"]) else r
                if len(r):
                    xx = np.array([-0.1, xmax]); ax.plot(xx, r["intercept"].iloc[0] + r["slope"].iloc[0] * xx, color=COL[e], lw=1.6, zorder=1)
        ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6); ax.set_xlim(-0.1 if xmax < 0.2 else -0.1, xmax)
        if ylim: ax.set_ylim(*ylim)
        ax.set_title(g, fontweight="bold"); ax.set_xlabel("relatedness"); ax.set_ylabel("mean $y_i y_j$")
    fig.tight_layout(rect=(0, 0.05, 1, 1)); legend(fig)
    fig.savefig(f"{OUT}/{name}.png"); plt.close(fig)

# 3. height: slope (h2) and intercept by group
rows = []
for g in groups:
    est = table(g, "kin_pheno_cov", "estimates.tsv")
    if est is None or g in SMALL: continue
    r = est[(est["pheno"] == TRAIT) & (est["transform"] == "raw") & est["estimator"].isin(list(COL))]
    r = r[r["covset"] == "base_pcs"] if "base_pcs" in set(r["covset"]) else r
    rows.append(r.assign(run=g))
if rows:
    t = pd.concat(rows); gs = list(dict.fromkeys(t["run"])); y = np.arange(len(gs))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.4), sharey=True)
    for ax, (v, lab) in zip(axes, [("slope", "slope = ${\\hat{h}}^{2}$"), ("intercept", "intercept")]):
        for yy in y[::2]: ax.axhspan(yy - 0.5, yy + 0.5, color="0.95", zorder=0)
        for k, e in enumerate(COL):
            for scope, fill, off in [("all pairs", True, -0.15), (t["pairs"][~t["pairs"].eq("all pairs")].iloc[0] if (~t["pairs"].eq("all pairs")).any() else None, False, 0.15)]:
                if scope is None: continue
                d = t[(t["estimator"] == e) & (t["pairs"] == scope)].set_index("run").reindex(gs)
                yy = y + (k - 0.5) * 0.36 + off * 0.5
                ax.errorbar(d[v], yy, xerr=[(d[v] - d[v + "_lo"]).clip(lower=0), (d[v + "_hi"] - d[v]).clip(lower=0)], fmt=MK[e], ms=7,
                            color=COL[e], mfc=COL[e] if fill else "white", mew=1.5, lw=1.4, capsize=0, zorder=3)
        ax.axvline(0, color="k", lw=0.6); ax.set_xlabel(lab); ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y, gs); axes[0].set_ylim(-0.6, len(gs) - 0.4)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    legend(fig, [plt.Line2D([], [], color="0.4", marker="o", ls="", label="all pairs"),
                 plt.Line2D([], [], color="0.4", marker="o", mfc="white", ls="", label="relatedness < 0.1")])
    fig.savefig(f"{OUT}/height_fits.png"); plt.close(fig)

print("wrote", OUT, sorted(os.listdir(OUT)))
