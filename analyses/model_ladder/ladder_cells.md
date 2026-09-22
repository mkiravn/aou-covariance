# Model ladder

Reads one run's merged cross-products and asks whether the relatedness slope
differs between unrelated and related pairs, before and after allowing an
offset for each degree of relatedness. Set `SAMPLE_SET` and `MERGED_DIR`;
everything else follows.

All slopes pass through the origin. U = a < 0.02, T = 0.02–0.05, R = 0.05–0.7.
Degree bands: deg4 0.05–0.09, deg3 0.09–0.18, deg2 0.18–0.36, deg1 0.36–0.7
(KING cutoffs rounded to bin edges).

| estimator | pairs | model | reported |
|---|---|---|---|
| `h2_Unrel` | U | y = β a | β |
| `h2_UnrelInt` | U | y = β a + c | β; `Unrel.intercept` = c |
| `h2_OneSlope` | U+T+R | y = β a | β |
| `h2_Rel` | U+T+R | y = β_U a·[U] + β_R a·[T+R] | β_R |
| `h2_OneSlopeOffsets` | U+R | y = β a + Σ o_k·[band k] | β; `OneSlopeOffsets.degK` |
| `h2_RelOffsets` | U+R | y = β_U a·[U] + β_R a·[R] + Σ o_k·[band k] | β_R; `RelOffsets.degK` |
| `diff_Rel-Unrel` | | β_R − β_U, no offsets | |
| `diff_RelOffsets-Unrel` | | β_R − β_U, with offsets | |
| `diff_Rel-RelOffsets` | | how much of β_R the offsets remove | |

- U and R share no bins, so β_U equals `h2_Unrel` exactly in both two-slope
  models.
- T is left out of the offset models because offsets are defined only above 0.05.
- With an offset in every band, the related slope of `h2_RelOffsets` comes only
  from variation in a within each degree. That is robust to shared environment
  but has a large SE, and within full sibs it still absorbs some dominance
  (IBD2 rises with a).
- Differences are formed inside each jackknife replicate, so their SEs carry
  the covariance of the two terms. With 50 blocks, read z against t with 49 df.
- Jackknife blocks are random over individuals and split families, so SEs of
  the related-range estimators are somewhat optimistic.

Compute: 4 vCPU is enough; the estimators run in seconds per model.

---

## Cell 1 — config

```python
import os, sys, glob, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

SAMPLE_SET = "eur_D2"
MERGED_DIR = os.path.expanduser("~/grm_pheno_cov_eur_D2")     # *_merged.{full,jk}.tsv
BINS       = os.path.expanduser("~/bins/bins_wide.txt")        # the bins the run was merged with
REPO       = os.path.expanduser("~/repos/AOU-covariance")
OUT        = os.path.expanduser(f"~/model_ladder/{SAMPLE_SET}")
DEST       = ("gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
              f"/phenotypic_covariance_v9/analyses/model_ladder/{SAMPLE_SET}")
os.makedirs(f"{OUT}/plots", exist_ok=True)

BIN_TAG, CLASS_TAG, NBLOCKS = "wide", "deg1", 50
COVSETS    = ["base", "base_pcs", "base_pcs_zip3", "base_pcs_zip3_ses"]
TRANSFORMS = ["raw", "invnorm"]
REF_TF, REF_CS = "invnorm", "base_pcs"
FIT_TF, MIN_N  = "invnorm", 20

CATEGORY = {
    **dict.fromkeys(["height", "weight", "bmi", "waist_circumference",
                     "hip_circumference", "waist_hip_ratio"], "anthropometric"),
    **dict.fromkeys(["systolic_bp", "diastolic_bp", "heart_rate"], "cardiovascular"),
    **dict.fromkeys(["glucose", "hba1c", "hdl_cholesterol", "ldl_cholesterol",
                     "triglycerides", "total_cholesterol"], "metabolic"),
    **dict.fromkeys(["hemoglobin", "wbc_count", "platelet_count", "rbc_count",
                     "hematocrit", "mch", "mcv", "mchc", "neutrophil_count",
                     "lymphocyte_count", "monocyte_count", "basophil_count",
                     "eosinophil_count"], "blood cells"),
    **dict.fromkeys(["creatinine", "alt"], "liver & kidney"),
    **dict.fromkeys(["alcohol_audit_c_score", "cigarettes_per_day", "pss_score",
                     "social_support_score", "loneliness_score", "eds_score"],
                    "behavioural"),
}
CAT_ORDER = ["anthropometric", "cardiovascular", "metabolic",
             "blood cells", "liver & kidney", "behavioural"]
CAT_COL = {"anthropometric": "#1f77b4", "cardiovascular": "#d62728",
           "metabolic": "#2ca02c", "blood cells": "#9467bd",
           "liver & kidney": "#8c564b", "behavioural": "#ff7f0e",
           "uncategorised": "#7f7f7f"}
CS_MARK = {"base": "o", "base_pcs": "s", "base_pcs_zip3": "D", "base_pcs_zip3_ses": "^"}

sys.path.insert(0, f"{REPO}/04_process_shards/src")
import importlib, he_estimators as HE
importlib.reload(HE)
assert hasattr(HE, "DEG_BANDS"), "he_estimators.py predates the ladder -- git pull the repo"
print(SAMPLE_SET, MERGED_DIR, OUT)
```

## Cell 2 — estimators for every merged model

```python
MID = HE.load_grid(BINS)
tags = sorted(os.path.basename(p)[:-len("_merged.full.tsv")]
              for p in glob.glob(f"{MERGED_DIR}/*__{CLASS_TAG}__{BIN_TAG}_merged.full.tsv"))
print(f"{len(tags)} merged models")

parts = []
for i, tag in enumerate(tags, 1):
    jk = f"{MERGED_DIR}/{tag}_merged.jk.tsv"
    if not os.path.isfile(jk):
        print(f"  skip {tag}: no .jk.tsv")
        continue
    S, N = HE.load_arrays(f"{MERGED_DIR}/{tag}_merged.full.tsv", jk, len(MID), NBLOCKS)
    r = HE.jackknife(S, N, MID, NBLOCKS)
    r["phenotype"], r["transform"], r["covset"] = tag.split("__")[:3]
    parts.append(r)
    if i % 25 == 0 or i == len(tags):
        print(f"  [{i}/{len(tags)}]")

E = pd.concat(parts, ignore_index=True)
E["category"] = E["phenotype"].map(CATEGORY).fillna("uncategorised")
E["z"] = E["est"] / E["se"]
E.to_csv(f"{OUT}/estimators_{SAMPLE_SET}.tsv", sep="\t", index=False)
print(f"{E['phenotype'].nunique()} phenotypes -> {OUT}/estimators_{SAMPLE_SET}.tsv")

LADDER = ["h2_Unrel", "h2_OneSlope", "h2_Rel", "h2_OneSlopeOffsets", "h2_RelOffsets",
          "diff_Rel-Unrel", "diff_RelOffsets-Unrel", "diff_Rel-RelOffsets"]
(E[(E["transform"] == REF_TF) & (E["covset"] == REF_CS) & (E["classes"] == "noPO")]
 .pivot_table(index="phenotype", columns="estimator", values="est")[LADDER].round(3))
```

## Cell 3 — style and plotting helpers

Colour does one job per figure: covariate set is an ordered blue ramp, pair
class is the first three categorical hues, trait category is spatial grouping
rather than colour. Hollow markers are estimates within one SE of zero.

```python
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CLASS_COL = {"other": "#2a78d6", "FS": "#eb6834", "PO": "#1baf7a"}
CLASS_MARK = {"other": "o", "FS": "s", "PO": "^"}
RAMP = ["#86b6ef", "#3987e5", "#256abf", "#104281"]          # ordered: light to dark
CS_COL = dict(zip(COVSETS, RAMP))
DIVERGING = LinearSegmentedColormap.from_list("bpr", ["#104281", "#f0efec", "#e34948"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
    "axes.titlecolor": INK, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": INK, "font.size": 9,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.titlesize": 11,
})

ref = (E[(E["estimator"] == "h2_Unrel") & (E["transform"] == REF_TF) & (E["covset"] == REF_CS)]
       .set_index("phenotype")["est"])
ORDER = sorted(E["phenotype"].unique(),
               key=lambda p: (CAT_ORDER.index(CATEGORY[p]) if p in CATEGORY else 99,
                              -ref.get(p, -np.inf)))
YPOS = {p: len(ORDER) - 1 - i for i, p in enumerate(ORDER)}


def dress_y(ax, ylabels=True):
    """Phenotype rows, grouped by trait category with a hairline between groups."""
    ax.set_yticks([YPOS[p] for p in ORDER])
    ax.set_yticklabels(ORDER if ylabels else [""] * len(ORDER))
    ax.set_ylim(-0.7, len(ORDER) - 0.3)
    ax.tick_params(axis="y", length=0)
    for prev, p in zip(ORDER, ORDER[1:]):
        if CATEGORY.get(p) != CATEGORY.get(prev):
            ax.axhline(YPOS[p] + 0.5, color=GRID, lw=0.8, zorder=0)


def label_categories(ax):
    """Category names in the free right margin, one per group."""
    for cat in CAT_ORDER:
        rows = [YPOS[p] for p in ORDER if CATEGORY.get(p) == cat]
        if rows:
            ax.text(1.015, (min(rows) + max(rows)) / 2, cat, transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=8, color=MUTED, clip_on=False)


def cs_handles():
    return [Line2D([], [], lw=0, marker="o", ms=6, color=CS_COL[c], label=c) for c in COVSETS]


def forest(D, title, xlabel, fname):
    """One row per phenotype, one panel per transform, covariate sets dodged and
    coloured light to dark. Filled markers are more than one SE from zero."""
    dodge = np.linspace(-0.28, 0.28, len(COVSETS))
    fig, axes = plt.subplots(1, len(TRANSFORMS), figsize=(13, 0.26 * len(ORDER) + 2.2),
                             sharey=True, sharex=True)
    axes = np.atleast_1d(axes)
    for ax, tf in zip(axes, TRANSFORMS):
        ax.grid(axis="x", zorder=0)
        ax.set_axisbelow(True)
        for cs, off in zip(COVSETS, dodge):
            s = D[(D["transform"] == tf) & (D["covset"] == cs)
                  & D["phenotype"].isin(YPOS) & np.isfinite(D["est"])]
            if s.empty:
                continue
            y = s["phenotype"].map(YPOS) + off
            se = s["se"].fillna(0)
            ax.hlines(y, s["est"] - se, s["est"] + se, colors=CS_COL[cs], lw=1.1, alpha=0.55)
            sig = (s["est"].abs() > se).to_numpy()
            ax.scatter(s["est"][sig], y[sig], color=CS_COL[cs], s=18, lw=0, zorder=3)
            ax.scatter(s["est"][~sig], y[~sig], facecolors=SURFACE, edgecolors=CS_COL[cs],
                       s=18, lw=0.9, zorder=3)
        ax.axvline(0, color=AXIS, lw=0.9)
        dress_y(ax)
        ax.set_title(tf)
        ax.set_xlabel(xlabel)
    label_categories(axes[-1])
    fig.legend(handles=cs_handles(), loc="upper center", bbox_to_anchor=(0.5, 0.02),
               ncol=len(COVSETS), title="covariate set", title_fontsize=8)
    fig.suptitle(f"{title} — {SAMPLE_SET}", y=1.0)
    plt.tight_layout()
    plt.savefig(f"{OUT}/plots/{fname}.png", dpi=150, bbox_inches="tight")
    plt.show()
```

## Cell 4 — slope differences

```python
NOPO = E[E["classes"] == "noPO"]
forest(NOPO[NOPO["estimator"] == "diff_Rel-Unrel"],
       r"related minus unrelated slope, no offsets ($\hat h^2_{Rel} - \hat h^2_{Unrel}$)",
       "difference (± jackknife SE)", "diff_rel_unrel")
forest(NOPO[NOPO["estimator"] == "diff_RelOffsets-Unrel"],
       r"related minus unrelated slope, with degree offsets ($\hat h^2_{RelOffsets} - \hat h^2_{Unrel}$)",
       "difference (± jackknife SE)", "diff_reloffsets_unrel")

(NOPO[(NOPO["transform"] == REF_TF) & (NOPO["covset"] == REF_CS)
      & NOPO["estimator"].str.startswith("diff_")]
 .pivot_table(index="phenotype", columns="estimator", values="z").round(1))
```

## Cell 5 — degree offsets

One panel per trait category, one line per phenotype, distant relatives on the
left. Each point is that band's intercept: covariance above what the shared
slope predicts.

```python
DEG = list(HE.DEG_BANDS)
O = NOPO[(NOPO["transform"] == REF_TF) & (NOPO["covset"] == REF_CS)
         & NOPO["estimator"].isin([f"RelOffsets.{d}" for d in DEG])].copy()
O["x"] = O["estimator"].str.split(".").str[1].map(DEG.index)
cats = [c for c in CAT_ORDER if c in set(O["category"])]

fig, axes = plt.subplots(1, len(cats), figsize=(2.9 * len(cats), 3.4), sharey=True)
for ax, cat in zip(np.atleast_1d(axes), cats):
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    for ph, s in O[O["category"] == cat].groupby("phenotype"):
        s = s.sort_values("x")
        ax.errorbar(s["x"], s["est"], yerr=s["se"], color=RAMP[2], lw=1.1, alpha=0.55,
                    elinewidth=0.6, marker="o", ms=3, capsize=0)
    ax.axhline(0, color=AXIS, lw=0.9)
    ax.set_xticks(range(len(DEG)))
    ax.set_xticklabels(DEG)
    ax.set_xlim(-0.35, len(DEG) - 0.65)
    ax.set_title(cat)
np.atleast_1d(axes)[0].set_ylabel("offset (cross-product scale)")
fig.suptitle(f"degree offsets of $h^2_{{RelOffsets}}$ — {SAMPLE_SET}, {REF_TF}, {REF_CS}, noPO")
plt.tight_layout()
plt.savefig(f"{OUT}/plots/reloffsets_by_degree.png", dpi=150, bbox_inches="tight")
plt.show()
```

Offsets that fall with degree are shared environment or dominance that decays
with relatedness. Offsets near zero leave the related slope to carry everything.

## Cell 5b — statistical support for each rung

Two questions, two panels. Left: is the *coefficient* resolved — each contrast
over its jackknife SE. Right: does the *fit* improve — the change in chi-square
per bin when the term is added, also over a jackknife SE, refitting both rungs
inside every replicate.

The table adds each rung's absolute misfit and a cross-validated error, where
every replicate is scored on the pairs it left out. Read the CV column with
care: each replicate holds out only 2/`NBLOCKS` of the pairs, so it has little
power to separate rungs.

```python
CONTRASTS = ["diff_Rel-Unrel", "diff_RelOffsets-Unrel", "diff_Rel-RelOffsets"]

Z = (NOPO[(NOPO["transform"] == REF_TF) & (NOPO["covset"] == REF_CS)
          & NOPO["estimator"].isin(CONTRASTS)]
     .pivot_table(index="phenotype", columns="estimator", values="z").reindex(ORDER)[CONTRASTS])

fit_rows, fit_z = [], {}
for ph in ORDER:
    tag = f"{ph}__{REF_TF}__{REF_CS}__{CLASS_TAG}__{BIN_TAG}"
    S, N = HE.load_arrays(f"{MERGED_DIR}/{tag}_merged.full.tsv",
                          f"{MERGED_DIR}/{tag}_merged.jk.tsv", len(MID), NBLOCKS)
    f_, c_ = HE.compare_models(S, N, MID, NBLOCKS)
    fit_rows.append(f_.assign(phenotype=ph))
    fit_z[ph] = c_.set_index("contrast")["z"]
FIT = pd.concat(fit_rows, ignore_index=True)
FITZ = pd.DataFrame(fit_z).T.reindex(ORDER)
FIT.to_csv(f"{OUT}/model_fit.tsv", sep="\t", index=False)
FITZ.to_csv(f"{OUT}/model_fit_z.tsv", sep="\t")


def zmap(ax, D, title, ylabels):
    v = np.nanmax(np.abs(D.to_numpy())) or 1
    im = ax.imshow(D.to_numpy(), aspect="auto", cmap=DIVERGING,
                   norm=TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v))
    ax.set_xticks(range(D.shape[1]))
    ax.set_xticklabels([c.replace("diff_", "") for c in D.columns], rotation=25, ha="right")
    ax.set_yticks(range(len(D)))
    ax.set_yticklabels(D.index if ylabels else [""] * len(D))
    ax.tick_params(length=0)
    ax.set_title(title)
    for k in range(1, len(D)):
        ax.axhline(k - 0.5, color=SURFACE, lw=0.6)
    return im


fig, axes = plt.subplots(1, 2, figsize=(11, 0.26 * len(ORDER) + 2.4))
zmap(axes[0], Z, "coefficient: estimate / SE", True)
im = zmap(axes[1], FITZ, "fit: change in chi-square / SE", False)
fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, label="z")
fig.suptitle(f"statistical support — {SAMPLE_SET}, {REF_TF}, {REF_CS}, noPO")
plt.savefig(f"{OUT}/plots/model_support.png", dpi=150, bbox_inches="tight")
plt.show()

print(FIT.pivot(index="phenotype", columns="model", values="chi2_per_bin")
      .reindex(ORDER).round(2).to_string())
```

## Cell 6 — per-phenotype fits with a coefficient panel

`FIT_TF` only, all covariate sets. Colour is the covariate set, marker shape the
pair class. Fitted lines are drawn for one covariate set so they stay readable.
Small vertical numbers above the bins are pair counts, from `LABEL_FROM` upward;
bins below `MIN_N` are not plotted at all.

```python
MS, XLIM = 4.5, (-0.06, 1.12)
FIT_CS = REF_CS
LABEL_FROM = HE.U_HI            # -np.inf to label every bin
H2_ROWS = ["h2_Unrel", "h2_OneSlope", "h2_Rel", "h2_OneSlopeOffsets",
           "h2_RelOffsets", "h2_FS", "h2_PedW25"]
B2_ROWS = ["b2_FS", "b2_step"]
EF = NOPO[NOPO["transform"] == FIT_TF].set_index(["phenotype", "covset", "estimator"])[["est", "se"]]


def coef(ph, cs, name):
    k = (ph, cs, name)
    return tuple(EF.loc[k]) if k in EF.index else (np.nan, np.nan)


def coef_panel(ax, ph, rows):
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)
    for cs, dy in zip(COVSETS, np.linspace(-0.26, 0.26, len(COVSETS))):
        for i, name in enumerate(rows):
            est, se = coef(ph, cs, name)
            if not np.isfinite(est):
                continue
            y = len(rows) - 1 - i + dy
            sig = not (np.isfinite(se) and abs(est) <= se)
            ax.hlines(y, est - (se or 0), est + (se or 0), colors=CS_COL[cs], lw=1.1, alpha=0.55)
            ax.scatter(est, y, s=18, zorder=3, lw=0 if sig else 0.9,
                       color=CS_COL[cs] if sig else SURFACE,
                       edgecolors=CS_COL[cs])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows[::-1])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.axvline(0, color=AXIS, lw=0.9)
    ax.tick_params(axis="y", length=0)


for ph in ORDER:
    tables = {}
    for cs in COVSETS:
        f = f"{MERGED_DIR}/{ph}__{FIT_TF}__{cs}__{CLASS_TAG}__{BIN_TAG}_merged.full.tsv"
        if os.path.isfile(f):
            t = pd.read_csv(f, sep="\t")
            tables[cs] = t[(t.full_n >= MIN_N) & (t["class"] != "pooled")]
    if not tables:
        continue

    fig = plt.figure(figsize=(14, 6.2))
    gs = GridSpec(2, 2, width_ratios=[3, 1.1], height_ratios=[len(H2_ROWS), len(B2_ROWS) + 0.6],
                  wspace=0.3, hspace=0.45)
    ax, axh, axb = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)

    for cs, t in tables.items():
        for cls in ("other", "FS", "PO"):
            s = t[t["class"] == cls]
            if len(s):
                ax.errorbar(s.bin_midpoint, s.full_mean, yerr=s.jk_se, fmt=CLASS_MARK[cls],
                            ms=MS, color=CS_COL[cs], lw=0, elinewidth=0.6, capsize=0,
                            ecolor=CS_COL[cs], alpha=0.9, mec=SURFACE, mew=0.4,
                            zorder=2 if cls == "other" else 4)

    # pair counts per bin, over the fitted range, from the covariate set the fits use
    ref_t = tables.get(FIT_CS, next(iter(tables.values())))
    lab = ref_t[(ref_t.bin_midpoint >= LABEL_FROM) & (ref_t.bin_midpoint < HE.R_HI)]
    span = np.ptp(np.r_[[t.full_mean.to_numpy() for t in tables.values()]])
    for _, row in lab.iterrows():
        ax.text(row.bin_midpoint, row.full_mean + row.jk_se + 0.015 * span, f"{int(row.full_n):,}",
                rotation=90, ha="center", va="bottom", fontsize=5, color=MUTED, zorder=6)

    g = lambda name: coef(ph, FIT_CS, name)[0]
    xs = np.array(XLIM)
    fits = []
    if np.isfinite(g("h2_Unrel")):
        fits.append((xs, g("h2_Unrel") * xs, "--", r"$h^2_{Unrel}$"))
    if np.isfinite(g("h2_Rel")):
        hi = np.array([HE.U_HI, HE.R_HI])
        fits.append((hi, g("h2_Rel") * hi, "-.", r"$h^2_{Rel}$"))
    if np.isfinite(g("h2_PedW25")):
        a = np.linspace(0.05, HE.R_HI, 100)
        fits.append((a, g("h2_PedW25") * a + g("PedW25.quad") * a ** 2, (0, (3, 1, 1, 1)),
                     r"$h^2_{Ped,W25}$"))
    for k, (d, (lo, hi)) in enumerate(HE.DEG_BANDS.items()):
        off = g(f"RelOffsets.{d}")
        if np.isfinite(g("h2_RelOffsets")) and np.isfinite(off):
            a = np.array([lo, hi])
            fits.append((a, g("h2_RelOffsets") * a + off, "-",
                         r"$h^2_{RelOffsets}$ + band" if k == 0 else None))
    for x, y, ls, lab_ in fits:
        ax.plot(x, y, ls=ls, lw=1.0, color=INK2, alpha=0.7, zorder=5, label=lab_)

    for lo, _ in HE.DEG_BANDS.values():
        ax.axvline(lo, color=GRID, lw=0.7, zorder=0)
    ax.axhline(0, color=AXIS, lw=0.9)
    ax.set_xlim(*XLIM)
    ax.set_xlabel(r"GRM relatedness  $a_{ij}$")
    ax.set_ylabel(r"mean phenotype cross-product  $\overline{y_i y_j}$")

    coef_panel(axh, ph, H2_ROWS)
    coef_panel(axb, ph, B2_ROWS)
    axh.set_title("estimates ± jackknife SE")

    cls_h = [Line2D([], [], lw=0, marker=CLASS_MARK[c], ms=MS + 1, color=INK2, label=c)
             for c in ("other", "FS", "PO")]
    fit_h = [Line2D([], [], ls=ls, lw=1.0, color=INK2, label=lab_)
             for _, _, ls, lab_ in fits if lab_]
    ax.legend(handles=cls_h + fit_h, loc="upper left", ncol=2)
    fig.legend(handles=cs_handles(), loc="lower center", bbox_to_anchor=(0.42, -0.04),
               ncol=len(COVSETS), title="covariate set", title_fontsize=8)
    fig.suptitle(f"{ph} — {SAMPLE_SET}, {FIT_TF}, noPO; bins with n ≥ {MIN_N}; "
                 f"fits for {FIT_CS}; small numbers are pair counts")
    plt.savefig(f"{OUT}/plots/{ph}_fits_{FIT_TF}.png", dpi=150, bbox_inches="tight")
    plt.show()
```

## Cell 7 — copy to the bucket

`gcloud storage cp` does not fail loudly on a partial object, so every file is
checked against its remote size.

```python
def remote_sizes(prefix):
    out = subprocess.run(["gcloud", "storage", "ls", "-l", f"{prefix}/"],
                         capture_output=True, text=True).stdout
    return {f[-1].rsplit("/", 1)[1]: int(f[0]) for f in map(str.split, out.splitlines())
            if len(f) >= 3 and f[0].isdigit() and f[-1].startswith("gs://")}

for sub, files in (("", glob.glob(f"{OUT}/*.tsv")), ("plots", glob.glob(f"{OUT}/plots/*.png"))):
    dst = f"{DEST}/{sub}".rstrip("/")
    subprocess.run(["gcloud", "storage", "cp", "-I", f"{dst}/"], input="\n".join(files),
                   text=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    remote = remote_sizes(dst)
    bad = [os.path.basename(f) for f in files if remote.get(os.path.basename(f)) != os.path.getsize(f)]
    print(f"{sub or 'tables':<7} {len(files) - len(bad)}/{len(files)} verified" + (f"  MISMATCH {bad}" if bad else ""))
```
