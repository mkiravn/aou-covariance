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

## Cell 3 — plotting helpers

```python
ref = (E[(E["estimator"] == "h2_Unrel") & (E["transform"] == REF_TF) & (E["covset"] == REF_CS)]
       .set_index("phenotype")["est"])
ORDER = sorted(E["phenotype"].unique(),
               key=lambda p: (CAT_ORDER.index(CATEGORY[p]) if p in CATEGORY else 99,
                              -ref.get(p, -np.inf)))
YPOS = {p: len(ORDER) - 1 - i for i, p in enumerate(ORDER)}


def dress_y(ax):
    ax.set_yticks([YPOS[p] for p in ORDER])
    ax.set_yticklabels(ORDER, fontsize=8)
    for lab in ax.get_yticklabels():
        lab.set_color(CAT_COL[CATEGORY.get(lab.get_text(), "uncategorised")])
    ax.set_ylim(-0.7, len(ORDER) - 0.3)
    for prev, p in zip(ORDER, ORDER[1:]):
        if CATEGORY.get(p) != CATEGORY.get(prev):
            ax.axhline(YPOS[p] + 0.5, color="0.85", lw=0.8, zorder=0)


def cat_handles():
    present = [c for c in CAT_ORDER if c in set(E["category"])]
    return [Line2D([], [], lw=0, marker="s", ms=8, color=CAT_COL[c], label=c) for c in present]


def forest(D, title, xlabel, fname):
    """One row per phenotype, one panel per transform, covariate sets dodged
    and marker-coded, colour = trait category."""
    offsets = np.linspace(-0.3, 0.3, len(COVSETS))
    fig, axes = plt.subplots(1, len(TRANSFORMS), figsize=(15, 0.3 * len(ORDER) + 2),
                             sharey=True, sharex=True)
    axes = np.atleast_1d(axes)
    for ax, tf in zip(axes, TRANSFORMS):
        for cs, off in zip(COVSETS, offsets):
            s = D[(D["transform"] == tf) & (D["covset"] == cs)
                  & D["phenotype"].isin(YPOS) & np.isfinite(D["est"])]
            if s.empty:
                continue
            y = s["phenotype"].map(YPOS) + off
            cols = s["category"].map(CAT_COL).tolist()
            se = s["se"].fillna(0)
            ax.hlines(y, s["est"] - se, s["est"] + se, colors=cols, lw=0.8, alpha=0.5)
            ax.scatter(s["est"], y, c=cols, marker=CS_MARK[cs], s=24,
                       edgecolors="0.25", linewidths=0.4, zorder=3)
        ax.axvline(0, color="grey", lw=0.6)
        dress_y(ax)
        ax.set_title(tf)
        ax.set_xlabel(xlabel)
    cs_h = [Line2D([], [], lw=0, marker=CS_MARK[c], ms=6, color="0.4", label=c) for c in COVSETS]
    fig.legend(handles=cat_handles() + cs_h, loc="upper center",
               bbox_to_anchor=(0.5, 0.0), ncol=5, fontsize=8, frameon=False)
    fig.suptitle(f"{title} — {SAMPLE_SET}")
    plt.tight_layout()
    plt.savefig(f"{OUT}/plots/{fname}.png", dpi=130, bbox_inches="tight")
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

```python
DEG = list(HE.DEG_BANDS)                 # deg4 .. deg1
O = NOPO[(NOPO["transform"] == REF_TF) & (NOPO["covset"] == REF_CS)
         & NOPO["estimator"].isin([f"RelOffsets.{d}" for d in DEG])].copy()
O["x"] = O["estimator"].str.split(".").str[1].map(DEG.index)
cats = [c for c in CAT_ORDER if c in set(O["category"])]
fig, axes = plt.subplots(1, len(cats), figsize=(3.2 * len(cats), 3.6), sharey=True)
for ax, cat in zip(np.atleast_1d(axes), cats):
    for ph, s in O[O["category"] == cat].groupby("phenotype"):
        s = s.sort_values("x")
        ax.errorbar(s["x"], s["est"], yerr=s["se"], color=CAT_COL[cat], lw=0.9,
                    elinewidth=0.5, alpha=0.8, marker="o", ms=3)
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_xticks(range(len(DEG)))
    ax.set_xticklabels(DEG)
    ax.set_title(cat, color=CAT_COL[cat], fontsize=10)
np.atleast_1d(axes)[0].set_ylabel("offset (cross-product scale)")
fig.suptitle(f"degree offsets of $h^2_{{RelOffsets}}$ — {SAMPLE_SET}, {REF_TF}, {REF_CS}, noPO")
plt.tight_layout()
plt.savefig(f"{OUT}/plots/reloffsets_by_degree.png", dpi=130, bbox_inches="tight")
plt.show()
```

Offsets that fall with degree are shared environment or dominance that decays
with relatedness. Offsets near zero leave the related slope to carry everything.

## Cell 6 — per-phenotype fits with a coefficient panel

`FIT_TF` only, all covariate sets. Left: binned means (hue = pair class, shade
= covariate set), the `h2_Unrel` additive line dashed, and the `h2_RelOffsets`
fit drawn band by band. Right: every estimator for every covariate set.

```python
CMAP  = {"other": plt.cm.Blues, "FS": plt.cm.Oranges, "PO": plt.cm.Purples}
MARK  = {"other": "o", "FS": "s", "PO": "^"}
SHADE = dict(zip(COVSETS, np.linspace(0.45, 0.95, len(COVSETS))))
MS, XLIM = 5, (-0.06, 1.12)
H2_ROWS = ["h2_Unrel", "h2_OneSlope", "h2_Rel", "h2_OneSlopeOffsets",
           "h2_RelOffsets", "h2_FS", "h2_PedW25"]
B2_ROWS = ["b2_FS", "b2_step"]
EF = NOPO[NOPO["transform"] == FIT_TF].set_index(["phenotype", "covset", "estimator"])[["est", "se"]]


def coef(ph, cs, name):
    k = (ph, cs, name)
    return tuple(EF.loc[k]) if k in EF.index else (np.nan, np.nan)


def coef_panel(ax, ph, rows):
    for cs, dy in zip(COVSETS, np.linspace(-0.27, 0.27, len(COVSETS))):
        col = plt.cm.Greys(SHADE[cs])
        for i, name in enumerate(rows):
            est, se = coef(ph, cs, name)
            if np.isfinite(est):
                ax.errorbar(est, len(rows) - 1 - i + dy, xerr=se if np.isfinite(se) else None,
                            fmt="o", ms=4, color=col, ecolor=col, elinewidth=0.8, mec="0.3", mew=0.3)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows[::-1], fontsize=8)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.axvline(0, color="grey", lw=0.5)
    ax.tick_params(axis="x", labelsize=8)


for ph in ORDER:
    tables = {}
    for cs in COVSETS:
        p = f"{MERGED_DIR}/{ph}__{FIT_TF}__{cs}__{CLASS_TAG}__{BIN_TAG}_merged.full.tsv"
        if os.path.isfile(p):
            t = pd.read_csv(p, sep="\t")
            tables[cs] = t[(t.full_n >= MIN_N) & (t["class"] != "pooled")]
    if not tables:
        continue

    fig = plt.figure(figsize=(15, 6.5))
    gs = GridSpec(2, 2, width_ratios=[3, 1.15], height_ratios=[len(H2_ROWS), len(B2_ROWS) + 0.6],
                  wspace=0.28, hspace=0.35)
    ax, axh, axb = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])

    for cs, t in tables.items():
        for cls in ("other", "FS", "PO"):
            s = t[t["class"] == cls]
            if len(s):
                col = CMAP[cls](SHADE[cs])
                ax.errorbar(s.bin_midpoint, s.full_mean, yerr=s.jk_se, fmt=MARK[cls], ms=MS,
                            color=col, lw=0, elinewidth=0.6, capsize=0, ecolor=col, alpha=0.85,
                            mec="0.3", mew=0.4, zorder=2 if cls == "other" else 4)
        h2u = coef(ph, cs, "h2_Unrel")[0]
        if np.isfinite(h2u):
            xs = np.array(XLIM)
            ax.plot(xs, h2u * xs, "--", lw=1, color=CMAP["other"](SHADE[cs]), alpha=0.8, zorder=3)
        beta_r = coef(ph, cs, "h2_RelOffsets")[0]
        for d, (lo, hi) in HE.DEG_BANDS.items():
            off = coef(ph, cs, f"RelOffsets.{d}")[0]
            if np.isfinite(beta_r) and np.isfinite(off):
                xs = np.array([lo, hi])
                ax.plot(xs, beta_r * xs + off, "-", lw=1.6, color=plt.cm.Greys(SHADE[cs]), zorder=5)

    for lo, _ in HE.DEG_BANDS.values():
        ax.axvline(lo, color="0.85", lw=0.6, zorder=0)
    ax.axhline(0, color="grey", lw=.4)
    ax.set_xlim(*XLIM)
    ax.set_xlabel(r"GRM relatedness  $a_{ij}$")
    ax.set_ylabel(r"mean phenotype cross-product  $\overline{y_i y_j}$")

    coef_panel(axh, ph, H2_ROWS)
    coef_panel(axb, ph, B2_ROWS)
    axh.set_title("estimates (± jackknife SE)", fontsize=9)

    cls_h = [Line2D([], [], lw=0, marker=MARK[c], ms=MS + 1, mfc=CMAP[c](0.75), mec="0.3", label=c)
             for c in ("other", "FS", "PO")]
    fit_h = [Line2D([], [], ls="--", color=CMAP["other"](0.7), label=r"$h^2_{Unrel}\,a$"),
             Line2D([], [], ls="-", lw=1.6, color="0.4", label="RelOffsets fit")]
    cov_h = [Line2D([], [], color=plt.cm.Greys(SHADE[cs]), lw=4, label=cs) for cs in tables]
    leg = ax.legend(handles=cls_h + fit_h, fontsize=8, loc="upper left", framealpha=.9)
    ax.add_artist(leg)
    ax.legend(handles=cov_h, fontsize=8, loc="lower right", title="covariate set",
              title_fontsize=8, framealpha=.9)
    fig.suptitle(f"{ph} — {SAMPLE_SET}, {FIT_TF}, noPO; bins with n$\\geq${MIN_N}; "
                 "grey rules = degree bands",
                 color=CAT_COL[CATEGORY.get(ph, "uncategorised")])
    plt.savefig(f"{OUT}/plots/{ph}_fits_{FIT_TF}.png", dpi=110, bbox_inches="tight")
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
