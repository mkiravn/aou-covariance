# eur_D2 — estimators and meta-analysis plots

Applies every relatedness-regression estimator from `notes.md` to each merged
model, with delete-block jackknife SEs, then plots across phenotypes coloured
by trait category.

Estimators are in `04_process_shards/src/he_estimators.py`, tested by
`make test` against synthetic data with known generating values.

| estimator | range | model | readout |
|---|---|---|---|
| `h2_Unrel` | a < 0.02 | y = β a | h² = β |
| `h2_FS` | 0.4–0.6 | y = β a | h² = β (conflates h² and shared env) |
| `h2_PedW25` | 0.05–0.7 | y = β₀ a + β₁ a² | h² = β₀ |
| `h2_Pedf` | 0.05–0.7 | y = β₀ a + class offsets | h² = β₀ |
| `b2_FS` | 0.4–0.6 | y = β₀ a + β₁ | b² = 2β₁ |
| `b2_step` | 0.2–0.6 | y = β₀ a + β_FS + β_HS | b² = 4(β_FS − β_HS) |
| `excess` | 0.4–0.6 | mean(y − h2_Unrel · a) | per class, FS and PO |

Each is computed on three class sets: **noPO** (other + FS, the primary result),
**pooled** (as an unclassified run would give) and, for the FS-range
estimators, **PO** alone as a contrast.

Run after the monitor notebook has merged models. Re-runnable as more finish.

---

## Cell 1 — config

```python
import os, sys, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = os.path.expanduser("~/repos/AOU-covariance")
WORK = os.path.expanduser("~/grm_pheno_cov_eur_D2")
BINS = os.path.expanduser("~/bins/bins_wide.txt")
META = f"{WORK}/meta"
os.makedirs(f"{META}/plots", exist_ok=True)

BIN_TAG, CLASS_TAG, NBLOCKS = "wide", "deg1", 50
COVSETS    = ["base", "base_pcs", "base_pcs_zip3", "base_pcs_zip3_ses"]
TRANSFORMS = ["raw", "invnorm"]
REF_TF, REF_CS = "invnorm", "base_pcs"      # reference model for single-model plots

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
print("ok")
```

## Cell 2 — compute every estimator for every merged model

```python
MID = HE.load_grid(BINS)
tags = sorted(os.path.basename(p)[:-len("_merged.full.tsv")]
              for p in glob.glob(f"{WORK}/*__{CLASS_TAG}__{BIN_TAG}_merged.full.tsv"))
print(f"{len(tags)} merged models")

parts = []
for i, tag in enumerate(tags, 1):
    jk = f"{WORK}/{tag}_merged.jk.tsv"
    if not os.path.isfile(jk):
        print(f"  skip {tag}: no .jk.tsv")
        continue
    S, N = HE.load_arrays(f"{WORK}/{tag}_merged.full.tsv", jk, len(MID), NBLOCKS)
    r = HE.jackknife(S, N, MID, NBLOCKS)
    ph, tf, cs = tag.split("__")[:3]
    r["phenotype"], r["transform"], r["covset"] = ph, tf, cs
    parts.append(r)
    if i % 25 == 0 or i == len(tags):
        print(f"  [{i}/{len(tags)}]")

E = pd.concat(parts, ignore_index=True)
E["category"] = E["phenotype"].map(CATEGORY).fillna("uncategorised")
E.to_csv(f"{META}/estimators_{CLASS_TAG}_{BIN_TAG}.tsv", sep="\t", index=False)

uncat = sorted(E.loc[E["category"] == "uncategorised", "phenotype"].unique())
print(f"\n{E['phenotype'].nunique()} phenotypes, {len(E)} estimate rows"
      + (f"\nuncategorised (add to CATEGORY): {uncat}" if uncat else ""))

(E[(E["transform"] == REF_TF) & (E["covset"] == REF_CS) & (E["classes"].isin(["noPO", "FS"]))]
 .pivot_table(index="phenotype", columns="estimator", values="est")
 [["h2_Unrel", "h2_FS", "h2_PedW25", "h2_Pedf", "b2_FS", "b2_step", "excess"]]
 .round(3))
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


def category_rules(ax):
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
        category_rules(ax)
        ax.set_title(tf)
        ax.set_xlabel(xlabel)
    dress_y(axes[0])
    cs_h = [Line2D([], [], lw=0, marker=CS_MARK[c], ms=6, color="0.4", label=c) for c in COVSETS]
    fig.legend(handles=cat_handles() + cs_h, loc="upper center",
               bbox_to_anchor=(0.5, 0.0), ncol=5, fontsize=8, frameon=False)
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(f"{META}/plots/{fname}.png", dpi=130, bbox_inches="tight")
    plt.show()
```

## Cell 4 — SNP heritability across covariate models

```python
forest(E[(E["estimator"] == "h2_Unrel") & (E["classes"] == "noPO")],
       r"$\hat h^2_{Unrel}$ — slope through the origin over unrelated pairs ($a_{ij}<0.02$)",
       r"$\hat h^2_{Unrel}$  (± jackknife SE)", "h2_unrel_forest")
```

A phenotype whose estimate falls as covariates are added (markers drifting
left along a row) was carrying ancestry, geography or SES structure that the
covariates absorbed.

## Cell 5 — sibling excess and shared environment

```python
forest(E[(E["estimator"] == "excess") & (E["classes"] == "FS")],
       "full-sib excess over the additive prediction  (mean $y_iy_j - \\hat h^2_{Unrel}\\,a_{ij}$, $0.4\\leq a<0.6$)",
       "FS excess", "fs_excess_forest")

forest(E[(E["estimator"] == "excess") & (E["classes"] == "PO")],
       "parent-offspring excess over the additive prediction",
       "PO excess", "po_excess_forest")

forest(E[(E["estimator"] == "b2_FS") & (E["classes"] == "noPO")],
       r"$\hat b^2_{FS}$ — shared environment among first-degree non-PO pairs",
       r"$\hat b^2_{FS}$", "b2_fs_forest")

# FS against PO excess, reference model
R = E[(E["transform"] == REF_TF) & (E["covset"] == REF_CS) & (E["estimator"] == "excess")]
fs = R[R["classes"] == "FS"].set_index("phenotype")
po = R[R["classes"] == "PO"].set_index("phenotype")
J = fs[["est", "se", "category"]].join(po[["est", "se"]], lsuffix="_fs", rsuffix="_po")
J = J.dropna(subset=["est_fs", "est_po"])

fig, ax = plt.subplots(figsize=(8, 8))
for ph, row in J.iterrows():
    c = CAT_COL[row["category"]]
    ax.errorbar(row["est_po"], row["est_fs"], xerr=row["se_po"], yerr=row["se_fs"],
                fmt="o", ms=6, color=c, ecolor=c, elinewidth=0.7, alpha=0.9,
                mec="0.25", mew=0.4)
    ax.annotate(ph, (row["est_po"], row["est_fs"]), fontsize=6, color=c,
                xytext=(3, 3), textcoords="offset points")
ext = np.r_[J["est_fs"] - J["se_fs"], J["est_fs"] + J["se_fs"],
            J["est_po"] - J["se_po"], J["est_po"] + J["se_po"], 0.0]
lo, hi = np.nanmin(ext), np.nanmax(ext)
pad = 0.06 * (hi - lo or 1)
lo, hi = lo - pad, hi + pad
ax.plot([lo, hi], [lo, hi], "k--", lw=1)
ax.axhline(0, color="grey", lw=.5); ax.axvline(0, color="grey", lw=.5)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("PO excess over additive"); ax.set_ylabel("FS excess over additive")
ax.set_title(f"FS vs PO excess — {REF_TF}, {REF_CS}\n"
             "above the diagonal: dominance + shared sibling environment")
ax.legend(handles=cat_handles(), fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig(f"{META}/plots/fs_vs_po_excess.png", dpi=130, bbox_inches="tight")
plt.show()
```

FS and PO share expected additive relatedness, so under a purely additive model
both excesses are zero. FS carry 0.25 dominance relatedness and a sibling
environment; PO carry neither — points above the diagonal are that difference.

## Cell 6 — all estimators, reference model

Filled = noPO (primary), hollow = pooled (what an unclassified run gives),
triangle = PO pairs alone. The gap between filled and hollow is the effect of
separating parent-offspring pairs.

```python
PANELS = [("h2_Unrel", r"$h^2_{Unrel}$"), ("h2_FS", r"$h^2_{FS}$"),
          ("h2_PedW25", r"$h^2_{Ped,W25}$"), ("h2_Pedf", r"$h^2_{Ped*f}$"),
          ("b2_FS", r"$b^2_{FS}$"), ("b2_step", r"$b^2_{step}$")]
STYLE = {"noPO": (0.0, "o", True), "pooled": (0.22, "o", False), "PO": (-0.22, "^", True)}

R = E[(E["transform"] == REF_TF) & (E["covset"] == REF_CS)]
fig, axes = plt.subplots(1, len(PANELS), figsize=(22, 0.3 * len(ORDER) + 2), sharey=True)
for ax, (est, label) in zip(axes, PANELS):
    for cl, (off, mk, filled) in STYLE.items():
        s = R[(R["estimator"] == est) & (R["classes"] == cl)
              & R["phenotype"].isin(YPOS) & np.isfinite(R["est"])]
        if s.empty:
            continue
        y = s["phenotype"].map(YPOS) + off
        cols = s["category"].map(CAT_COL).tolist()
        se = s["se"].fillna(0)
        ax.hlines(y, s["est"] - se, s["est"] + se, colors=cols, lw=0.8,
                  alpha=0.35 if cl == "PO" else 0.55)
        ax.scatter(s["est"], y, marker=mk, s=26, zorder=3, linewidths=0.8,
                   c=cols if filled else "none", edgecolors=cols if not filled else "0.25",
                   alpha=0.55 if cl == "PO" else 1.0)
    ax.axvline(0, color="grey", lw=0.6)
    category_rules(ax)
    ax.set_title(label)
dress_y(axes[0])
cls_h = [Line2D([], [], lw=0, marker="o", ms=6, mfc="0.4", mec="0.25", label="noPO"),
         Line2D([], [], lw=0, marker="o", ms=6, mfc="none", mec="0.4", label="pooled"),
         Line2D([], [], lw=0, marker="^", ms=6, mfc="0.6", mec="0.25", label="PO only")]
fig.legend(handles=cat_handles() + cls_h, loc="upper center", bbox_to_anchor=(0.5, 0.0),
           ncol=5, fontsize=8, frameon=False)
fig.suptitle(f"all estimators — {REF_TF}, {REF_CS}")
plt.tight_layout()
plt.savefig(f"{META}/plots/all_estimators_{REF_TF}_{REF_CS}.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 7 — estimator agreement with h2_Unrel

```python
R = E[(E["transform"] == REF_TF) & (E["classes"] == "noPO")]
base = R[R["estimator"] == "h2_Unrel"].set_index(["phenotype", "covset"])

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
for ax, est in zip(axes, ["h2_FS", "h2_PedW25", "h2_Pedf"]):
    other = R[R["estimator"] == est].set_index(["phenotype", "covset"])
    J = base[["est", "se", "category"]].join(other[["est", "se"]], lsuffix="_u", rsuffix="_e").dropna()
    for cs in COVSETS:
        s = J.xs(cs, level="covset") if cs in J.index.get_level_values("covset") else None
        if s is None or s.empty:
            continue
        cols = s["category"].map(CAT_COL).tolist()
        ax.errorbar(s["est_u"], s["est_e"], xerr=s["se_u"], yerr=s["se_e"], fmt="none",
                    ecolor="0.75", elinewidth=0.5, zorder=1)
        ax.scatter(s["est_u"], s["est_e"], c=cols, marker=CS_MARK[cs], s=26,
                   edgecolors="0.25", linewidths=0.4, zorder=3)
    vals = J[["est_u", "est_e"]].to_numpy()
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_xlabel(r"$\hat h^2_{Unrel}$"); ax.set_ylabel(est); ax.set_title(f"{est} vs h2_Unrel")
cs_h = [Line2D([], [], lw=0, marker=CS_MARK[c], ms=6, color="0.4", label=c) for c in COVSETS]
fig.legend(handles=cat_handles() + cs_h, loc="upper center", bbox_to_anchor=(0.5, 0.0),
           ncol=5, fontsize=8, frameon=False)
fig.suptitle(f"pedigree-based h² against SNP h² — {REF_TF}, noPO")
plt.tight_layout()
plt.savefig(f"{META}/plots/estimator_agreement_{REF_TF}.png", dpi=130, bbox_inches="tight")
plt.show()
```

`h2_FS` sitting above the identity line is expected — its slope through the
origin absorbs shared environment. `h2_PedW25` and `h2_Pedf` exist to remove
that, so they should sit closer to the line.

## Cell 8 — copy everything to the bucket

`gcloud storage cp` does not fail loudly on a partial object, so every file is
checked against its remote size afterwards. Re-run to retry anything flagged.

```python
import subprocess

DEST = ("gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
        "/phenotypic_covariance_v9/03_grm_shards/eur_D2/crossproducts")
COPY_ACCUMULATORS = False     # ~4,600 files, several GB; only needed to re-merge

GROUPS = {
    "merged":     sorted(glob.glob(f"{WORK}/*_merged.full.tsv") + glob.glob(f"{WORK}/*_merged.jk.tsv")),
    "summaries":  sorted(glob.glob(f"{WORK}/summary_*.tsv")),
    "plots":      sorted(glob.glob(f"{WORK}/plots/*.png")),
    "meta":       sorted(glob.glob(f"{META}/*.tsv")),
    "meta/plots": sorted(glob.glob(f"{META}/plots/*.png")),
}
if COPY_ACCUMULATORS:
    GROUPS["accumulators"] = sorted(glob.glob(f"{WORK}/*_shard*.acc.tsv"))


def remote_sizes(prefix):
    out = subprocess.run(["gcloud", "storage", "ls", "-l", f"{prefix}/"],
                         capture_output=True, text=True).stdout
    sizes = {}
    for line in out.splitlines():
        f = line.split()
        if len(f) >= 3 and f[0].isdigit() and f[-1].startswith("gs://"):
            sizes[f[-1].rsplit("/", 1)[1]] = int(f[0])
    return sizes


for sub, files in GROUPS.items():
    if not files:
        print(f"{sub:<13} nothing to copy")
        continue
    dst = f"{DEST}/{sub}"
    subprocess.run(["gcloud", "storage", "cp", "-I", f"{dst}/"],
                   input="\n".join(files), text=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    remote = remote_sizes(dst)
    bad = [os.path.basename(f) for f in files
           if remote.get(os.path.basename(f)) != os.path.getsize(f)]
    print(f"{sub:<13} {len(files) - len(bad):>5}/{len(files)} verified"
          + (f"   MISMATCH: {bad[:5]}{' ...' if len(bad) > 5 else ''}" if bad else ""))
```

Layout in the bucket: `crossproducts/{merged, summaries, plots, meta, meta/plots}`.
Files an earlier monitor run pushed directly under `crossproducts/` are left
where they are.

---

## Caveats

- **Pair-weighted fits on bin means** use each bin's midpoint as `a_ij`. With
  0.001-wide bins below 0.02 this is negligible for `h2_Unrel`; in the 0.02-wide
  related bins it is a small approximation to pair-level OLS.
- **PO-only fits are noisy** — a few thousand pairs spread over 0.4–0.6. Read
  them as a contrast against FS, not as standalone estimates.
- **`excess` assumes the additive prediction extrapolates linearly** from the
  unrelated region to a = 0.5.
- **`b2_step_ratio`** is in the estimator table for the geometric-decay test
  but not plotted; it is unstable whenever the HS intercept is near zero.
