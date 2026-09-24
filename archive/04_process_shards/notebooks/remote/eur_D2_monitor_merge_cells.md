# eur_D2 — merge and plot whatever has finished

Companion to `eur_D2_accumulate_cells.md`. Safe to run repeatedly **while the
accumulate run is still going**: it only reads accumulators, merges the models
whose 16 shards are all present, and plots what exists.

Read-only with respect to the running job — merged output goes to different
filenames, and merging is idempotent.

Run cells in order; re-run cells 2–5 whenever you want an update.

---

## Cell 1 — config

```python
import os, re, glob, subprocess, time
from collections import defaultdict

REPO = os.path.expanduser("~/repos/AOU-covariance")
SRC  = f"{REPO}/04_process_shards/src/grm_class_tool"
WORK = os.path.expanduser("~/grm_pheno_cov_eur_D2")
BINS = os.path.expanduser("~/bins/bins_wide.txt")
BUCKET_DIR_GS = ("gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
                 "/phenotypic_covariance_v9/03_grm_shards/eur_D2")

N_SHARDS, NBLOCKS = 16, 50
BIN_TAG, CLASS_TAG = "wide", "deg1"

# An accumulator younger than this may still be mid-write. accumulate writes
# all its outputs at the very end of the shard pass, so the window is
# milliseconds, but there is no reason to race it.
SETTLE_SEC = 15

os.makedirs(f"{WORK}/plots", exist_ok=True)
print(WORK)
```

## Cell 2 — progress

```python
ACC_RE = re.compile(r"^(?P<tag>.+)_shard(?P<k>\d+)\.acc\.tsv$")

def scan():
    by_tag = defaultdict(set)
    now = time.time()
    for p in glob.glob(f"{WORK}/*_shard*.acc.tsv"):
        m = ACC_RE.match(os.path.basename(p))
        if not m:
            continue
        st = os.stat(p)
        if st.st_size == 0 or now - st.st_mtime < SETTLE_SEC:
            continue
        by_tag[m.group("tag")].add(int(m.group("k")))
    return by_tag

by_tag = scan()
complete = {t for t, ks in by_tag.items() if len(ks) == N_SHARDS}
partial  = {t: len(ks) for t, ks in by_tag.items() if len(ks) < N_SHARDS}

# phenotype = tag with the trailing __transform__covset__class__bin stripped
pheno_of = lambda t: t.split("__")[0]
phenos = {pheno_of(t) for t in by_tag}

print(f"models seen      : {len(by_tag)}")
print(f"models complete  : {len(complete)}")
print(f"models partial   : {len(partial)}")
print(f"phenotypes touched: {len(phenos)}")
print(f"shards written    : {sum(len(k) for k in by_tag.values())}")
if partial:
    print("\nin flight:")
    for t, n in sorted(partial.items(), key=lambda x: -x[1])[:10]:
        print(f"  {n:>2}/{N_SHARDS}  {t}")
```

## Cell 3 — merge newly-complete models

Skips anything already merged from the same accumulators, so re-running is
cheap.

```python
def needs_merge(tag, shards):
    out = f"{WORK}/{tag}_merged.full.tsv"
    if not os.path.isfile(out):
        return True
    t_out = os.path.getmtime(out)
    return any(os.path.getmtime(f"{WORK}/{tag}_shard{k}.acc.tsv") > t_out for k in shards)

todo = sorted(t for t in complete if needs_merge(t, by_tag[t]))
print(f"{len(todo)} to merge ({len(complete) - len(todo)} already up to date)\n")

os.makedirs(f"{WORK}/lists", exist_ok=True)
merged_now, merge_fail = [], []
for tag in todo:
    accs = [f"{WORK}/{tag}_shard{k}.acc.tsv" for k in sorted(by_tag[tag])]
    lst = f"{WORK}/lists/{tag}_acc.txt"
    with open(lst, "w") as f:
        f.write("\n".join(accs) + "\n")
    r = subprocess.run([SRC, "merge", "--acc-list", lst, "--bins", BINS,
                        "--nblocks", str(NBLOCKS),
                        "--out-prefix", f"{WORK}/{tag}_merged"],
                       capture_output=True, text=True)
    (merged_now if r.returncode == 0 else merge_fail).append(tag)
    if r.returncode != 0:
        print(f"  FAIL {tag}: {r.stderr.strip().splitlines()[-1] if r.stderr else '?'}")

print(f"merged {len(merged_now)}, failed {len(merge_fail)}")
all_merged = sorted(os.path.basename(p)[:-len("_merged.full.tsv")]
                    for p in glob.glob(f"{WORK}/*_merged.full.tsv"))
print(f"{len(all_merged)} merged models available")
```

## Cell 4 — summary table

```python
import pandas as pd, numpy as np

def summarise(tag):
    t = pd.read_csv(f"{WORK}/{tag}_merged.full.tsv", sep="\t")
    t = t[(t.full_n > 0) & (t["class"] != "pooled")]
    f_ = t[(t["class"] == "other") & t.bin_midpoint.between(-0.02, 0.02) & (t.jk_se > 0)]
    if len(f_) < 3:
        return None
    slope, icept = np.polyfit(f_.bin_midpoint, f_.full_mean, 1, w=1 / f_.jk_se)
    parts = tag.split("__")
    rec = {"phenotype": parts[0], "transform": parts[1], "covset": parts[2],
           "h2_unrel": slope, "additive_at_0.5": icept + slope * 0.5}
    for cls in ("PO", "FS"):
        s = t[(t["class"] == cls) & t.bin_midpoint.between(0.35, 0.7)]
        n = int(s.full_n.sum()) if len(s) else 0
        rec[f"{cls}_n"] = n
        rec[f"{cls}_mean"] = s.full_sum.sum() / n if n else np.nan
    rec["n_unrelated"] = int(t[(t["class"] == "other") &
                               t.bin_midpoint.between(-0.02, 0.02)].full_n.sum())
    return rec

rows = [r for r in (summarise(t) for t in all_merged) if r]
S = pd.DataFrame(rows)
S["FS_excess"] = S.FS_mean - S["additive_at_0.5"]
S["PO_excess"] = S.PO_mean - S["additive_at_0.5"]
S = S.sort_values("h2_unrel", ascending=False)
S.to_csv(f"{WORK}/summary_{CLASS_TAG}_{BIN_TAG}.tsv", sep="\t", index=False)

print(f"{len(S)} models; {S['phenotype'].nunique()} phenotypes\n")

# NB: bracket access, not S.transform -- `transform` is a DataFrame method,
# so attribute access returns the method rather than the column.
view = S[S["covset"].eq("base_pcs") & S["transform"].eq("invnorm")]
print(view[["phenotype", "h2_unrel", "PO_n", "PO_excess", "FS_n", "FS_excess"]]
      .round(4).to_string(index=False))
```

## Cell 5 — plot every merged phenotype

One figure per phenotype, one panel per transform, with all four covariate
sets overlaid. Colour encodes the covariate set (light → dark follows the
nested staircase `base` → `+PCs` → `+zip3` → `+SES`); marker encodes the pair
class. Error bars are deliberately faint — with four models overlaid they are
context, not the message.

The comparison to read is whether the PO and FS markers move relative to their
own dashed additive line as covariates are added. If the FS excess shrinks once
zip3 and SES enter, part of what looked like shared sibling environment was
geography and socioeconomic confounding.

```python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MIN_N = 20
COVSETS    = ["base", "base_pcs", "base_pcs_zip3", "base_pcs_zip3_ses"]
TRANSFORMS = ["raw", "invnorm"]
XLIM = (-0.06, 1.12)

CSCOL = {cs: plt.cm.viridis(v) for cs, v in zip(COVSETS, np.linspace(0.15, 0.85, len(COVSETS)))}
MARK  = {"other": ("o", 2.5), "PO": ("^", 7.0), "FS": ("s", 7.0)}

by_pheno = defaultdict(dict)
for tag in all_merged:
    parts = tag.split("__")
    if len(parts) >= 3:
        by_pheno[parts[0]][(parts[1], parts[2])] = tag
print(f"{len(by_pheno)} phenotypes merged\n")


def load(tag):
    t = pd.read_csv(f"{WORK}/{tag}_merged.full.tsv", sep="\t")
    t = t[(t.full_n > 0) & (t["class"] != "pooled")]
    f_ = t[(t["class"] == "other") & t.bin_midpoint.between(-0.02, 0.02) & (t.jk_se > 0)]
    if len(f_) < 3:
        return None, None, None
    slope, icept = np.polyfit(f_.bin_midpoint, f_.full_mean, 1, w=1 / f_.jk_se)
    return t, slope, icept


for pheno in sorted(by_pheno):
    loaded = {k: v for k, v in ((k, load(t)) for k, t in by_pheno[pheno].items())
              if v[0] is not None}
    if not loaded:
        continue

    vals = np.concatenate([d[d.full_n >= MIN_N].full_mean.values for d, _, _ in loaded.values()])
    pad = 0.08 * (np.nanmax(vals) - np.nanmin(vals) or 1)
    ylim = (np.nanmin(vals) - pad, np.nanmax(vals) + pad)

    fig, axes = plt.subplots(1, len(TRANSFORMS), figsize=(17, 6),
                             sharex=True, sharey=True)
    slopes = {}

    for ax, tf in zip(np.atleast_1d(axes), TRANSFORMS):
        for cs in COVSETS:
            got = loaded.get((tf, cs))
            if got is None:
                continue
            t, slope, icept = got
            slopes[(tf, cs)] = slope
            col = CSCOL[cs]
            d = t[t.full_n >= MIN_N]

            for cls, (mk, ms) in MARK.items():
                s = d[d["class"] == cls]
                if not len(s):
                    continue
                ax.errorbar(s.bin_midpoint, s.full_mean, yerr=s.jk_se,
                            fmt=mk, ms=ms, color=col,
                            lw=0, elinewidth=0.6, capsize=0, ecolor=col,
                            alpha=0.45 if cls == "other" else 0.95,
                            mec="none" if cls == "other" else "0.25", mew=0.5,
                            zorder=2 if cls == "other" else 4)

            xs = np.array(XLIM)
            ax.plot(xs, icept + slope * xs, "--", lw=1, color=col, alpha=0.7, zorder=3)

        ax.axhline(0, color="grey", lw=.4)
        ax.axvline(0.5, color="grey", lw=.4, ls=":")
        ax.set_xlim(*XLIM); ax.set_ylim(*ylim)
        ax.set_xlabel(r"GRM relatedness  $a_{ij}$")
        ax.set_title(tf)
    np.atleast_1d(axes)[0].set_ylabel(r"mean phenotype cross-product  $\overline{y_i y_j}$")

    cov_handles = [Line2D([], [], color=CSCOL[cs], lw=3,
                          label=f"{cs}  (slope {slopes.get((TRANSFORMS[-1], cs), float('nan')):.3f})")
                   for cs in COVSETS if any((tf, cs) in slopes for tf in TRANSFORMS)]
    cls_handles = [Line2D([], [], color="0.35", lw=0, marker=mk, ms=ms if cls != "other" else 5,
                          label=cls) for cls, (mk, ms) in MARK.items()]
    leg1 = np.atleast_1d(axes)[0].legend(handles=cov_handles, fontsize=8,
                                         loc="upper left", title="covariate set",
                                         title_fontsize=8, framealpha=.9)
    np.atleast_1d(axes)[0].add_artist(leg1)
    np.atleast_1d(axes)[-1].legend(handles=cls_handles, fontsize=8, loc="upper left",
                                   title="pair class", title_fontsize=8, framealpha=.9)

    fig.suptitle(f"{pheno} — bins with n$\\geq${MIN_N}; dashed = additive prediction "
                 "fitted on each model's own unrelated region", y=1.00)
    plt.tight_layout()
    plt.savefig(f"{WORK}/plots/{pheno}_models.png", dpi=110, bbox_inches="tight")
    plt.show()
```

## Cell 6 — push what exists to the bucket

```bash
%%bash -s "$WORK" "$BUCKET_DIR_GS"
set -e
WORK=$1; DEST="$2/crossproducts"
gcloud storage cp "${WORK}"/*_merged.full.tsv "${DEST}/" 2>/dev/null || true
gcloud storage cp "${WORK}"/*_merged.jk.tsv   "${DEST}/" 2>/dev/null || true
gcloud storage cp "${WORK}"/summary_*.tsv     "${DEST}/" 2>/dev/null || true
gcloud storage cp "${WORK}"/plots/*.png "${DEST}/plots/" 2>/dev/null || true
echo "objects: $(gcloud storage ls "${DEST}/" | wc -l)"
```

Safe to re-run — it overwrites with whatever is current.

---

## Notes

- **`SETTLE_SEC`** skips accumulators written in the last 15 s. `accumulate`
  writes all its outputs at the end of the shard pass, so the write window is
  milliseconds, but a partial read would produce a malformed merge rather than
  an error.
- **Merging is idempotent** and re-merges only when an accumulator is newer
  than the merged file, so leaving this notebook open and re-running cells 2–6
  costs almost nothing.
- Every model is merged independently, so a phenotype whose batch is still
  running simply doesn't appear yet.
