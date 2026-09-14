# Binned phenotype cross-products — eur_D2

Runs `grm_class_tool` over the 16 eur_D2 shards for every residualized
phenotype, keeping PO and FS pairs in their own bin sets.

Batched by phenotype: one invocation per (phenotype, shard) covers that
phenotype's 8 models (2 transforms × 4 covariate sets) in a single pass, since
the 4-byte read, the bin lookup and the pair-class cursor are shared across
them. 576 invocations for all 288 models, ~2 h at 32-way concurrency.

**eur_D2 actuals**: N = 221,992; 24,640,335,028 entries = 91.8 GiB across 16
shards of 5.74 GiB; 8,887 classified 1st-degree pairs; 153 bins (`--wide`).

No FID alignment step is needed — `grm_class_tool` keys on IID alone, so it
reads the residualized `.pheno` files straight from the bucket.

Run cells in order.

---

## Cell 1 — pull, build, verify the binary

```python
import os, subprocess

REPO = os.path.expanduser("~/repos/AOU-covariance")
SRC  = f"{REPO}/04_process_shards/src/grm_class_tool"

subprocess.run(["git", "pull"], cwd=REPO, check=True)
subprocess.run(["make"], cwd=f"{REPO}/04_process_shards/src", check=True)

probe = subprocess.run([SRC, "accumulate"], capture_output=True, text=True)
assert "--pheno-list" in probe.stderr, \
    "stale binary -- multi-phenotype support missing; check git pull && make"
print("binary OK\n")
print(probe.stderr.strip())
```

A stale binary fails in argument parsing after 0.2 s, 576 times over. Check it
once instead.

## Cell 2 — config

```python
SHARDS  = os.path.expanduser("~/scratch_grm/shards")
CLASSES = os.path.expanduser("~/scratch_grm/relatedness_screen/deg1_classified.tsv")
BINS    = os.path.expanduser("~/bins/bins_wide.txt")
WORK    = os.path.expanduser("~/grm_pheno_cov_eur_D2")
PHENO_DIR = os.path.expanduser(
    "~/workspace/Data from All of Us Controlled Tier /shared-env-pilot"
    "/phenotypic_covariance_v9/02_phenotype/eur_D2/residualized")
BUCKET_DIR_GS = ("gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
                 "/phenotypic_covariance_v9/03_grm_shards/eur_D2")

N_IDS, N_SHARDS = 221992, 16
N_CONCURRENT, NBLOCKS, SEED = 32, 50, 1

# In every accumulator filename: the accumulator stores bin INDICES, not edges,
# so merging a 153-bin accumulator against a 236-bin file would silently
# attribute sums to the wrong relatedness. Same for the classification.
BIN_TAG, CLASS_TAG = "wide", "deg1"

GRM_ID = f"{SHARDS}/grm_shard_1_of_{N_SHARDS}.grm.id"
for d in ("logs", "lists", "plots"):
    os.makedirs(f"{WORK}/{d}", exist_ok=True)

print(f"shards {SHARDS}\nwork   {WORK}\ntags   {CLASS_TAG}/{BIN_TAG}")
```

## Cell 3 — verify the shards

`gcloud storage cp` does not fail loudly on a partial object, and a truncated
shard is not detectable by eye — sizes differ by a few hundred MB out of 5.8 GB.
GRM shard sizes are exactly computable, so check rather than assume.

```python
import math

def _tdo(t):
    if t == 0: return 1
    v = int(math.sqrt(t)) + 2
    while v > 1 and (v-1)*(v-2) >= t: v -= 1
    while v*(v-1) < t: v += 1
    return v

def row_start(idx):                       # idx is 0-based
    v = _tdo((N_IDS * (N_IDS - 1) * idx) // N_SHARDS)
    return 0 if v == 1 else v

C = lambda i: i * (i + 1) // 2

bad = []
print(f"{'shard':>5}  {'rows':>20}  {'expected':>15}  {'actual':>15}")
for k in range(1, N_SHARDS + 1):
    a = row_start(k - 1)
    b = N_IDS if k == N_SHARDS else row_start(k)
    exp = (C(b) - C(a)) * 4
    p = f"{SHARDS}/grm_shard_{k}_of_{N_SHARDS}.grm.bin.{k}"
    act = os.path.getsize(p) if os.path.isfile(p) else -1
    if act != exp:
        bad.append(k)
    print(f"{k:>5}  [{a:>6}, {b:>6})  {exp:>15,}  {act:>15,}  "
          f"{'OK' if act == exp else 'BAD %+d' % (act - exp)}")

for f in (GRM_ID, BINS, CLASSES):
    print(f"{'OK  ' if os.path.isfile(f) else 'MISS'} {f}")

assert not bad, f"re-copy these shards: {bad}"
print("\nall shards verified")
```

If any are BAD, re-copy only those:

```bash
%%bash
GS="gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9/03_grm_shards/eur_D2/shards"
DEST=~/scratch_grm/shards
for k in 1 6; do          # <-- the bad shard numbers
  rm -f "$DEST/grm_shard_${k}_of_16.grm.bin.${k}"
  gcloud storage cp "$GS/grm_shard_${k}_of_16.grm.bin.${k}" "$DEST/"
done
```

## Cell 4 — discover phenotypes and write per-phenotype lists

```python
from collections import defaultdict

models = defaultdict(list)
for f in sorted(os.listdir(PHENO_DIR)):
    if f.endswith(".pheno"):
        stem = f[:-len(".pheno")]
        models[stem.rsplit("__", 2)[0]].append(stem)
models = dict(sorted(models.items()))

n_models = sum(len(v) for v in models.values())
print(f"{len(models)} phenotypes, {n_models} models "
      f"({n_models // max(1, len(models))} per phenotype)\n")
for p in list(models)[:5]:
    print(f"  {p}: {[m.split('__', 1)[1] for m in models[p]]}")

for pheno, combos in models.items():
    with open(f"{WORK}/lists/{pheno}.tsv", "w") as f:
        for c in combos:
            f.write(f"{c}__{CLASS_TAG}__{BIN_TAG}\t{PHENO_DIR}/{c}.pheno\n")
print(f"\nwrote {len(models)} list files to {WORK}/lists/")
```

To run a subset, filter `models` here — e.g. keep only `invnorm__base_pcs`:

```python
# models = {p: [m for m in ms if m.endswith("__invnorm__base_pcs")]
#           for p, ms in models.items()}
# models = {p: ms for p, ms in models.items() if ms}
```

## Cell 5 — accumulate

```python
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def acc_path(combo, k):
    return f"{WORK}/{combo}__{CLASS_TAG}__{BIN_TAG}_shard{k}.acc.tsv"

def run(pheno, k):
    outs = [acc_path(c, k) for c in models[pheno]]
    if all(os.path.isfile(p) and os.path.getsize(p) > 0 for p in outs):
        return pheno, k, 0.0, "skip", ""

    log = f"{WORK}/logs/{pheno}_shard{k}.log"
    t0 = time.monotonic()
    with open(log, "w") as lf:
        r = subprocess.run(
            [SRC, "accumulate",
             "--grm-id", GRM_ID,
             "--shard", f"{SHARDS}/grm_shard_{k}_of_{N_SHARDS}.grm.bin.{k}",
             "--parallel", str(k), str(N_SHARDS),
             "--pheno-list", f"{WORK}/lists/{pheno}.tsv",
             "--bins", BINS,
             "--pair-classes", CLASSES,
             "--nblocks", str(NBLOCKS), "--seed", str(SEED),
             # {name} must reach the tool literally -- not an f-string
             "--out-pattern", WORK + "/{name}_shard" + str(k) + ".acc.tsv"],
            stdout=lf, stderr=subprocess.STDOUT, text=True)
    el = time.monotonic() - t0

    if r.returncode != 0:
        for p in outs:                     # drop partials so a rerun redoes them
            if os.path.isfile(p):
                os.remove(p)
        return pheno, k, el, f"FAIL rc={r.returncode}", log
    return pheno, k, el, "ok", ""


tasks = [(p, k) for p in models for k in range(1, N_SHARDS + 1)]
print(f"{len(tasks)} invocations, {N_CONCURRENT} concurrent\n")

fails, times, t0 = [], [], time.monotonic()
with ThreadPoolExecutor(max_workers=N_CONCURRENT) as pool:
    futs = [pool.submit(run, p, k) for p, k in tasks]
    for i, f in enumerate(as_completed(futs), 1):
        pheno, k, el, status, info = f.result()
        if status.startswith("FAIL"):
            fails.append((pheno, k, info))
        elif status == "ok":
            times.append(el)
        if i % 25 == 0 or status.startswith("FAIL") or i == len(tasks):
            eta = (len(tasks) - i) * (sum(times)/len(times) if times else 0) / N_CONCURRENT
            print(f"[{i:>4}/{len(tasks)}] {pheno} shard {k:>2}: {el:6.1f}s {status} {info}"
                  f"   ~{eta/60:.0f} min left")

print(f"\ndone in {(time.monotonic()-t0)/60:.1f} min, {len(fails)} failed")
for p, k, log in fails:
    print(f"  FAILED {p} shard {k}: {log}")
```

Progress goes to the per-invocation logs, not the notebook — `tail -f` one to
watch a shard advance.

## Cell 6 — merge

```python
import glob

merged, incomplete = [], []
for pheno, combos in models.items():
    for c in combos:
        tag = f"{c}__{CLASS_TAG}__{BIN_TAG}"
        accs = sorted(glob.glob(f"{WORK}/{tag}_shard*.acc.tsv"),
                      key=lambda p: int(p.rsplit("shard", 1)[1].split(".")[0]))
        if len(accs) != N_SHARDS:
            incomplete.append((tag, len(accs)))
            continue
        lst = f"{WORK}/lists/{tag}_acc.txt"
        with open(lst, "w") as f:
            f.write("\n".join(accs) + "\n")
        subprocess.run([SRC, "merge", "--acc-list", lst, "--bins", BINS,
                        "--nblocks", str(NBLOCKS),
                        "--out-prefix", f"{WORK}/{tag}_merged"],
                       check=True, capture_output=True)
        merged.append(tag)

print(f"{len(merged)} merged")
for tag, n in incomplete:
    print(f"  INCOMPLETE {tag}: {n}/{N_SHARDS} shards")
```

## Cell 7 — plot one

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt

TAG = f"height__invnorm__base_pcs__{CLASS_TAG}__{BIN_TAG}"   # <-- pick any merged tag
MIN_N = 20

raw = pd.read_csv(f"{WORK}/{TAG}_merged.full.tsv", sep="\t")
raw = raw[(raw.full_n > 0) & (raw["class"] != "pooled")]
d = raw[raw.full_n >= MIN_N]

COL = {"other": "#4C78A8", "FS": "#F58518", "PO": "#E45756"}
fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

ax = axes[0]
for cls, c in COL.items():
    s = d[d["class"] == cls]
    if len(s):
        ax.errorbar(s.bin_midpoint, s.full_mean, yerr=s.jk_se, fmt="o",
                    ms=3.5 if cls == "other" else 8, lw=.8, capsize=2, color=c,
                    zorder=1 if cls == "other" else 3,
                    label=f"{cls} ({int(s.full_n.sum()):,})")
fit = d[(d["class"] == "other") & d.bin_midpoint.between(-0.02, 0.02) & (d.jk_se > 0)]
slope, icept = np.polyfit(fit.bin_midpoint, fit.full_mean, 1, w=1/fit.jk_se)
xs = np.array([raw.bin_midpoint.min(), raw.bin_midpoint.max()])
ax.plot(xs, icept + slope*xs, "k--", lw=1.2, zorder=2,
        label=f"additive (slope={slope:.3f})")
ax.axhline(0, color="grey", lw=.5)
ax.set_xlabel(r"$a_{ij}$"); ax.set_ylabel(r"mean $y_i y_j$")
ax.set_title(f"cross-product by class (n$\\geq${MIN_N})"); ax.legend(fontsize=8)

ax = axes[1]
for cls, c in COL.items():
    s = raw[raw["class"] == cls].sort_values("bin_midpoint")
    if len(s):
        ax.step(s.bin_midpoint, s.full_n, where="mid", color=c, lw=1.4, label=cls)
ax.set_yscale("log")
ax.set_xlabel(r"$a_{ij}$"); ax.set_ylabel("pairs per bin (log)")
ax.set_title("relatedness distribution"); ax.legend(fontsize=8)

ax = axes[2]
for cls, c in COL.items():
    s = raw[(raw["class"] == cls) & raw.bin_midpoint.between(0.2, 1.15)].sort_values("bin_midpoint")
    if len(s) and s.full_n.sum() > 0:
        ax.step(s.bin_midpoint, s.full_n / s.full_n.sum(), where="mid",
                color=c, lw=1.6, label=f"{cls} ({int(s.full_n.sum()):,})")
ax.axvline(0.5, color="grey", lw=.5, ls=":")
ax.set_xlabel(r"$a_{ij}$"); ax.set_ylabel("fraction of class")
ax.set_title("shape in the related region"); ax.legend(fontsize=8)

plt.suptitle(TAG, y=1.02)
plt.tight_layout(); plt.savefig(f"{WORK}/plots/{TAG}.png", dpi=120, bbox_inches="tight")
plt.show()

print(f"h2_Unrel = {slope:.4f}\n")
print("1st-degree band (0.35-0.7), all bins:")
for cls in ["PO", "FS", "other"]:
    s = raw[(raw["class"] == cls) & raw.bin_midpoint.between(0.35, 0.7)]
    if len(s):
        n = s.full_n.sum(); m = s.full_sum.sum()/n
        print(f"  {cls:<6} n={int(n):>7}  mean={m:.4f}  "
              f"excess over additive {m - (icept + slope*0.5):+.4f}")
```

## Cell 8 — summary across phenotypes

```python
rows = []
for tag in merged:
    t = pd.read_csv(f"{WORK}/{tag}_merged.full.tsv", sep="\t")
    t = t[(t.full_n > 0) & (t["class"] != "pooled")]
    f_ = t[(t["class"] == "other") & t.bin_midpoint.between(-0.02, 0.02) & (t.jk_se > 0)]
    if len(f_) < 3:
        continue
    sl, ic = np.polyfit(f_.bin_midpoint, f_.full_mean, 1, w=1/f_.jk_se)
    rec = {"tag": tag, "h2_unrel": sl}
    for cls in ["PO", "FS"]:
        s = t[(t["class"] == cls) & t.bin_midpoint.between(0.35, 0.7)]
        rec[f"{cls}_n"] = int(s.full_n.sum()) if len(s) else 0
        rec[f"{cls}_mean"] = s.full_sum.sum()/s.full_n.sum() if len(s) and s.full_n.sum() else np.nan
    rec["additive_at_0.5"] = ic + sl*0.5
    rows.append(rec)

summary = pd.DataFrame(rows).sort_values("h2_unrel", ascending=False)
summary["FS_excess"] = summary["FS_mean"] - summary["additive_at_0.5"]
summary["PO_excess"] = summary["PO_mean"] - summary["additive_at_0.5"]
summary.to_csv(f"{WORK}/summary_{CLASS_TAG}_{BIN_TAG}.tsv", sep="\t", index=False)
summary.head(30)
```

## Cell 9 — persist to the bucket

```bash
%%bash -s "$WORK" "$BUCKET_DIR_GS" "$CLASS_TAG" "$BIN_TAG"
set -e
WORK=$1; BUCKET=$2; CLASS_TAG=$3; BIN_TAG=$4
DEST="${BUCKET}/crossproducts"

gcloud storage cp "${WORK}"/*_merged.full.tsv "${DEST}/" 
gcloud storage cp "${WORK}"/*_merged.jk.tsv   "${DEST}/"
gcloud storage cp "${WORK}/summary_${CLASS_TAG}_${BIN_TAG}.tsv" "${DEST}/"
gcloud storage cp "${WORK}"/plots/*.png "${DEST}/plots/" 2>/dev/null || true
gcloud storage ls "${DEST}/" | wc -l
```

Accumulators stay local — there are `n_models × 16` of them and they're cheap
to regenerate. Merged results, the summary and plots go to the bucket.

---

## Checks before trusting the output

1. **Cell 3** passes — a truncated shard is the one failure that produces
   plausible-looking wrong numbers rather than an error.
2. **`pair classes: 8887 read, 8887 mapped, 0 unmapped`** in any log. Zero
   mapped looks identical to "no classed pairs existed".
3. **PO and FS both non-zero** in the per-phenotype class counts.
4. **`h2_Unrel`** for height in a plausible range (~0.4–0.7 for a SNP-based
   estimate in a European-ancestry sample). If height is wrong, nothing else
   is interpretable.
5. **FS excess > PO excess** over the additive line — FS carry 0.25 dominance
   relatedness and a shared sibling environment; PO have neither.

## Notes

- `EXCL_TAG`/`BIN_TAG` in every filename: bins are stored as indices, so
  mixing schemes corrupts silently rather than erroring.
- The skip check requires all of a phenotype's models to exist, so an
  interrupted run redoes that phenotype's batch rather than leaving a partial
  set.
- Accumulators and `.pheno` files are participant-derived — gitignored, never
  commit.
