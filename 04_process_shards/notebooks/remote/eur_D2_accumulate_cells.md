# Binned phenotype cross-products — eur_D2 (single-phenotype test)

Runs `grm_shard_tool accumulate` + `merge` over the 16 eur_D2 Batch shards for
**one** `(phenotype, transform, covariate_set)` combo, to validate the whole
chain end-to-end before committing to the full phenotype list.

Nothing has been accumulated for eur_D2 yet, so this is a clean start — no
stale `.acc.tsv` files to invalidate.

**eur_D2 actuals**: N = 221,992; 24,640,335,028 GRM entries = 91.8 GiB across
16 shards of 5.74 GiB. A numpy scan of one shard took ~59 s, so expect roughly
**10–20 min** for one combo across all 16 shards.

Filenames carry an `EXCL_TAG` so that a later PO-exclusion-aware run writes to
different paths instead of silently reusing these accumulators.

Run cells in order.

---

## Cell 1 — build the tool

```python
import os
import subprocess

REPO_DIR = os.path.expanduser("~/repos/aou-covariance")
TOOL_DIR = f"{REPO_DIR}/GRM-pairs/grm_bin_sharded"
TOOL_BIN = f"{TOOL_DIR}/grm_shard_tool"

subprocess.run(["make"], cwd=TOOL_DIR, check=True)
assert os.path.isfile(TOOL_BIN), f"missing {TOOL_BIN}"
print("tool:", TOOL_BIN)
```

## Cell 2 — paths

```python
WORKSPACE_BUCKET = os.path.expanduser(
    "~/workspace/Data from All of Us Controlled Tier /shared-env-pilot"
)
PROJECT_DIR = "phenotypic_covariance_v9"
SAMPLE_SET = "eur_D2"
N_SHARDS = 16

BUCKET_DIR = f"{WORKSPACE_BUCKET}/{PROJECT_DIR}/03_grm_shards/{SAMPLE_SET}"
SHARDS_DIR = f"{BUCKET_DIR}/shards"
GRM_ID = f"{SHARDS_DIR}/grm_shard_1_of_{N_SHARDS}.grm.id"

PHENO_DIR = (f"{WORKSPACE_BUCKET}/{PROJECT_DIR}/02_phenotype/{SAMPLE_SET}/residualized")

BINS_FILE = f"{REPO_DIR}/GRM-pairs/full_grm_bin/bins.txt"
NBLOCKS = 50
SEED = 1

# No pair exclusion applied yet. Bump this when a PO-exclusion list is wired
# in, so those runs write to separate files rather than reusing these.
EXCL_TAG = "noexcl"

WORK_DIR = os.path.expanduser(f"~/grm_pheno_cov_{SAMPLE_SET}")
PLOTS_DIR = f"{WORK_DIR}/plots"
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

print(f"shards : {SHARDS_DIR}")
print(f"pheno  : {PHENO_DIR}")
print(f"work   : {WORK_DIR}")
```

## Cell 3 — verify inputs

```python
import glob

shard_files = {}
for k in range(1, N_SHARDS + 1):
    p = f"{SHARDS_DIR}/grm_shard_{k}_of_{N_SHARDS}.grm.bin.{k}"
    if os.path.isfile(p):
        shard_files[k] = p

print(f"shards      : {len(shard_files)}/{N_SHARDS}")
print(f"grm.id      : {os.path.isfile(GRM_ID)}")
print(f"bins        : {os.path.isfile(BINS_FILE)}")
print(f"pheno dir   : {os.path.isdir(PHENO_DIR)}")

if os.path.isdir(PHENO_DIR):
    n_pheno = len(glob.glob(f"{PHENO_DIR}/*.pheno"))
    print(f"pheno files : {n_pheno}")
else:
    print("\n*** PHENO_DIR missing — 02_phenotype has not been run for "
          f"SAMPLE_SET={SAMPLE_SET}. Residualize first; this notebook needs "
          "its .pheno output. ***")
```

If `PHENO_DIR` doesn't exist, stop here — eur_D2 is a new sample set and the
residualization (which uses its own keep list and PCs) has to run first.

## Cell 4 — FID alignment

`grm_shard_tool`'s pheno lookup keys on `(FID, IID)`. `.grm.id` has `FID = 0`;
`write_grm_pheno()` emits `FID = IID = person_id`. Rewrite the phenotype
file's FID from the real `.grm.id`, keyed on IID, so it works regardless of
which convention either side uses.

Header lines are `#`-prefixed with exactly two whitespace tokens —
`grm_shard_tool` skips any line that doesn't parse to three, so it ignores them.

```python
import hashlib
import pandas as pd
from datetime import datetime, timezone

grm_ids = pd.read_csv(GRM_ID, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
id_map = dict(zip(grm_ids["IID"], grm_ids["FID"]))
print(f"{len(grm_ids)} GRM ids, {grm_ids['FID'].nunique()} distinct FID -> "
      f"{sorted(grm_ids['FID'].unique())[:3]}")


def build_aligned_pheno(pheno_path, out_path):
    pheno = pd.read_csv(pheno_path, sep=r"\s+", dtype={"FID": str, "IID": str})
    n_source = len(pheno)
    pheno["FID"] = pheno["IID"].map(id_map)
    keep = pheno["FID"].notna()
    out = pheno[keep]
    fingerprint = hashlib.sha256(
        "\n".join(sorted(out["FID"] + "\t" + out["IID"])).encode()
    ).hexdigest()
    with open(out_path, "w") as f:
        f.write(f"# source={os.path.basename(pheno_path)}\n")
        f.write(f"# generated={datetime.now(timezone.utc).isoformat()}\n")
        f.write(f"# n_source={n_source}\n")
        f.write(f"# n_matched={len(out)}\n")
        f.write(f"# n_dropped_unmatched={int((~keep).sum())}\n")
        f.write(f"# id_set_sha256={fingerprint}\n")
        out.to_csv(f, sep=" ", index=False, na_rep="NA")
    return out_path, len(out), int((~keep).sum())
```

## Cell 5 — pick one combo

```python
combos = []
for fname in sorted(os.listdir(PHENO_DIR)):
    if not fname.endswith(".pheno"):
        continue
    parts = fname[:-len(".pheno")].rsplit("__", 2)
    if len(parts) == 3:
        combos.append(tuple(parts) + (fname,))

phenotypes = sorted({c[0] for c in combos})
transforms = sorted({c[1] for c in combos})
covsets = sorted({c[2] for c in combos})

print(f"{len(combos)} files: {len(phenotypes)} phenotypes x "
      f"{len(transforms)} transforms x {len(covsets)} covariate sets")
print(f"\nphenotypes: {phenotypes}")
print(f"transforms: {transforms}")
print(f"covsets   : {covsets}")
```

```python
# <-- set these three from the lists above
PHENOTYPE = "height"
TRANSFORM = "invnorm"
COVSET = "base_pcs"

PHENO_FILE = f"{PHENOTYPE}__{TRANSFORM}__{COVSET}.pheno"
PHENO_PATH = f"{PHENO_DIR}/{PHENO_FILE}"
assert os.path.isfile(PHENO_PATH), f"no such combo: {PHENO_FILE}"

TAG = f"{PHENOTYPE}__{TRANSFORM}__{COVSET}__{EXCL_TAG}"
print(f"testing: {TAG}")
```

`height` with `invnorm` and `base_pcs` is the sensible first test — high
heritability, well behaved, and PCs included so residual ancestry structure
isn't driving the unrelated-region slope. Swap if it isn't in the list.

## Cell 6 — align the phenotype

```python
aligned, n_matched, n_dropped = build_aligned_pheno(
    PHENO_PATH, f"{WORK_DIR}/{TAG}_aligned.pheno")

print(f"{n_matched} matched, {n_dropped} dropped (IID not in .grm.id)")
print(f"-> {aligned}\n")
print(open(aligned).read(500))
```

A large `n_dropped` means the phenotype was residualized against a different
sample set than the GRM was built on — check before continuing.

## Cell 7 — accumulate

One call per shard. Each is a linear scan over 5.74 GiB.

```python
import time

acc_files = []
t0 = time.monotonic()

for k in sorted(shard_files):
    acc_out = f"{WORK_DIR}/{TAG}_shard{k}.acc.tsv"
    if os.path.isfile(acc_out) and os.path.getsize(acc_out) > 0:
        print(f"[{k}/{N_SHARDS}] skip (exists)")
    else:
        t1 = time.monotonic()
        subprocess.run(
            [TOOL_BIN, "accumulate",
             "--grm-id", GRM_ID,
             "--shard", shard_files[k],
             "--parallel", str(k), str(N_SHARDS),
             "--pheno", aligned,
             "--bins", BINS_FILE,
             "--nblocks", str(NBLOCKS),
             "--seed", str(SEED),
             "--out", acc_out],
            check=True,
        )
        print(f"[{k}/{N_SHARDS}] {time.monotonic() - t1:.0f}s")
    acc_files.append(acc_out)

print(f"\n{len(acc_files)} accumulators, {(time.monotonic() - t0) / 60:.1f} min")
```

`accumulate` prints each shard's resolved row range to stderr — check the
first is `[0, 55499)` and the last `[214943, 221992)`.

## Cell 8 — merge

```python
acc_list = f"{WORK_DIR}/{TAG}_acc_list.txt"
with open(acc_list, "w") as f:
    f.write("\n".join(acc_files) + "\n")

out_prefix = f"{WORK_DIR}/{TAG}_merged"
subprocess.run(
    [TOOL_BIN, "merge",
     "--acc-list", acc_list,
     "--bins", BINS_FILE,
     "--nblocks", str(NBLOCKS),
     "--out-prefix", out_prefix],
    check=True,
)

res = pd.read_csv(f"{out_prefix}.full.tsv", sep="\t")
print(f"{len(res)} bins, {int(res['full_n'].sum()):,} pairs binned")
res.head()
```

Total pairs binned will be below 24.6e9 — pairs are dropped when either
phenotype is NaN or `a_ij` falls outside every bin. A big shortfall in the
latter is worth checking: `bins.txt` starts at −0.05, so any more-negative
`a_ij` is silently discarded.

## Cell 9 — inspect

```python
d = res[res["full_n"] > 0].copy()

print("=== unrelated region ===")
print(d[(d.bin_midpoint > -0.01) & (d.bin_midpoint < 0.01)]
      [["bin_midpoint", "full_n", "full_mean", "jk_se"]].to_string(index=False))

print("\n=== related bins (n >= 20) ===")
print(d[(d.bin_midpoint > 0.1) & (d.full_n >= 20)]
      [["bin_midpoint", "full_n", "full_mean", "jk_se"]].to_string(index=False))

print(f"\npairs in 1st-degree band (0.35-0.7): "
      f"{int(d[(d.bin_midpoint >= 0.35) & (d.bin_midpoint < 0.7)]['full_n'].sum()):,}")
print(f"pairs at a ~ 1.0 (0.9-1.1): "
      f"{int(d[(d.bin_midpoint >= 0.9) & (d.bin_midpoint < 1.1)]['full_n'].sum()):,}")
```

The 1st-degree count here should roughly match what the relatedness screen
found — a cross-check that both are reading the same shards the same way.

## Cell 10 — plot

```python
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (lo, hi, title) in zip(axes, [
    (-0.02, 0.02, "unrelated region"),
    (-0.05, 1.10, "full range"),
]):
    s = d[(d.bin_midpoint >= lo) & (d.bin_midpoint <= hi)]
    ax.errorbar(s["bin_midpoint"], s["full_mean"], yerr=s["jk_se"],
                fmt="o", ms=3, lw=0.8, capsize=2)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("relatedness (bin midpoint)")
    ax.set_ylabel("mean phenotype cross-product")
    ax.set_title(f"{TAG}\n{title}")

# weighted slope over the unrelated region only
fit = d[(d.bin_midpoint > -0.01) & (d.bin_midpoint < 0.01) & (d.jk_se > 0)]
slope, intercept = np.polyfit(fit["bin_midpoint"], fit["full_mean"],
                              1, w=1 / fit["jk_se"])
xs = np.array([-0.05, 1.10])
axes[1].plot(xs, intercept + slope * xs, "--", lw=1,
             label=f"unrelated slope = {slope:.3f}")
axes[1].legend()

plt.tight_layout()
plt.savefig(f"{PLOTS_DIR}/{TAG}.png", dpi=120)
plt.show()
print(f"h2_Unrel (unrelated-region slope) ~= {slope:.4f}")
```

## Cell 11 — persist to the bucket

```bash
%%bash -s "$WORK_DIR" "$BUCKET_DIR" "$TAG"
set -e
WORK_DIR=$1; BUCKET_DIR=$2; TAG=$3

DEST="${BUCKET_DIR}/crossproducts"
mkdir -p "$DEST/plots"
cp "${WORK_DIR}/${TAG}_merged".*.tsv "$DEST/"
cp "${WORK_DIR}/plots/${TAG}.png" "$DEST/plots/" 2>/dev/null || true
ls -lh "$DEST"
```

Accumulators stay local — they're per-shard intermediates and cheap to
regenerate. Merged results and plots go to the bucket.

---

## What to check before scaling to all phenotypes

1. **Row ranges** in Cell 7's stderr match `[0, 55499)` … `[214943, 221992)`.
2. **`n_dropped_unmatched`** in Cell 6 is small.
3. **1st-degree pair count** in Cell 9 is consistent with the relatedness screen.
4. **Unrelated-region slope** in Cell 10 is a plausible h² for the phenotype.
5. **Timing** — multiply Cell 7's total by the number of combos to size the
   full run, and use `01b_grm_shard_accumulate_parallel.ipynb`'s
   `ProcessPoolExecutor` pattern if it's too long serially.

## Notes

- `EXCL_TAG` is in every filename. When PO exclusion lands, set it to
  something like `expo` so those runs don't collide with these.
- `bins.txt` spans −0.05 to 1.5. Pairs outside are silently dropped; Cell 8's
  total is the check.
- Accumulators and `.pheno` files are participant-derived — gitignored, never
  commit.
