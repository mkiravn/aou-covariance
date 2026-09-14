# Relatedness screen — candidate individuals for KING

Scans the GRM shards for pairs above a relatedness threshold and writes the
union of individuals involved, as a `--keep` list for a much smaller KING /
`plink2 --make-king-table` run.

Running KING on the full cohort is O(N²) over 2.46e10 pairs and doesn't
finish. Roughly 10-12% of AoU participants have a 1st- or 2nd-degree relative,
so screening first drops N from ~222K to maybe ~25K — about 1.3% of the pair
space.

The GRM can't separate PO from FS (both a_ij ≈ 0.5), which is what KING is
still needed for. It can cheaply say *which pairs are worth asking about*.

**Threshold**: `0.2`. The screen only has to avoid false negatives — a missed
PO pair stays in the FS bin and contaminates exactly what this is meant to
clean. Extra individuals cost only a slightly larger KING run. 0.2 also
captures most 2nd-degree pairs, useful if the a_ij ≈ 0.25 bin needs the same
treatment later.

**eur_D2 actuals (2026-09-12)**: N = 221,992; 24,640,335,028 GRM entries =
91.8 GiB across 16 Batch shards of 5.74 GiB each.

Run cells in order. **Cell 4 tests one shard** — check it before the full scan.

---

## Cell 1 — config

```python
import os
import glob

SHARD_DIR = os.path.expanduser(
    "~/workspace/Data from All of Us Controlled Tier /shared-env-pilot"
    "/phenotypic_covariance_v9/03_grm_shards/eur_D2/shards"
)

# Must match the run that produced the shards: 16 for the Batch run,
# 300 for the local 03a run. Row-range recovery depends on it.
N_SHARDS = 16

GRM_ID_PATH = os.path.join(SHARD_DIR, f"grm_shard_1_of_{N_SHARDS}.grm.id")

THRESHOLD = 0.2

OUT_DIR = os.path.expanduser("~/scratch_grm/relatedness_screen")
os.makedirs(OUT_DIR, exist_ok=True)

def shard_path(k):
    return os.path.join(SHARD_DIR, f"grm_shard_{k}_of_{N_SHARDS}.grm.bin.{k}")

missing = [k for k in range(1, N_SHARDS + 1) if not os.path.isfile(shard_path(k))]
total_gb = sum(os.path.getsize(shard_path(k)) for k in range(1, N_SHARDS + 1)
               if os.path.isfile(shard_path(k))) / 1024**3

print(f"shard dir : {SHARD_DIR}")
print(f"n_shards  : {N_SHARDS}  ({N_SHARDS - len(missing)} present, {total_gb:.1f} GiB)")
print(f"threshold : {THRESHOLD}")
if missing:
    print(f"MISSING   : {missing}")
```

## Cell 2 — load IDs

Row index `i` in the GRM corresponds to line `i` of `.grm.id`.

```python
ids = []
with open(GRM_ID_PATH) as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            ids.append((parts[0], parts[1]))

N_IDS = len(ids)
fids = {f for f, _ in ids}

print(f"{N_IDS} individuals")
print(f"expected GRM entries: {N_IDS * (N_IDS + 1) // 2:,} "
      f"({N_IDS * (N_IDS + 1) // 2 * 4 / 1024**3:.1f} GiB)")
print(f"first 3 rows: {ids[:3]}")
print(f"distinct FIDs: {len(fids)} -> {sorted(fids)[:5]}")
```

Confirmed 2026-09-12: **FID is `0` for everyone**, IID is `person_id`. So
`--keep` lists must be written as `0<TAB><person_id>` (matching the panel
`.fam`), and any pair/phenotype join should key on IID alone.

> This is also why `grm_shard_tool`'s phenotype join needs fixing before
> accumulate: `read_pheno_aligned()` keys on `(FID, IID)` but
> `write_grm_pheno()` emits `FID = IID = person_id`, so every phenotype
> would read as NaN. Fix in `02_phenotype/scripts/local/residualize_lib.R`
> (`transmute(FID = 0, ...)`) or key on IID in the tool.

## Cell 3 — row-range recovery

Port of `grm_shard_tool.cpp`'s `plink_parallel_row_start` /
`resolve_shard_range`. plink 1.9 picks the `--parallel` split by balancing
off-diagonal pair counts (`v*(v-1) >= target`), then the file stores each
row's off-diagonals followed by one diagonal entry.

```python
import math

def triangle_divide_off_diag(target):
    if target == 0:
        return 1
    v = int(math.sqrt(target) + 1.0) + 1
    if v < 1:
        v = 1
    while v > 1 and (v - 1) * (v - 2) >= target:
        v -= 1
    while v * (v - 1) < target:
        v += 1
    return v

def plink_parallel_row_start(n_ids, parallel_idx, n_shards):
    """parallel_idx is 0-based."""
    ct_tot = n_ids * (n_ids - 1)
    target = (ct_tot * parallel_idx) // n_shards
    v = triangle_divide_off_diag(target)
    return 0 if v == 1 else v

def try_row_start(candidate_start, n_ids, shard_floats):
    consumed = 0
    i = candidate_start
    while consumed < shard_floats and i < n_ids:
        consumed += i + 1
        i += 1
    return i if consumed == shard_floats else None

def resolve_shard_range(n_ids, k, n_shards, shard_floats):
    """Returns (row_start, row_end_exclusive), verified against file size."""
    guess = plink_parallel_row_start(n_ids, k - 1, n_shards)
    for off in (0, -1, 1, -2, 2, -3, 3):
        cand = guess + off
        if cand < 0 or cand >= n_ids:
            continue
        row_end = try_row_start(cand, n_ids, shard_floats)
        if row_end is not None:
            return cand, row_end
    raise RuntimeError(
        f"could not resolve row range: n_ids={n_ids} parallel={k}/{n_shards} "
        f"shard_floats={shard_floats} guess={guess}"
    )

def shard_range_from_file(path, k):
    nbytes = os.path.getsize(path)
    if nbytes % 4 != 0:
        raise RuntimeError(f"{path}: size not a multiple of 4")
    return resolve_shard_range(N_IDS, k, N_SHARDS, nbytes // 4)

print("expected ranges:")
for k in (1, 2, N_SHARDS):
    print(f"  shard {k:>2}: {shard_range_from_file(shard_path(k), k)}")
```

For eur_D2 shard 1 should be `(0, 55499)` and shard 16 `(214943, 221992)`.

## Cell 4 — test on one shard

Shard 16 is the better test: fewest rows but each is ~222K floats, so it
exercises the large-block read path over the densest part of the triangle.

```python
import numpy as np
import time

def scan_shard(path, k, threshold, row_block=256):
    """Returns (row_start, row_end, [(i, j, a_ij), ...]) above threshold."""
    row_start, row_end = shard_range_from_file(path, k)
    hits = []
    with open(path, "rb") as f:
        i = row_start
        while i < row_end:
            stop = min(i + row_block, row_end)
            n_floats = sum(r + 1 for r in range(i, stop))
            buf = f.read(4 * n_floats)
            if len(buf) != 4 * n_floats:
                raise RuntimeError(f"{path}: short read at row {i}")
            block = np.frombuffer(buf, dtype=np.float32)
            pos = 0
            for r in range(i, stop):
                row = block[pos:pos + r + 1]
                pos += r + 1
                offdiag = row[:r]                      # row[r] is the diagonal
                idx = np.nonzero(offdiag > threshold)[0]
                for j in idx:
                    hits.append((r, int(j), float(offdiag[j])))
            i = stop
    return row_start, row_end, hits

TEST_K = 16

t0 = time.monotonic()
rs, re_, hits = scan_shard(shard_path(TEST_K), TEST_K, THRESHOLD)
elapsed = time.monotonic() - t0

print(f"shard {TEST_K}/{N_SHARDS}: rows [{rs}, {re_}), {re_ - rs} rows")
print(f"{elapsed:.0f}s  ->  full scan ~{elapsed * N_SHARDS / 60:.0f} min")
print(f"{len(hits)} pairs above {THRESHOLD}")

if hits:
    vals = np.array([h[2] for h in hits])
    print(f"a_ij range: {vals.min():.3f} - {vals.max():.3f}")
    print("\ncounts by band:")
    for lo, hi, label in [(0.2, 0.35, "2nd deg"), (0.35, 0.7, "1st deg"),
                          (0.7, 1.5, "dup/MZ")]:
        print(f"  [{lo}, {hi}) {label:>8}: {((vals >= lo) & (vals < hi)).sum()}")
```

Expect clustering near 0.25, 0.5 and possibly 1.0. A smooth continuum with no
structure means the row-range recovery is off — stop before the full scan.

If this shard takes more than a few minutes, copy the shards to local disk
first (92 GiB) and repoint `SHARD_DIR`:

```
gcloud storage cp "gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9/03_grm_shards/eur_D2/shards/*.grm.bin.*" ~/scratch_grm/shards/
```

## Cell 5 — full scan

```python
all_hits = []
total_floats = 0
t0 = time.monotonic()

for k in range(1, N_SHARDS + 1):
    path = shard_path(k)
    if not os.path.isfile(path):
        print(f"MISSING shard {k}")
        continue
    total_floats += os.path.getsize(path) // 4
    _, _, hits = scan_shard(path, k, THRESHOLD)
    all_hits.extend(hits)
    print(f"[{k}/{N_SHARDS}] {len(hits):>7} hits  "
          f"({len(all_hits):>8} total, {time.monotonic() - t0:.0f}s)")

expected = N_IDS * (N_IDS + 1) // 2
print(f"\nfloats read : {total_floats:,}")
print(f"expected    : {expected:,}")
print("MATCH" if total_floats == expected else
      f"MISMATCH (diff {total_floats - expected:+,}) — shards missing or N_SHARDS wrong")
print(f"\n{len(all_hits)} pairs above {THRESHOLD}")
```

A float-count mismatch means the scan did not cover the whole matrix — do not
trust the pair list until it reconciles.

## Cell 6 — write outputs

```python
import pandas as pd

pairs = pd.DataFrame(
    [(ids[i][0], ids[i][1], ids[j][0], ids[j][1], a) for i, j, a in all_hits],
    columns=["FID1", "IID1", "FID2", "IID2", "grm_a"],
)

pairs_path = os.path.join(OUT_DIR, f"candidate_pairs_gt{THRESHOLD}.tsv")
pairs.to_csv(pairs_path, sep="\t", index=False)

involved = set(zip(pairs["FID1"], pairs["IID1"])) | set(zip(pairs["FID2"], pairs["IID2"]))
keep_path = os.path.join(OUT_DIR, f"candidate_individuals_gt{THRESHOLD}.keep")
with open(keep_path, "w") as f:
    for fid, iid in sorted(involved):
        f.write(f"{fid}\t{iid}\n")

print(f"{len(pairs)} pairs        -> {pairs_path}")
print(f"{len(involved)} individuals -> {keep_path}")
print(f"\nreduction: {len(involved)} / {N_IDS} = {100 * len(involved) / N_IDS:.1f}% of cohort")
print(f"KING pair-space: {(len(involved) / N_IDS) ** 2 * 100:.2f}% of original")
```

## Cell 7 — distribution check

```python
vals = pairs["grm_a"].values
print(f"{len(vals)} pairs, a_ij from {vals.min():.3f} to {vals.max():.3f}\n")

edges = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.55, 0.7, 0.9, 1.1, 2.0]
for lo, hi in zip(edges[:-1], edges[1:]):
    n = ((vals >= lo) & (vals < hi)).sum()
    bar = "#" * min(60, int(60 * n / max(1, len(vals)) * 4))
    print(f"[{lo:.2f}, {hi:.2f}) {n:>8}  {bar}")

deg1 = ((vals >= 0.35) & (vals < 0.7)).sum()
print(f"\n1st-degree band (0.35-0.7): {deg1} pairs -- these are what KING must split PO vs FS")
```

## Cell 8 — install plink2, build the KING panel

```bash
%%bash
set -e
BIN_DIR="$HOME/bin"
mkdir -p "$BIN_DIR"

if [ ! -x "$BIN_DIR/plink2" ]; then
  # URL is dated; if it 404s get the current one from
  # https://www.cog-genomics.org/plink/2.0/
  PLINK2_URL="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
  cd /tmp
  wget -q -O plink2.zip "$PLINK2_URL"
  unzip -o -q plink2.zip plink2 -d "$BIN_DIR"
  chmod +x "$BIN_DIR/plink2"
fi

export PATH="$BIN_DIR:$PATH"
plink2 --version
```

```python
import os
PANEL = os.path.expanduser(
    "~/workspace/Data from All of Us Controlled Tier /shared-env-pilot"
    "/phenotypic_covariance_v9/03_grm_shards/eur_D2/grm_input/eur_D2_GRM_QC"
)
os.environ["PANEL"] = PANEL
os.environ["SCREEN_DIR"] = OUT_DIR
os.environ["KEEP_PATH"] = keep_path
os.environ["KING_N_SNPS_TARGET"] = "50000"
print(PANEL)
```

```bash
%%bash
set -e
export PATH="$HOME/bin:$PATH"
SCREEN_DIR="$HOME/scratch_grm/relatedness_screen"
KEEP_PATH="${SCREEN_DIR}/candidate_individuals_gt0.2.keep"
PANEL="$HOME/workspace/Data from All of Us Controlled Tier /shared-env-pilot/phenotypic_covariance_v9/03_grm_shards/eur_D2/grm_input/eur_D2_GRM_QC"
KING_N_SNPS_TARGET=50000

N_KEEP=$(wc -l < "$KEEP_PATH")
N_CURRENT=$(wc -l < "${PANEL}.bim")
THIN_P=$(python3 -c "print(min(1.0, ${KING_N_SNPS_TARGET} / ${N_CURRENT}))")
echo "keep $N_KEEP individuals; thin $N_CURRENT -> ~${KING_N_SNPS_TARGET} (p=$THIN_P)"

plink2 --bfile "$PANEL" \
  --keep "$KEEP_PATH" \
  --thin "$THIN_P" --seed 1 \
  --make-bed --threads "$(nproc)" \
  --out "${SCREEN_DIR}/king_panel"

echo "variants: $(wc -l < "${SCREEN_DIR}/king_panel.bim")"
echo "samples : $(wc -l < "${SCREEN_DIR}/king_panel.fam")"
```

This reads the full 61 GB panel over gcsfuse (~40 min). The "Writing ... 0%"
progress counter tracks the tiny output, not the large input, so it looks
stalled while it is in fact working. Verify completion by exact byte count:
`3 + n_variants * ceil(n_samples / 4)`.

## Cell 9 — kinship + IBS0 on the subset

`plink2 --make-king-table` gives kinship *and* IBS0 and is already installed.
It has no `InfType`, so classify PO vs FS on `IBS0/NSNP` (expect bimodal:
PO ≈ 0, FS clearly above).

```bash
%%bash
set -e
export PATH="$HOME/bin:$PATH"
SCREEN_DIR="$HOME/scratch_grm/relatedness_screen"

plink2 --bfile "${SCREEN_DIR}/king_panel" \
  --make-king-table --king-table-filter 0.15 \
  --threads "$(nproc)" \
  --out "${SCREEN_DIR}/king_subset"

echo "related pairs: $(( $(wc -l < "${SCREEN_DIR}/king_subset.kin0") - 1 ))"
head -3 "${SCREEN_DIR}/king_subset.kin0"
```

Reads only the ~281 MB local panel, so this takes seconds.

## Cell 10 — classify PO vs FS

```python
kin = pd.read_csv(os.path.join(OUT_DIR, "king_subset.kin0"), sep=r"\s+")
kin.columns = [c.lstrip("#") for c in kin.columns]
# plink2 reports IBS0 and HETHET as PROPORTIONS, not counts (unlike the KING
# binary) -- do not divide by NSNP.
kin["ibs0_rate"] = kin["IBS0"]

deg1 = kin[(kin["KINSHIP"] >= 0.177) & (kin["KINSHIP"] < 0.354)].copy()
print(f"{len(deg1)} 1st-degree pairs")
print("\nibs0_rate quantiles (expect bimodal: PO near 0, FS clearly above):")
print(deg1["ibs0_rate"].quantile([0, .1, .25, .5, .75, .9, 1]).round(5))

edges = [0, 1e-4, 3e-4, 1e-3, 3e-3, 5e-3, 8e-3, 1.2e-2, 2e-2, 1]
for lo, hi in zip(edges[:-1], edges[1:]):
    n = ((deg1["ibs0_rate"] >= lo) & (deg1["ibs0_rate"] < hi)).sum()
    bar = "#" * min(60, int(60 * n / max(1, len(deg1))))
    print(f"[{lo:7.1e}, {hi:7.1e}) {n:>6}  {bar}")
```

Expected modes: PO near 1e-4 (IBS0 is genetically impossible for PO, so this
is just genotyping error), FS around 5e-3 to 1.7e-2. Pick the cut from the gap.

```python
IBS0_CUT = 1e-3   # <-- set from the histogram above

po = deg1[deg1["ibs0_rate"] < IBS0_CUT]
fs = deg1[deg1["ibs0_rate"] >= IBS0_CUT]
print(f"{len(po)} PO, {len(fs)} FS")

po_path = os.path.join(OUT_DIR, "po_pairs_exclude.tsv")
po[["IID1", "IID2", "KINSHIP", "IBS0", "NSNP", "ibs0_rate"]].to_csv(
    po_path, sep="\t", index=False)
print(f"-> {po_path}")
```

---

## Cross-check

Join `candidate_pairs_gt0.2.tsv` against `king_subset.kin0` on the IID pair
and compare `grm_a` against `2 * KINSHIP`. They should track closely;
disagreement points at the GRM or the ID alignment.

## Notes

- `.keep` and `.tsv` are gitignored — outputs are participant-level, never commit.
- `N_SHARDS` must match the run that produced the shards. Wrong value gives
  wrong row ranges, which the Cell 5 float-count check catches.
- Cell 8 reads the 61 GB panel over gcsfuse; `--keep` + `--thin` make the
  output tiny but the input read is still full-size.
