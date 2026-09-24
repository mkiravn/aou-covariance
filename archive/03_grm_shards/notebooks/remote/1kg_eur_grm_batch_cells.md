# GRM shard submission — 1kg_eur (Google Batch / dsub)

Sharded GRM (`plink --make-grm-bin --parallel k N_SHARDS`) for the 1kg_eur
round-2 sample set (n ≈ 223 K), run as Google Batch tasks via `dsub`. Each
task localizes the panel once and runs its assigned shards sequentially.

Panel: `1kg_eur_GRM_QC`, QC'd and filtered to round-2 keep list, staged in
`03_grm_shards/1kg_eur/grm_input/` (with precomputed `.frq`) — no variant
filtering in the shard jobs.

**Dimensions**: ~1.69M variants × 223,209 individuals.
Full GRM: `N(N+1)/2 × 4 bytes` ≈ **99.6 GB** across all shards.
BED file: ≈ **94 GB** (estimate — update `BED_SIZE_GB` from `gcloud storage ls`
before submitting).

**Interactive VM**: this notebook only submits and polls. Use the smallest
environment (2 vCPU / 8 GB); all work happens on Batch workers.

## The wrapper-image patch — read this first

`dsub` hardcodes `CLOUD_SDK_IMAGE` in `<dsub>/providers/google_utils.py` for
its five wrapper runnables. No dsub release works out of the box: 0.5.3/0.5.4
pin `499.0.0-slim`, which gcr.io has purged. Only recent tags (~581+) are still
published. `--image` does NOT fix this — it only sets the user-command runnable.

The patch lives in site-packages and ANY `pip install dsub` reverts it. Every
submission cell below therefore patches and submits **in the same shell**.

`581.0.0-slim` confirmed working 2026-09-10. Check current tags:
```
gcloud container images list-tags gcr.io/google.com/cloudsdktool/cloud-sdk --limit 20
```

---

## Cell 1 — prereqs

```bash
%%bash
set -e
pip install --quiet 'dsub==0.5.4'
dsub --version

DSUB_DIR=$(python -c 'import dsub,os;print(os.path.dirname(dsub.__file__))')
echo "--- wrapper copy tool (must be 'gcloud storage cp', NOT 'gsutil') ---"
grep -o 'gcloud storage cp\|gsutil .*cp' "$DSUB_DIR/providers/google_utils.py" | sort -u

echo "--- gcloud config ---"
gcloud config list --format='text(core.project,compute.region)' 2>&1 || true
```

## Cell 2 — config

```python
import math, os

PROJECT_ID = "wb-swift-sprout-7231"
REGION     = "us-central1"
SERVICE_ACCOUNT = "pet-27799165194323faf22e2@wb-swift-sprout-7231.iam.gserviceaccount.com"
NETWORK    = f"projects/{PROJECT_ID}/global/networks/network"
SUBNETWORK = f"projects/{PROJECT_ID}/regions/{REGION}/subnetworks/subnetwork"

WORKSPACE_BUCKET_GS = "gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
PROJECT_DIR = "phenotypic_covariance_v9"
SAMPLE_SET  = "1kg_eur"

BUCKET_DIR_GS    = f"{WORKSPACE_BUCKET_GS}/{PROJECT_DIR}/03_grm_shards/{SAMPLE_SET}"
GRM_INPUT_DIR_GS = f"{BUCKET_DIR_GS}/grm_input"
SHARD_OUT_DIR_GS = f"{BUCKET_DIR_GS}/shards"
PLINK_BIN_GS     = f"{BUCKET_DIR_GS}/bin/plink"

BED_NAME = "1kg_eur_GRM_QC"   # update if the staged panel has a different name

LIVE_TAG = "581.0.0-slim"

# BED file is ~94 GB (1.69M variants × 223K samples); confirm with:
#   gcloud storage ls -l "${GRM_INPUT_DIR_GS}/${BED_NAME}.bed"
# and update before submitting.
BED_SIZE_GB   = 94            # PLACEHOLDER — verify before submitting
MACHINE_VCPUS = 16
MEMORY_MB     = math.ceil((BED_SIZE_GB * 2 + 4) * 1024 / 256) * 256
PLINK_MEM_MB  = MEMORY_MB - 8192

_mem_per_vcpu = MEMORY_MB // MACHINE_VCPUS
MACHINE_TYPE  = (f"n1-custom-{MACHINE_VCPUS}-{MEMORY_MB}-ext" if _mem_per_vcpu > 8192
                 else f"n1-custom-{MACHINE_VCPUS}-{MEMORY_MB}")

# N scales computation as N²; 223K vs 155K (eur_D2) ≈ 2× harder per shard.
# N_SHARDS=20/N_TASKS=5 keeps 4 shards/task (same as eur_D2) but adds capacity.
# Watch the first shard's wall-clock and recalibrate if needed.
N_SHARDS     = 20
N_TASKS      = 5
DISK_SIZE_GB = 200   # 94 GB panel + ~5 GB per shard output; .grm.N.bin dropped

for k, v in dict(
    PROJECT_ID=PROJECT_ID, REGION=REGION, SERVICE_ACCOUNT=SERVICE_ACCOUNT,
    NETWORK=NETWORK, SUBNETWORK=SUBNETWORK,
    BUCKET_DIR_GS=BUCKET_DIR_GS, GRM_INPUT_DIR_GS=GRM_INPUT_DIR_GS,
    SHARD_OUT_DIR_GS=SHARD_OUT_DIR_GS, PLINK_BIN_GS=PLINK_BIN_GS,
    BED_NAME=BED_NAME, LIVE_TAG=LIVE_TAG,
    MACHINE_TYPE=MACHINE_TYPE, MACHINE_VCPUS=str(MACHINE_VCPUS),
    MEMORY_MB=str(MEMORY_MB), PLINK_MEM_MB=str(PLINK_MEM_MB),
    N_SHARDS=str(N_SHARDS), N_TASKS=str(N_TASKS), DISK_SIZE_GB=str(DISK_SIZE_GB),
).items():
    os.environ[k] = v

print(f"machine    : {MACHINE_TYPE}  ({_mem_per_vcpu} MB/vCPU)")
print(f"plink      : --memory {PLINK_MEM_MB} --threads {MACHINE_VCPUS}")
print(f"sharding   : {N_SHARDS} shards / {N_TASKS} tasks "
      f"({N_SHARDS / N_TASKS:.1f} shards/task)")
print(f"GRM size   : {223209 * 223210 // 2 * 4 / 1e9:.1f} GB expected")
print(f"output     : {SHARD_OUT_DIR_GS}")
```

## Cell 3 — verify inputs

`.bim` line count must equal `.frq` line count minus 1 (header). If the panel
is not staged yet, run Cell 3b first.

```bash
%%bash
set -e
echo "=== panel ==="
for ext in bed bim fam; do
  gcloud storage ls -l "${GRM_INPUT_DIR_GS}/${BED_NAME}.${ext}" 2>/dev/null \
    || echo "  MISSING: ${BED_NAME}.${ext}"
done
echo
echo "=== precomputed frequencies ==="
gcloud storage ls -l "${GRM_INPUT_DIR_GS}/${BED_NAME}_freq.frq" 2>/dev/null \
  || echo "  MISSING: ${BED_NAME}_freq.frq — run Cell 3b"
echo
echo "=== plink binary ==="
gcloud storage ls -l "$PLINK_BIN_GS" 2>/dev/null \
  || echo "  not staged yet — run Cell 4"
echo
echo "=== sample count in .fam ==="
gcloud storage cat "${GRM_INPUT_DIR_GS}/${BED_NAME}.fam" 2>/dev/null | wc -l \
  || echo "  (fam not readable)"
```

## Cell 3b — stage GRM panel from genome-wide QC output (one-time)

Run this only if the panel is not already in GCS. Pulls the genome-wide panel
from `01_ancestry_filtering/`, applies per-sample-set QC (MAF/HWE/missingness),
filters to the round-2 keep list, and uploads. Run Cell 4 first so that
`$HOME/bin/plink` exists for the frequency step (plink2 outputs .afreq not .frq,
but the shard jobs need plink-1.9-format .frq via `--read-freq`).

```bash
%%bash
set -e
PANEL_GS="gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9/01_ancestry_filtering/genome_wide_panel/genome_wide_panel_v9"
KEEP_GS="gs://cloned-shared-env-pilot-wb-swift-sprout-7231/phenotypic_covariance_v9/eur_r2/01_ancestry/round2/eur_r2_keep_ids.txt"

LOCAL_STAGE="$HOME/scratch_1kg_eur_grm_stage"
mkdir -p "$LOCAL_STAGE"

PLINK2="$HOME/bin/plink2"   # placed by env.ensure_plink2() in gate.ipynb / covariate_pcs.ipynb
PLINK1="$HOME/bin/plink"    # placed by Cell 4 — run Cell 4 first

if [ ! -x "$PLINK2" ]; then echo "ERROR: plink2 not found — run env.ensure_plink2() first"; exit 1; fi
if [ ! -x "$PLINK1" ]; then echo "ERROR: plink not found — run Cell 4 first"; exit 1; fi

# localize the genome-wide panel (bed format)
for ext in bed bim fam; do
  [ -f "$LOCAL_STAGE/gwpanel.$ext" ] || \
    gcloud storage cp "${PANEL_GS}.$ext" "$LOCAL_STAGE/gwpanel.$ext"
done
gcloud storage cp "$KEEP_GS" "$LOCAL_STAGE/keep.txt"

# QC in the round-2 cohort and filter to keep list
"$PLINK2" \
  --bfile "$LOCAL_STAGE/gwpanel" \
  --keep "$LOCAL_STAGE/keep.txt" --nonfounders \
  --maf 0.01 --hwe 1e-6 0.001 keep-fewhet --geno 0.05 \
  --max-alleles 2 \
  --threads $(nproc) \
  --make-bed --out "$LOCAL_STAGE/1kg_eur_GRM_QC"

# compute frequencies — plink 1.9 so the shard jobs can --read-freq it
"$PLINK1" \
  --bfile "$LOCAL_STAGE/1kg_eur_GRM_QC" \
  --freq --out "$LOCAL_STAGE/1kg_eur_GRM_QC_freq"

echo "samples: $(wc -l < "$LOCAL_STAGE/1kg_eur_GRM_QC.fam")"
echo "variants: $(wc -l < "$LOCAL_STAGE/1kg_eur_GRM_QC.bim")"

# upload
for f in "$LOCAL_STAGE/1kg_eur_GRM_QC.bed" \
          "$LOCAL_STAGE/1kg_eur_GRM_QC.bim" \
          "$LOCAL_STAGE/1kg_eur_GRM_QC.fam" \
          "$LOCAL_STAGE/1kg_eur_GRM_QC_freq.frq"; do
  gcloud storage cp "$f" "${GRM_INPUT_DIR_GS}/"
done
echo "staged to ${GRM_INPUT_DIR_GS}/"
```

## Cell 4 — stage the plink 1.9 binary (one-time)

Must be PLINK 1.9: `GRM-pairs`' row-range recovery is calibrated to 1.9's
`--parallel` split algorithm.

```bash
%%bash
set -e
BIN_DIR="$HOME/bin"; mkdir -p "$BIN_DIR"
if [ ! -x "$BIN_DIR/plink" ]; then
  # Try copying from an already-staged copy in the bucket (e.g. from eur_D2).
  # VPC-SC blocks external downloads (s3.amazonaws.com is not reachable).
  EXISTING_GS="$(echo "$PLINK_BIN_GS" | sed 's|/1kg_eur/|/eur_D2/|')"
  gcloud storage cp "$EXISTING_GS" "$BIN_DIR/plink" 2>/dev/null && \
    chmod +x "$BIN_DIR/plink" || {
    echo "Could not find a staged binary. Upload plink 1.9 manually:"
    echo "  gcloud storage cp gs://<path-to-plink1.9> $PLINK_BIN_GS"
    exit 1
  }
fi
"$BIN_DIR/plink" --version
gcloud storage cp "$BIN_DIR/plink" "$PLINK_BIN_GS"
gcloud storage ls -l "$PLINK_BIN_GS"
```

## Cell 5 — smoke test

Patches dsub, submits an `echo` job on a small machine. Confirms wrapper images
pull, logs delocalize, and a no-external-IP worker can reach
`storage.googleapis.com`. Do not proceed to the shard batch until SUCCESS.

```bash
%%bash
set -e
LIVE_TAG="581.0.0-slim"
DSUB_DIR=$(python -c 'import dsub,os;print(os.path.dirname(dsub.__file__))')
sed -i -E "s|cloud-sdk:[0-9]+\.[0-9]+\.[0-9]+-slim|cloud-sdk:${LIVE_TAG}|g" \
  "$DSUB_DIR/providers/google_utils.py"
find "$DSUB_DIR" -name '*.pyc' -delete
grep CLOUD_SDK_IMAGE "$DSUB_DIR/providers/google_utils.py"

dsub \
  --provider google-batch --project "$PROJECT_ID" --regions "$REGION" \
  --logging "${BUCKET_DIR_GS}/dsub_logs" \
  --service-account "$SERVICE_ACCOUNT" \
  --network "$NETWORK" --subnetwork "$SUBNETWORK" --use-private-address \
  --name "eur1kg-smoke" \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:${LIVE_TAG}" \
  --env BUCKET_DIR_GS="$BUCKET_DIR_GS" \
  --command '
    set -x
    hostname; date
    echo ok > /tmp/m.txt
    gcloud storage cp /tmp/m.txt "${BUCKET_DIR_GS}/dsub_logs/smoke_$(date +%s).txt" \
      && echo WROTE_OK || echo WRITE_FAILED
  ' \
  > /tmp/eur1kg_smoke_job_id.txt
cat /tmp/eur1kg_smoke_job_id.txt
```

Wait ~3 minutes, then check:

```bash
%%bash
JOB_ID=$(cat /tmp/eur1kg_smoke_job_id.txt)
dstat --provider google-batch --project "$PROJECT_ID" --location "$REGION" \
  --jobs "$JOB_ID" --users '*' --status '*'
echo "--- wrapper images ---"
gcloud batch jobs describe "${JOB_ID}-0-0" --location "$REGION" --project "$PROJECT_ID" \
  --format="value(taskGroups[0].taskSpec.runnables[].container.imageUri)"
echo "--- marker ---"
gcloud storage ls "${BUCKET_DIR_GS}/dsub_logs/smoke_"'*' 2>/dev/null || echo "(none)"
```

## Cell 6 — build the task list

`k` values assigned round-robin: `--parallel` balances by off-diagonal pair
count (not row count), so round-robin spreads cost variance across tasks.

```python
task_shards = {i: [] for i in range(N_TASKS)}
for k in range(1, N_SHARDS + 1):
    task_shards[(k - 1) % N_TASKS].append(k)

TASKS_PATH = "/tmp/eur1kg_grm_tasks.tsv"
with open(TASKS_PATH, "w") as f:
    f.write("--env K_LIST\n")
    for i in range(N_TASKS):
        f.write(",".join(str(k) for k in task_shards[i]) + "\n")

os.environ["TASKS_PATH"] = TASKS_PATH
print(open(TASKS_PATH).read())
for i, ks in task_shards.items():
    print(f"task {i}: {len(ks)} shards -> {ks}")
```

## Cell 7 — submit the GRM shard batch

Each task localizes the panel once and runs its shards sequentially.
`.grm.N.bin` dropped immediately after each shard (only `.grm.bin` is kept).

**Watch task 0's first shard wall-clock before assuming total runtime.**
If a shard takes > 2h, increase `N_SHARDS` and resubmit from Cell 6.

```bash
%%bash
set -e
DRY_RUN=""   # set to "--dry-run" to inspect without submitting
LIVE_TAG="581.0.0-slim"
DSUB_DIR=$(python -c 'import dsub,os;print(os.path.dirname(dsub.__file__))')
sed -i -E "s|cloud-sdk:[0-9]+\.[0-9]+\.[0-9]+-slim|cloud-sdk:${LIVE_TAG}|g" \
  "$DSUB_DIR/providers/google_utils.py"
find "$DSUB_DIR" -name '*.pyc' -delete
grep CLOUD_SDK_IMAGE "$DSUB_DIR/providers/google_utils.py"

echo "machine-type: $MACHINE_TYPE"

dsub \
  --provider google-batch --project "$PROJECT_ID" --regions "$REGION" \
  --logging "${BUCKET_DIR_GS}/dsub_logs" \
  --service-account "$SERVICE_ACCOUNT" \
  --network "$NETWORK" --subnetwork "$SUBNETWORK" --use-private-address \
  --machine-type "$MACHINE_TYPE" --disk-size "$DISK_SIZE_GB" \
  --name "eur1kg-grm-shards" \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:${LIVE_TAG}" \
  --input-recursive BED_DIR="$GRM_INPUT_DIR_GS" \
  --input PLINK_BIN="$PLINK_BIN_GS" \
  --env BED_NAME="$BED_NAME" \
  --env N_SHARDS="$N_SHARDS" \
  --output-recursive OUT_DIR="$SHARD_OUT_DIR_GS" \
  --command '
    set -e
    chmod +x "$PLINK_BIN"
    BED_PREFIX="${BED_DIR}/${BED_NAME}"
    FREQ_PATH="${BED_PREFIX}_freq.frq"

    IFS="," read -ra KS <<< "$K_LIST"
    for K in "${KS[@]}"; do
      echo "=== shard $K of $N_SHARDS : $(date -u +%H:%M:%S) ==="
      "$PLINK_BIN" \
        --bfile "$BED_PREFIX" \
        --make-grm-bin \
        --parallel "$K" "$N_SHARDS" \
        --read-freq "$FREQ_PATH" \
        --memory '"$PLINK_MEM_MB"' \
        --threads '"$MACHINE_VCPUS"' \
        --out "${OUT_DIR}/grm_shard_${K}_of_${N_SHARDS}"
      rm -f "${OUT_DIR}/grm_shard_${K}_of_${N_SHARDS}.grm.N.bin.${K}"
      echo "=== shard $K done : $(date -u +%H:%M:%S) ==="
    done
  ' \
  --tasks "$TASKS_PATH" \
  $DRY_RUN \
  > /tmp/eur1kg_shards_job_id.txt
cat /tmp/eur1kg_shards_job_id.txt
```

## Cell 8 — poll status

```bash
%%bash
JOB_ID=$(cat /tmp/eur1kg_shards_job_id.txt)
dstat --provider google-batch --project "$PROJECT_ID" --location "$REGION" \
  --jobs "$JOB_ID" --users '*' --status '*'
echo "--- shards written so far ---"
gcloud storage ls "${SHARD_OUT_DIR_GS}/" 2>/dev/null | grep -c 'grm.bin' || echo 0
```

## Cell 9 — per-task status and wall-clock

```python
import json, subprocess
from datetime import datetime

JOB_ID = open("/tmp/eur1kg_shards_job_id.txt").read().strip()
out = subprocess.run(
    ["dstat", "--provider", "google-batch", "--project", PROJECT_ID,
     "--location", REGION, "--jobs", JOB_ID, "--users", "*",
     "--status", "*", "--full", "--format", "json"],
    capture_output=True, text=True, check=True).stdout
tasks = json.loads(out)

def ts(x):
    return datetime.fromisoformat(x.replace("Z", "+00:00")) if x else None

print(f"{'task':>4}  {'status':<10} {'wall_h':>7}  k_list")
walls = []
for t in sorted(tasks, key=lambda t: int(t.get("task-id", 0) or 0)):
    c, e = ts(t.get("create-time")), ts(t.get("end-time"))
    w = (e - c).total_seconds() / 3600 if c and e else None
    if w:
        walls.append(w)
    print(f"{t.get('task-id'):>4}  {t.get('status'):<10} "
          f"{(f'{w:.2f}' if w else '-'):>7}  {t.get('envs', {}).get('K_LIST', '?')}")

n_ok = sum(1 for t in tasks if t.get("status") == "SUCCESS")
print(f"\n{n_ok}/{len(tasks)} tasks SUCCESS")
if walls:
    print(f"wall-clock {min(walls):.2f}h – {max(walls):.2f}h")
```

## Cell 10 — verify outputs

Total `.grm.bin` bytes must equal `N(N+1)/2 × 4`. Shortfall = shard failed;
excess = rows double-counted.

```python
import re, subprocess

listing = subprocess.run(
    ["gcloud", "storage", "ls", "-l", f"{SHARD_OUT_DIR_GS}/"],
    capture_output=True, text=True, check=True).stdout

sizes = {}
for line in listing.splitlines():
    m = re.match(r"\s*(\d+)\s+\S+\s+(\S+)", line)
    if m and ".grm.bin." in m.group(2):
        sizes[m.group(2)] = int(m.group(1))

print(f"shard .grm.bin files: {len(sizes)} / {N_SHARDS}")
missing = [k for k in range(1, N_SHARDS + 1)
           if not any(f"grm_shard_{k}_of_{N_SHARDS}.grm.bin.{k}" in p for p in sizes)]
if missing:
    print(f"MISSING shards: {missing}")

total = sum(sizes.values())
print(f"total .grm.bin: {total / 1e9:.2f} GB")

n_ind = int(subprocess.run(
    ["gcloud", "storage", "cat", f"{GRM_INPUT_DIR_GS}/{BED_NAME}.fam"],
    capture_output=True, text=True, check=True).stdout.count("\n"))
expected = n_ind * (n_ind + 1) // 2 * 4
print(f"expected for N={n_ind}: {expected / 1e9:.2f} GB")
print("MATCH" if total == expected else f"MISMATCH (diff {(total - expected) / 1e9:+.2f} GB)")
```

```bash
%%bash
echo "=== .grm.id (first 3 rows) ==="
gcloud storage cat "${SHARD_OUT_DIR_GS}/grm_shard_1_of_${N_SHARDS}.grm.id" 2>/dev/null | head -3 \
  || gcloud storage ls "${SHARD_OUT_DIR_GS}/" | grep 'grm.id'
```

## Cell 11 — task logs

```bash
%%bash
JOB_ID=$(cat /tmp/eur1kg_shards_job_id.txt)
for f in $(gcloud storage ls "${BUCKET_DIR_GS}/dsub_logs/${JOB_ID}"'*' 2>/dev/null); do
  echo "=================== $f ==================="
  gcloud storage cat "$f"
done
```

---

## Next steps

1. Check `.grm.id` FID column matches `write_grm_pheno()` output (FID = IID = person_id).
2. Run `04_process_shards` — `grm_shard_tool accumulate` merges shards into the final GRM.
3. Relatedness handling: GRM-screen pairs near a_ij ≈ 0.5, then KING to separate PO vs FS; separate PO accumulator for PO mean phenotypic covariance output.

**If submission fails with exit 1 and no logs**: the patch reverted. Check
wrapper images with the `gcloud batch jobs describe` line in Cell 5 and re-run
the failing submission cell (each one re-patches itself).
