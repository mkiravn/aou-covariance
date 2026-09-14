# GRM shard submission — eur_D2 (Google Batch / dsub)

Sharded GRM (`plink --make-grm-bin --parallel k N_SHARDS`) for the eur_D2
sample set, run as Google Batch tasks via `dsub`. Each task gets its own VM,
localizes the panel once, and runs its assigned shards in sequence.

Panel: `eur_D2_GRM_QC`, already QC'd and staged in
`03_grm_shards/eur_D2/grm_input/` (with precomputed `.frq`) — no variant
filtering here.

**Dimensions**: ~1.69M variants x ~155K individuals. Full GRM is
`N(N+1)/2 x 4 bytes` ~= 48 GB across all shards.

**Interactive VM**: this notebook only submits and polls. Use the smallest
environment (2 vCPU / 8 GB); all work happens on Batch workers.

## The wrapper-image patch — read this first

`dsub` hardcodes `CLOUD_SDK_IMAGE` in `<dsub>/providers/google_utils.py` for
its five wrapper runnables (localize / log-stream / delocalize). **No dsub
release works out of the box here**: 0.5.2 pins `294.0.0-slim`, 0.5.3 and
0.5.4 pin `499.0.0-slim`, and gcr.io has purged both. Only recent tags (~581+)
are still published.

Two failure modes, identical symptom:

1. **Dead image tag** — nothing can pull.
2. **dsub 0.5.2 on a live image** — 0.5.2's wrappers call `gsutil` and a bare
   `python`, and current cloud-sdk images ship neither. Under `set -o errexit`
   the logging runnable dies on `gsutil: command not found`.

Both produce: task reaches RUNNING, dies after ~40-60s with exit code 1, and
writes **zero** logs — the logging runnable never gets far enough to upload
anything, so `dstat --full` shows only "exit code 1" with no stderr anywhere.

So the requirements collide: an old dsub needs an old image, but old images are
gone. **The only working combination is dsub >= 0.5.3 patched to a live tag.**

`--image` does **not** fix either mode; it only sets the sixth, user-command
runnable.

The patch lives in site-packages, so any `pip install dsub` silently reverts
it. Every submission cell below therefore patches and submits **in the same
shell**. Do not factor the patch out into its own cell.

`581.0.0-slim` is confirmed working with dsub 0.5.3 and 0.5.4 (2026-09-10).
Check what is still published with:

```
gcloud container images list-tags gcr.io/google.com/cloudsdktool/cloud-sdk --limit 20
```

## Cell 1 — prereqs

**Must be dsub >= 0.5.3.** 0.5.2's wrapper scripts call `gsutil` and a bare
`python`, neither of which exists in a current cloud-sdk image — and current
images are the only ones still published (see the next section). 0.5.3+ uses
`gcloud storage` and `python3` instead. Pinning 0.5.2 produces exit 1 with
zero logs even when the image tag is live.

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
import math
import os

PROJECT_ID = "wb-swift-sprout-7231"
REGION     = "us-central1"

# This VM's attached "pet" SA. The calling identity can actAs this one; it
# cannot actAs the project's default Compute Engine SA (PERMISSION_DENIED).
SERVICE_ACCOUNT = "pet-27799165194323faf22e2@wb-swift-sprout-7231.iam.gserviceaccount.com"

# VPC-SC is enabled: Batch tasks MUST run with no external IP, on the
# Workbench-standard network/subnetwork names. Dropping --use-private-address
# is rejected outright at submit time.
NETWORK    = f"projects/{PROJECT_ID}/global/networks/network"
SUBNETWORK = f"projects/{PROJECT_ID}/regions/{REGION}/subnetworks/subnetwork"

# Real bucket name, from `ps -eo pid,args | grep gcsfuse` -- not gs://<project-id>.
WORKSPACE_BUCKET_GS = "gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
PROJECT_DIR = "phenotypic_covariance_v9"
SAMPLE_SET  = "eur_D2"

BUCKET_DIR_GS    = f"{WORKSPACE_BUCKET_GS}/{PROJECT_DIR}/03_grm_shards/{SAMPLE_SET}"
GRM_INPUT_DIR_GS = f"{BUCKET_DIR_GS}/grm_input"
SHARD_OUT_DIR_GS = f"{BUCKET_DIR_GS}/shards"
PLINK_BIN_GS     = f"{BUCKET_DIR_GS}/bin/plink"

BED_NAME = "eur_D2_GRM_QC"

# Cloud SDK image tag used for BOTH the dsub wrappers (via the sed patch in
# each submission cell) and the user-command runnable (via --image).
LIVE_TAG = "581.0.0-slim"

# Machine sizing, same formula as the eur_base run: 2x bed size + 4 GB.
# plink's real RSS per --make-grm-bin shard with --read-freq is <1 GB -- the
# --memory flag is a virtual ceiling, not a reservation -- so this is generous
# on purpose. N1 custom machine types require memory to be an exact multiple
# of 256 MB, else Batch rejects the VM spec before any task or log runs.
BED_SIZE_GB = 65.45          # from `gcloud storage ls -l` on the staged .bed
MACHINE_VCPUS = 16
MEMORY_MB = math.ceil((BED_SIZE_GB * 2 + 4) * 1024 / 256) * 256
PLINK_MEM_MB = MEMORY_MB - 8192   # headroom for OS + input localization

# N1 custom needs an "-ext" suffix once memory exceeds 8192 MB/vCPU.
_mem_per_vcpu = MEMORY_MB // MACHINE_VCPUS
MACHINE_TYPE = (f"n1-custom-{MACHINE_VCPUS}-{MEMORY_MB}-ext" if _mem_per_vcpu > 8192
                else f"n1-custom-{MACHINE_VCPUS}-{MEMORY_MB}")

# Few, fat shards: per-shard cost is dominated by the fixed panel load, so more
# shards just pays that more times. Grouping into tasks amortizes the ~20 min
# per-task panel download across several shards of real work.
N_SHARDS = 16
N_TASKS  = 4

DISK_SIZE_GB = 150   # 61 GB panel + ~3 GB per shard output, .grm.N.bin dropped

for k, v in dict(
    PROJECT_ID=PROJECT_ID, REGION=REGION, SERVICE_ACCOUNT=SERVICE_ACCOUNT,
    NETWORK=NETWORK, SUBNETWORK=SUBNETWORK, BUCKET_DIR_GS=BUCKET_DIR_GS,
    GRM_INPUT_DIR_GS=GRM_INPUT_DIR_GS, SHARD_OUT_DIR_GS=SHARD_OUT_DIR_GS,
    PLINK_BIN_GS=PLINK_BIN_GS, BED_NAME=BED_NAME, LIVE_TAG=LIVE_TAG,
    MACHINE_TYPE=MACHINE_TYPE, MACHINE_VCPUS=str(MACHINE_VCPUS),
    MEMORY_MB=str(MEMORY_MB), PLINK_MEM_MB=str(PLINK_MEM_MB),
    N_SHARDS=str(N_SHARDS), N_TASKS=str(N_TASKS), DISK_SIZE_GB=str(DISK_SIZE_GB),
).items():
    os.environ[k] = v

print(f"machine    : {MACHINE_TYPE}  ({_mem_per_vcpu} MB/vCPU)")
print(f"plink      : --memory {PLINK_MEM_MB} --threads {MACHINE_VCPUS}")
print(f"sharding   : {N_SHARDS} shards / {N_TASKS} tasks "
      f"({N_SHARDS / N_TASKS:.1f} shards per task)")
print(f"output     : {SHARD_OUT_DIR_GS}")
```

## Cell 3 — verify inputs

`.bim` line count must equal `.frq` line count minus 1 (header), and the
`.frq` must have been computed on this same panel.

```bash
%%bash
set -e
echo "=== panel ==="
for ext in bed bim fam; do
  gcloud storage ls -l "${GRM_INPUT_DIR_GS}/${BED_NAME}.${ext}"
done
echo
echo "=== precomputed frequencies ==="
gcloud storage ls -l "${GRM_INPUT_DIR_GS}/${BED_NAME}_freq.frq"
echo
echo "=== plink binary ==="
gcloud storage ls -l "$PLINK_BIN_GS" 2>/dev/null || echo "  not staged yet -- run Cell 4"
```

## Cell 4 — stage the plink binary (one-time)

The Batch worker image has no `wget` or `curl`, so a task can't download plink
itself. Install it locally once and upload; `dsub` localizes it via `--input`.
Must be PLINK 1.9: `GRM-pairs`' row-range recovery is calibrated to 1.9's
`--parallel` split algorithm, and 1.9 is what produced every prior shard set.

```bash
%%bash
set -e
BIN_DIR="$HOME/bin"; mkdir -p "$BIN_DIR"
if [ ! -x "$BIN_DIR/plink" ]; then
  PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip"
  cd /tmp && wget -q -O plink.zip "$PLINK_URL"
  unzip -o -q plink.zip plink -d "$BIN_DIR"
  chmod +x "$BIN_DIR/plink"
fi
"$BIN_DIR/plink" --version
gcloud storage cp "$BIN_DIR/plink" "$PLINK_BIN_GS"
gcloud storage ls -l "$PLINK_BIN_GS"
```

## Cell 5 — smoke test

Cheapest possible real submission: patches dsub, submits an `echo` job on a
small machine, and has the task write a marker straight to GCS. Confirms in
one shot that the wrapper images pull, that logs delocalize, and that a
no-external-IP worker can reach `storage.googleapis.com`.

Do not proceed to the shard batch until this returns SUCCESS.

```bash
%%bash
set -e
# --- patch dsub's wrapper image, SAME SHELL as the submit ---
# dsub pins a cloud-sdk tag that gcr.io has purged (0.5.2 -> 294.0.0-slim,
# 0.5.3/0.5.4 -> 499.0.0-slim). --image only sets the user-command runnable,
# not the five wrappers, so an unpatched submission dies after ~40-60s with
# exit 1 and zero logs. The patch is in site-packages and ANY pip install of
# dsub reverts it, so it must never live in its own cell.
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
  --name "eurd2-smoke" \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:${LIVE_TAG}" \
  --env BUCKET_DIR_GS="$BUCKET_DIR_GS" \
  --command '
    set -x
    hostname; date
    echo ok > /tmp/m.txt
    gcloud storage cp /tmp/m.txt "${BUCKET_DIR_GS}/dsub_logs/smoke_$(date +%s).txt" \
      && echo WROTE_OK || echo WRITE_FAILED
  ' \
  > /tmp/eurd2_smoke_job_id.txt
cat /tmp/eurd2_smoke_job_id.txt
```

Wait ~3 minutes, then check. All six runnables must show the live tag.

```bash
%%bash
JOB_ID=$(cat /tmp/eurd2_smoke_job_id.txt)
dstat --provider google-batch --project "$PROJECT_ID" --location "$REGION" \
  --jobs "$JOB_ID" --users '*' --status '*'
echo "--- wrapper images ---"
gcloud batch jobs describe "${JOB_ID}-0-0" --location "$REGION" --project "$PROJECT_ID" \
  --format="value(taskGroups[0].taskSpec.runnables[].container.imageUri)"
echo "--- marker ---"
gcloud storage ls "${BUCKET_DIR_GS}/dsub_logs/smoke_"'*' 2>/dev/null || echo "(none)"
echo "--- logs ---"
gcloud storage ls "${BUCKET_DIR_GS}/dsub_logs/${JOB_ID}"'*' 2>/dev/null || echo "(none)"
```

## Cell 6 — build the task list

`k` values are assigned round-robin across tasks rather than in contiguous
blocks: `--parallel` balances shards by off-diagonal *pair* count, not row
count, so different `k` can cost somewhat differently. Round-robin spreads
that across all tasks instead of concentrating it in one.

```python
task_shards = {i: [] for i in range(N_TASKS)}
for k in range(1, N_SHARDS + 1):
    task_shards[(k - 1) % N_TASKS].append(k)

TASKS_PATH = "/tmp/eurd2_grm_tasks.tsv"
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

Each task localizes the panel once (~20 min for 61 GB) and then runs its
shards sequentially at `--threads 16`. `.grm.N.bin` is dropped immediately
after each shard: `grm_shard_tool` only reads `.grm.bin`, and keeping it would
double both output size and upload time.

Set `DRY_RUN="--dry-run"` to inspect the request without submitting.

Timing is an estimate extrapolated from the eur_base run — **watch task 0's
first shard** and recalculate before assuming the total.

```bash
%%bash
set -e
DRY_RUN=""
# --- patch dsub's wrapper image, SAME SHELL as the submit ---
# dsub pins a cloud-sdk tag that gcr.io has purged (0.5.2 -> 294.0.0-slim,
# 0.5.3/0.5.4 -> 499.0.0-slim). --image only sets the user-command runnable,
# not the five wrappers, so an unpatched submission dies after ~40-60s with
# exit 1 and zero logs. The patch is in site-packages and ANY pip install of
# dsub reverts it, so it must never live in its own cell.
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
  --name "eurd2-grm-shards" \
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
  > /tmp/eurd2_shards_job_id.txt
cat /tmp/eurd2_shards_job_id.txt
```

## Cell 8 — poll status

Re-run for a snapshot. Read the job id with `cat`, never `tail -1` — dsub's
stdout for a `--tasks` submission doesn't reliably end with a bare job id, and
a wrong id makes `dstat` silently return `[]`.

```bash
%%bash
JOB_ID=$(cat /tmp/eurd2_shards_job_id.txt)
dstat --provider google-batch --project "$PROJECT_ID" --location "$REGION" \
  --jobs "$JOB_ID" --users '*' --status '*'
echo "--- shards written so far ---"
gcloud storage ls "${SHARD_OUT_DIR_GS}/" 2>/dev/null | grep -c 'grm.bin' \
  || echo 0
```

## Cell 9 — per-task status and wall-clock

```python
import json
import subprocess
from datetime import datetime

JOB_ID = open("/tmp/eurd2_shards_job_id.txt").read().strip()
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
    print(f"wall-clock {min(walls):.2f}h - {max(walls):.2f}h")
```

## Cell 10 — verify outputs

Total `.grm.bin` bytes must equal `N(N+1)/2 x 4` regardless of `N_SHARDS`. A
shortfall means a shard silently failed; an excess means rows were
double-counted.

Also check the `.grm.id` FID column: `grm_shard_tool`'s pheno lookup keys on
the full `(FID, IID)` pair, while `write_grm_pheno()` sets
`FID = IID = person_id` — plink commonly defaults `FID` to `0`.

```python
import re
import subprocess

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

Dumps stdout/stderr for every task. The per-shard timestamps printed by the
`--command` loop give real per-shard wall-clock.

```bash
%%bash
JOB_ID=$(cat /tmp/eurd2_shards_job_id.txt)
for f in $(gcloud storage ls "${BUCKET_DIR_GS}/dsub_logs/${JOB_ID}"'*' 2>/dev/null); do
  echo "=================== $f ==================="
  gcloud storage cat "$f"
done
```

## Next steps

1. Confirm the `.grm.id` FID column above matches what `write_grm_pheno()`
   produces, before the phenotype cross-product trusts it.
2. Run `04_process_shards/notebooks/remote/01_grm_shard_processing.ipynb` —
   `grm_shard_tool accumulate` merges the shards into the final GRM.

To run another sample set, change `SAMPLE_SET` / `BED_NAME` / `BED_SIZE_GB` in
Cell 2 and re-run from Cell 3.

**If a submission fails with exit 1 and no logs**, the patch reverted — check
the wrapper images with the `gcloud batch jobs describe` line in Cell 5 and
re-run the failing submission cell (each one re-patches on its own).
