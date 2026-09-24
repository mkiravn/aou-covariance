# GRM shard construction — eur_D2 (local VM, production)

Complete cell set for running sharded GRM on a single `n1-highmem-64` VM via concurrent plink processes. No dsub, no Batch, no image registry headaches.

**Timing:** ~8 hours wall-clock. **Cost:** ~$24. **Machine:** n1-highmem-64, 250 GB pd-ssd boot disk, Debian 12.

**Run order:** cells 1–10 in sequence, on the VM (SSH or a notebook running on that VM).

---

## VM Setup (run locally, not on the VM)

```bash
#!/usr/bin/env bash
gcloud compute instances create eur-d2-grm-build \
  --project wb-swift-sprout-7231 \
  --zone us-central1-a \
  --machine-type n1-highmem-64 \
  --boot-disk-size 250 \
  --boot-disk-type pd-ssd \
  --image-family debian-12 \
  --image-project debian-cloud \
  --scopes https://www.googleapis.com/auth/cloud-platform \
  --metadata enable-oslogin=TRUE
```

Panel dims (eur_D2): ~1.69M variants × ~155K individuals. Peak local disk
~115 GB (61 GB bed + 48 GB `.grm.bin` shards; `.grm.N.bin` deleted per-shard).

Wait for it to finish, then SSH in:

```bash
#!/usr/bin/env bash
gcloud compute ssh eur-d2-grm-build --zone us-central1-a --project wb-swift-sprout-7231
```

---

## Cell 1 — Install plink

```bash
#!/usr/bin/env bash
set -e
BIN_DIR="$HOME/bin"
mkdir -p "$BIN_DIR"
if [ ! -x "$BIN_DIR/plink" ]; then
  PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip"
  cd /tmp
  wget -q -O plink.zip "$PLINK_URL"
  unzip -o -q plink.zip plink -d "$BIN_DIR"
  chmod +x "$BIN_DIR/plink"
fi
export PATH="$BIN_DIR:$PATH"
plink --version
nproc
```

## Cell 2 — Config paths (define once, reuse in all following cells)

Writes `$HOME/.grm_env.sh`; every later cell `source`s it.

```bash
#!/usr/bin/env bash
set -e
cat > ~/.grm_env.sh << 'EOF'
#!/usr/bin/env bash
export WORKSPACE_BUCKET="$HOME/workspace/Data from All of Us Controlled Tier /shared-env-pilot"
export CDR_VERSION="v9"
export PROJECT_DIR="phenotypic_covariance_v9"
export SAMPLE_SET="eur_D2"
# eur_D2's QC'd panel was staged directly into 03_grm_shards/eur_D2/grm_input/
# (not the 01_ancestry_filtering/genome_wide_panel_* location that eur_base uses)
export PANEL_DIR="${WORKSPACE_BUCKET}/${PROJECT_DIR}/03_grm_shards/${SAMPLE_SET}/grm_input"
export BUCKET_DIR="${WORKSPACE_BUCKET}/${PROJECT_DIR}/03_grm_shards/${SAMPLE_SET}"
export LOCAL_WORK_DIR="$HOME/scratch_grm"
export SHARD_OUT_DIR="${LOCAL_WORK_DIR}/shards"
export BED_NAME="eur_D2_GRM_QC"
export BED_PREFIX="${LOCAL_WORK_DIR}/${BED_NAME}"
export PATH="$HOME/bin:$PATH"
EOF

source ~/.grm_env.sh
mkdir -p "$SHARD_OUT_DIR"
echo "✓ Config loaded. PANEL_DIR=$PANEL_DIR"
ls -lh "$PANEL_DIR"/
```

## Cell 3 — Copy panel + frequencies from bucket to local scratch

The `.frq` is already computed and staged in `grm_input/` alongside the panel — copy it, don't recompute.

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
set -e

mkdir -p "$LOCAL_WORK_DIR"

for ext in bed bim fam; do
  src="${PANEL_DIR}/${BED_NAME}.${ext}"
  dest="${BED_PREFIX}.${ext}"
  if [ ! -f "$dest" ]; then
    echo "Copying ${ext} (~$(du -h "$src" | cut -f1))..."
    cp "$src" "$dest"
  else
    echo "✓ Already present: ${ext}"
  fi
done

# precomputed frequencies (staged with the panel)
if [ ! -f "${BED_PREFIX}_freq.frq" ]; then
  echo "Copying freq..."
  cp "${PANEL_DIR}/${BED_NAME}_freq.frq" "${BED_PREFIX}_freq.frq"
else
  echo "✓ Already present: freq"
fi

echo
ls -lh "${BED_PREFIX}"*
```

## Cell 4 — Confirm frequencies present

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
FREQ_PATH="${BED_PREFIX}_freq.frq"
if [ -f "$FREQ_PATH" ]; then
  echo "✓ $(wc -l < "$FREQ_PATH") lines in $FREQ_PATH"
  # line count minus 1 (header) should equal the .bim variant count
  echo "  .bim variants: $(wc -l < "${BED_PREFIX}.bim")"
else
  echo "MISSING — recompute with: plink --bfile $BED_PREFIX --freq --out ${BED_PREFIX}_freq"
fi
```

## Cell 5 — Calculate memory and save run params

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
set -e

# --memory ceiling: 2x bed size + 4 GB headroom. plink's actual RSS per
# --make-grm-bin shard with --read-freq is <1 GB, so this is just a cap.
BED_SIZE_BYTES=$(stat -c%s "${BED_PREFIX}.bed" 2>/dev/null || stat -f%z "${BED_PREFIX}.bed")
BED_SIZE_GB=$(echo "scale=1; $BED_SIZE_BYTES / 1000000000" | bc)
MEMORY_MB=$(echo "($BED_SIZE_GB * 2 + 4) * 1024" | bc | cut -d. -f1)

# Run params (proven on eur_base, n1-highmem-64)
N_SHARDS=300
N_CONCURRENT=20
N_THREADS=4

cat > ~/.grm_run_params.sh << EOF
#!/usr/bin/env bash
export FREQ_PATH="${BED_PREFIX}_freq.frq"
export MEMORY_MB=${MEMORY_MB}
export N_SHARDS=${N_SHARDS}
export N_CONCURRENT=${N_CONCURRENT}
export N_THREADS=${N_THREADS}
EOF

source ~/.grm_run_params.sh
echo "✓ Run params:"
echo "  BED_SIZE_GB=${BED_SIZE_GB}"
echo "  MEMORY_MB=${MEMORY_MB}"
echo "  N_SHARDS=${N_SHARDS}, N_CONCURRENT=${N_CONCURRENT}, N_THREADS=${N_THREADS}"
echo "  Expected wall-clock: ~8 hours (15 waves × ~1450s per shard)"
```

## Cell 6 — Timing check (k=1 only, real settings)

Runs one shard at the actual thread count. **Watch the number** — expect ~1400–1500 s. If it's 2–3× that, stop and diagnose (disk contention, wrong thread count) before launching all 300.

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh
set -e

echo "Timing shard k=1 with --threads=$N_THREADS..."
time plink --bfile "$BED_PREFIX" --make-grm-bin \
  --parallel 1 "$N_SHARDS" \
  --read-freq "$FREQ_PATH" --memory "$MEMORY_MB" \
  --threads "$N_THREADS" \
  --out "${SHARD_OUT_DIR}/grm_shard_1_of_${N_SHARDS}_test"

echo
echo "Expected total: $(python3 -c "import math; w=math.ceil($N_SHARDS/$N_CONCURRENT); print(f'{w} waves x 1450s = {w*1450/3600:.1f} h')")"
```

## Cell 7 — Write the shard-runner script

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh
set -e

cat > ~/run_all_shards.py << 'PYSCRIPT'
#!/usr/bin/env python3
import os
import subprocess
import time
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

BED_PREFIX    = os.path.expanduser("~/scratch_grm/eur_D2_GRM_QC")
SHARD_OUT_DIR = os.path.expanduser("~/scratch_grm/shards")
FREQ_PATH     = os.environ["FREQ_PATH"]
MEMORY_MB     = int(os.environ["MEMORY_MB"])
N_SHARDS      = int(os.environ["N_SHARDS"])
N_CONCURRENT  = int(os.environ["N_CONCURRENT"])
N_THREADS     = int(os.environ["N_THREADS"])

def run_shard(k):
    out_prefix = os.path.join(SHARD_OUT_DIR, f"grm_shard_{k}_of_{N_SHARDS}")
    local_bin = f"{out_prefix}.grm.bin.{k}"

    # skip if already done (re-run resilience)
    if os.path.isfile(local_bin) and os.path.getsize(local_bin) > 1000:
        return {"k": k, "elapsed_sec": 0, "status": "skipped", "rc": 0, "stderr_tail": ""}

    cmd = ["plink", "--bfile", BED_PREFIX, "--make-grm-bin",
           "--parallel", str(k), str(N_SHARDS),
           "--read-freq", FREQ_PATH, "--memory", str(MEMORY_MB),
           "--threads", str(N_THREADS),
           "--out", out_prefix]

    start = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.monotonic() - start

    # downstream grm_shard_tool only reads .grm.bin -- drop .grm.N.bin
    # immediately to keep peak local disk to ~115 GB, not ~160 GB
    if proc.returncode == 0:
        n_bin = f"{out_prefix}.grm.N.bin.{k}"
        if os.path.isfile(n_bin):
            os.remove(n_bin)

    return {
        "k": k,
        "elapsed_sec": elapsed,
        "status": "ok" if proc.returncode == 0 else "FAILED",
        "rc": proc.returncode,
        "stderr_tail": proc.stderr[-200:] if proc.returncode != 0 else "",
    }

results = []
start_all = time.monotonic()
print(f"Starting {N_SHARDS} shards, {N_CONCURRENT} concurrent, {N_THREADS} threads each", flush=True)
print(f"Expected: {(N_SHARDS + N_CONCURRENT - 1) // N_CONCURRENT} waves x ~1450s "
      f"= ~{(N_SHARDS + N_CONCURRENT - 1) // N_CONCURRENT * 1450 // 3600}h wall-clock\n", flush=True)

with ThreadPoolExecutor(max_workers=N_CONCURRENT) as pool:
    futures = {pool.submit(run_shard, k): k for k in range(1, N_SHARDS + 1)}
    for i, future in enumerate(as_completed(futures), 1):
        r = future.result()
        results.append(r)
        sym = {"skipped": "=", "ok": "+", "FAILED": "x"}[r["status"]]
        print(f"[{i:3}/{N_SHARDS}] k={r['k']:3}: {sym} {r['elapsed_sec']:7.1f}s", flush=True)

wall_h = (time.monotonic() - start_all) / 3600
n_ok   = sum(1 for r in results if r["status"] == "ok")
n_skip = sum(1 for r in results if r["status"] == "skipped")
n_fail = sum(1 for r in results if r["status"] == "FAILED")

print(f"\n{'='*60}", flush=True)
print(f"Completed in {wall_h:.1f} h wall-clock", flush=True)
print(f"  {n_ok} new + {n_skip} already present = {n_ok + n_skip}/{N_SHARDS} done", flush=True)
if n_fail:
    print(f"FAILED: {n_fail} shard(s)", flush=True)
    for r in results:
        if r["status"] == "FAILED":
            print(f"  k={r['k']}: rc={r['rc']}  {r['stderr_tail'][:120]}", flush=True)
    sys.exit(1)
print("All shards succeeded.", flush=True)
PYSCRIPT

chmod +x ~/run_all_shards.py
echo "✓ Wrote ~/run_all_shards.py"
```

## Cell 8 — Launch the run in the background

`nohup` + `&` so a browser/SSH disconnect doesn't kill the 8-hour run. Env
vars are exported into the script's process before it detaches.

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh

cd ~
nohup env \
  FREQ_PATH="$FREQ_PATH" \
  MEMORY_MB="$MEMORY_MB" \
  N_SHARDS="$N_SHARDS" \
  N_CONCURRENT="$N_CONCURRENT" \
  N_THREADS="$N_THREADS" \
  PATH="$PATH" \
  python3 ~/run_all_shards.py > ~/grm_run.log 2>&1 &

echo $! > ~/grm_run.pid
echo "✓ Launched PID $(cat ~/grm_run.pid)"
echo "  monitor:  tail -f ~/grm_run.log"
echo "  progress: grep -c '^\[' ~/grm_run.log   # shards completed"
echo "  stop:     kill \$(cat ~/grm_run.pid)"
```

## Cell 9 — Monitor progress

Re-run this cell whenever you want a status snapshot.

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh

PID=$(cat ~/grm_run.pid 2>/dev/null || echo "")
if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then
  echo "RUNNING (PID $PID)"
else
  echo "NOT RUNNING (finished or died — check tail below)"
fi

done_count=$(find "$SHARD_OUT_DIR" -name "grm_shard_*_of_${N_SHARDS}.grm.bin.*" 2>/dev/null | grep -v _test | wc -l)
echo "Shards done: ${done_count}/${N_SHARDS}"

echo "Disk:"
df -h "$HOME" | tail -1

echo "--- last 15 log lines ---"
tail -15 ~/grm_run.log
```

## Cell 10 — Verify shard outputs

Run once the log prints `All shards succeeded.`

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh

n_bin=$(find "$SHARD_OUT_DIR" -name "grm_shard_*_of_${N_SHARDS}.grm.bin.*" | grep -v _test | wc -l)
n_id=$(find "$SHARD_OUT_DIR" -name "grm_shard_*_of_${N_SHARDS}.grm.id" | grep -v _test | wc -l)
echo "grm.bin files: $n_bin / $N_SHARDS"
echo "grm.id  files: $n_id"

TOTAL_KB=$(find "$SHARD_OUT_DIR" -name "grm_shard_*_of_${N_SHARDS}.grm.bin.*" | grep -v _test \
  | xargs du -c | tail -1 | awk '{print $1}')
echo "Total .grm.bin size: $(numfmt --to=iec $((TOTAL_KB * 1024)))"
# Expected: N(N+1)/2 * 4 bytes. N from .fam:
N=$(wc -l < "${BED_PREFIX}.fam")
python3 -c "n=$N; print(f'Expected for N={n}: {n*(n+1)//2*4/1e9:.1f} GB')"

echo "--- .grm.id FID column (first 3 rows) ---"
head -3 "$(find "$SHARD_OUT_DIR" -name 'grm_shard_1_of_'"${N_SHARDS}"'.grm.id' | grep -v _test)"
```

## Cell 11 — Copy results back to bucket

```bash
#!/usr/bin/env bash
source ~/.grm_env.sh
source ~/.grm_run_params.sh
set -e

DEST="${BUCKET_DIR}/grm_shards"
mkdir -p "$DEST"
echo "Copying to $DEST"

n_copied=0
n_skipped=0
for k in $(seq 1 "$N_SHARDS"); do
  src="${SHARD_OUT_DIR}/grm_shard_${k}_of_${N_SHARDS}.grm.bin.${k}"
  dst="${DEST}/grm_shard_${k}_of_${N_SHARDS}.grm.bin.${k}"
  if [ -f "$dst" ]; then
    n_skipped=$((n_skipped + 1))
  else
    cp "$src" "$dst"
    n_copied=$((n_copied + 1))
  fi
done

# one canonical .grm.id (identical across all shards)
cp "${SHARD_OUT_DIR}/grm_shard_1_of_${N_SHARDS}.grm.id" "${DEST}/genome_wide.grm.id"

echo "✓ ${n_copied} new, ${n_skipped} already present, + 1 .grm.id"
gcloud storage du "${BUCKET_DIR}/grm_shards" --readable-sizes 2>/dev/null \
  || du -sh "$DEST"
```

## Cell 12 — Delete the VM

```bash
#!/usr/bin/env bash
gcloud compute instances delete eur-d2-grm-build \
  --zone us-central1-a --project wb-swift-sprout-7231 --quiet
echo "✓ VM deleted"
```

---

## Next Steps

Once shards are in `${BUCKET_DIR}/grm_shards/`:

1. **Verify the `.grm.id` FID column** — `grm_shard_tool` keys on `(FID, IID)`; `write_grm_pheno()` sets `FID = IID = person_id`, but plink often writes `FID = 0`. Check before the phenotype cross-product trusts it.
2. **Run `04_process_shards/notebooks/remote/01_grm_shard_processing.ipynb`** — `grm_shard_tool accumulate` merges the 300 shards into the final GRM.

---

## Bug Log

- Batch (03b) path abandoned for eur_D2: dsub hardcodes `cloud-sdk:294.0.0-slim`
  (purged from gcr.io); patching to a live tag got one SUCCESS but the state
  was too fragile to reproduce reliably. Local 03a is the fallback.
- `PANEL_DIR` for eur_D2 is `03_grm_shards/eur_D2/grm_input/`, not the
  `01_ancestry_filtering/genome_wide_panel_*` location eur_base uses.
- (add more as they arise)
