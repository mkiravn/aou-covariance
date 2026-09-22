# 1kg_eur — covariate PCs

Builds the covariate PCs for the round-2 sample set, using the WGS panel down
to low frequency. Variant QC, pruning and the PCA all run on Batch; the
notebook only reads results.

Three PCAs on the same participants: common variants only, low-frequency only,
and both together. The common one is the default; the other two say whether
low-frequency variants add structure or just artefacts.

Variants come from the unified panel, which was built with a pooled
`--maf 0.001`, so nothing below that is available without going back to raw
ACAF.

---

## Cell 1 — config

```python
import os, shutil, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SAMPLE_SET = "1kg_eur"
CDR_VERSION = "v9"
PROJECT_ID = "wb-swift-sprout-7231"
REGION = "us-central1"
SERVICE_ACCOUNT = "pet-27799165194323faf22e2@wb-swift-sprout-7231.iam.gserviceaccount.com"
NETWORK = f"projects/{PROJECT_ID}/global/networks/network"
SUBNETWORK = f"projects/{PROJECT_ID}/regions/{REGION}/subnetworks/subnetwork"

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
WS_GS = "gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
B, B_GS = f"{WS}/phenotypic_covariance_v9", f"{WS_GS}/phenotypic_covariance_v9"

KEEP = f"{B}/{SAMPLE_SET}/01_ancestry/round2/{SAMPLE_SET}_keep_ids.txt"
KEEP_GS = f"{B_GS}/{SAMPLE_SET}/01_ancestry/round2/{SAMPLE_SET}_keep_ids.txt"
PANEL_GS = f"{B_GS}/01_ancestry_filtering/unified_panel/unified_panel_{CDR_VERSION}"
PLINK2_GS = f"{B_GS}/01_ancestry_filtering/unified_panel/bin/plink2"

OUT = f"{B}/{SAMPLE_SET}/01_ancestry/covariate_pca"
OUT_GS = f"{B_GS}/{SAMPLE_SET}/01_ancestry/covariate_pca"
LOGS_GS = f"{B_GS}/{SAMPLE_SET}/dsub_logs"
LOCAL = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_covpca")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOCAL, exist_ok=True)

KG_DIR = f"{WS}/1000g_reference"
KG_BFILE = f"{KG_DIR}/1kg_all_qc"
KG_PANEL = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"
KG_ACOUNT_GS = f"{WS_GS}/1000g_reference/1kg_all_qc.acount"
EUR_POPS = ["CEU", "GBR", "FIN", "TSI", "IBS"]
ANCHOR_POPS = ["CEU", "GBR"]

# Each arm is a variant set to run the PCA on. hm3: True keeps only HapMap3
# variants, False keeps only the rest, None ignores the split. Add an arm by
# adding an entry.
ARMS = {
    "hm3":        {"maf": (0.01, 0.5),   "hm3": True},
    "panel":      {"maf": (0.01, 0.5),   "hm3": None},
    "non_hm3":    {"maf": (0.01, 0.5),   "hm3": False},
    # "lowfreq":  {"maf": (0.001, 0.01), "hm3": None},
}
THIN_TO = None                 # None keeps every pruned variant; set a number to cut
                               # each arm to it, which makes the arms differ only in
                               # which variants they hold, not how many
PRUNE = "1000kb 1 0.05"        # kb window: plink2 requires the step to be 1
PRUNE_PASSES = 2               # a second pass catches LD the first window missed
N_PCS_FIT = 20
N_PCS_COVARIATE = 20
SEED = 1

QC_MACHINE, QC_MEM_MB, QC_DISK = "n1-highmem-16", 96_000, 500
PCA_MACHINE, PCA_MEM_MB, PCA_DISK = "n1-highmem-32", 200_000, 500

assert os.path.isfile(KEEP), f"{KEEP} -- run round 2 first"
print(f"{sum(1 for _ in open(KEEP)):,} participants -> {OUT}")
```

## Cell 2 — dsub

`dsub` pins a purged cloud-sdk tag; patch it and submit in the same shell.

```python
%%bash
set -e
pip install --quiet --upgrade 'dsub>=0.5.3'
DSUB_DIR=$(python -c "import dsub, os; print(os.path.dirname(dsub.__file__))")
sed -i -E "s|cloud-sdk:[0-9]+\.[0-9]+\.[0-9]+-slim|cloud-sdk:581.0.0-slim|g" \
  "${DSUB_DIR}/providers/google_utils.py"
dsub --version
grep -o "cloud-sdk:[^'\"]*" "${DSUB_DIR}/providers/google_utils.py" | sort -u
```

## Cell 2b — plink2 on this VM

Only Cell 7 needs it, for the 1000G projection.

```python
%%bash
set -e
BIN_DIR="$HOME/bin"
mkdir -p "$BIN_DIR"
if [ ! -x "$BIN_DIR/plink2" ]; then
  cd /tmp
  wget -q -O plink2.zip "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
  unzip -o -q plink2.zip plink2 -d "$BIN_DIR"
  chmod +x "$BIN_DIR/plink2"
fi
"$BIN_DIR/plink2" --version
```

```python
os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"
```

## Cell 3 — long-range LD regions

```python
HIGH_LD_REGIONS_GRCH38 = """\
chr1 47761740 51761740 1
chr2 85919365 100517106 2
chr2 89917298 89917322 3
chr2 87416141 87416186 4
chr2 87417804 87417863 5
chr2 87418924 87418981 6
chr1 144106678 144106709 7
chr14 87391719 87391996 8
chr9 40365644 40365693 9
chr2 182427027 189427029 10
chr3 47483505 49987563 11
chr3 83368158 86868160 12
chr5 44464140 51168409 13
chr5 129636407 132636409 14
chr6 25391792 33424245 15
chr6 26726947 26726981 16
chr6 57788603 58453888 17
chr6 61109122 61357029 18
chr6 61424410 61424451 19
chr9 64198500 64200392 20
chr6 139637169 142137170 21
chr7 54964812 66897578 22
chr7 62182500 62277073 23
chr12 34639034 34639084 24
chr14 94658026 94658080 25
chr17 43159541 43159574 26
chr2 135275091 135275210 27
chr1 181955019 181955047 28
chr22 30060084 30060162 29
chr9 88958735 88959017 30
chr2 207609786 207609808 31
chr22 42980497 42980522 32
chr20 4031884 4032441 33
chr8 8105067 12105082 34
chr8 43025699 48924888 35
chr8 47303500 47317337 36
chr8 110918594 113918595 37
chr10 36671065 43184546 38
chr10 41693521 41885273 39
chr1 125169943 125170022 40
chr11 88127183 91127184 41
chr12 32955798 41319931 42
chr20 33948532 36438183 43
"""
PEAKS = f"{OUT}/ld_peaks.txt"            # written by Cell 7b, empty on the first pass
LD_REGIONS = f"{OUT}/exclude_regions.txt"
open(LD_REGIONS, "w").write(HIGH_LD_REGIONS_GRCH38
                            + (open(PEAKS).read() if os.path.isfile(PEAKS) else ""))
LD_REGIONS_GS = f"{OUT_GS}/exclude_regions.txt"
print(f"{sum(1 for _ in open(LD_REGIONS))} regions excluded -> {LD_REGIONS}")
```

## Cell 4 — QC and prune each arm (Batch)

QC once, then per arm: apply the frequency range, keep or drop the HapMap3
variants, remove the long-range LD regions, and prune. Thinning happens later,
once the pruned counts are known, so every arm can be cut to the same size.
`bed1` is the 1-based interval format the region file uses.

```python
MAF_FLOOR = min(a["maf"][0] for a in ARMS.values())
arm_case = "\n".join(
    f'        {name}) FLOOR={a["maf"][0]}; CEIL={a["maf"][1]}; '
    f'HM3={"extract" if a["hm3"] is True else "exclude" if a["hm3"] is False else "none"} ;;'
    for name, a in ARMS.items())

cmd = f"""
    set -e
    chmod +x "$PLINK_BIN"
    Q="${{TMPDIR:-/tmp}}/qc"

    "$PLINK_BIN" --pgen "$PANEL_PGEN" --pvar "$PANEL_PVAR" --psam "$PANEL_PSAM" \
      --keep "$KEEP_PATH" --nonfounders \
      --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all \
      --maf {MAF_FLOOR} --geno 0.01 --hwe 1e-6 0 keep-fewhet \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "$Q"
    echo "after QC: $(grep -vc '^##' "${{Q}}.pvar") variants"

    # HapMap3 membership: ID+REF+ALT agreement with the 1000G reference
    grep -v '^##' "${{Q}}.pvar" | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > "${{Q}}_ira.sorted"
    awk 'NR>1 {{print $2, $3, $4}}' "$KG_ACOUNT" | LC_ALL=C sort > "${{Q}}_kg.sorted"
    LC_ALL=C comm -12 "${{Q}}_ira.sorted" "${{Q}}_kg.sorted" | awk '{{print $1}}' > "${{OUT_DIR}}/hm3.ids"
    echo "HapMap3-agreeing: $(wc -l < "${{OUT_DIR}}/hm3.ids")"

    for ARM in {" ".join(ARMS)}; do
      case "$ARM" in
{arm_case}
      esac
      case "$HM3" in
        extract) SEL="--extract ${{OUT_DIR}}/hm3.ids" ;;
        exclude) SEL="--exclude ${{OUT_DIR}}/hm3.ids" ;;
        *)       SEL="" ;;
      esac

      "$PLINK_BIN" --pfile "$Q" --maf "$FLOOR" --max-maf "$CEIL" $SEL \
        --exclude bed1 "$LD_REGIONS" --nonfounders \
        --indep-pairwise {PRUNE} \
        --threads "$VCPUS" --memory "$MEM_MB" --out "${{TMPDIR:-/tmp}}/${{ARM}}"
      echo "$ARM pruned, pass 1: $(wc -l < "${{TMPDIR:-/tmp}}/${{ARM}}.prune.in")"

      for PASS in $(seq 2 {PRUNE_PASSES}); do
        "$PLINK_BIN" --pfile "$Q" --extract "${{TMPDIR:-/tmp}}/${{ARM}}.prune.in" --nonfounders \
          --indep-pairwise {PRUNE} \
          --threads "$VCPUS" --memory "$MEM_MB" --out "${{TMPDIR:-/tmp}}/${{ARM}}_p${{PASS}}"
        mv "${{TMPDIR:-/tmp}}/${{ARM}}_p${{PASS}}.prune.in" "${{TMPDIR:-/tmp}}/${{ARM}}.prune.in"
        echo "$ARM pruned, pass $PASS: $(wc -l < "${{TMPDIR:-/tmp}}/${{ARM}}.prune.in")"
      done
      cp "${{TMPDIR:-/tmp}}/${{ARM}}.prune.in" "${{OUT_DIR}}/${{ARM}}.prune.in"
    done

    cp "${{Q}}.pgen" "${{OUT_DIR}}/qc.pgen"
    cp "${{Q}}.pvar" "${{OUT_DIR}}/qc.pvar"
    cp "${{Q}}.psam" "${{OUT_DIR}}/qc.psam"
"""
open(f"{LOCAL}/qc_cmd.sh", "w").write(cmd)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "covpca-prune-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {QC_MACHINE} --disk-size {QC_DISK} \
  --input PANEL_PGEN="{PANEL_GS}.pgen" \
  --input PANEL_PVAR="{PANEL_GS}.pvar" \
  --input PANEL_PSAM="{PANEL_GS}.psam" \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input KEEP_PATH="{KEEP_GS}" \
  --input LD_REGIONS="{LD_REGIONS_GS}" \
  --input KG_ACOUNT="{KG_ACOUNT_GS}" \
  --env VCPUS={QC_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={QC_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/panels" \
  --script {LOCAL}/qc_cmd.sh
"""], capture_output=True, text=True)
print(job.stdout or job.stderr)
QC_JOB = job.stdout.strip().splitlines()[-1] if job.stdout else None
```

## Cell 5 — watch it

```python
print(subprocess.run(["bash", "-c",
    f"dstat --provider google-batch --project {PROJECT_ID} --location {REGION} "
    f"--jobs {QC_JOB} --users '*' --status '*' --full"],
    capture_output=True, text=True).stdout[-3000:])
```

## Cell 5b — pruned counts

Each arm keeps what pruning left it. Set `THIN_TO` in Cell 1 if you would
rather every arm carry the same number of variants, which separates "which
variants" from "how many" at the cost of discarding some.

```python
counts = {a: sum(1 for _ in open(f"{OUT}/panels/{a}.prune.in")) for a in ARMS}
C_ = pd.Series(counts, name="pruned").to_frame()
C_["used"] = [min(n, THIN_TO) if THIN_TO else n for n in C_["pruned"]]
print(C_.to_string())
if THIN_TO is None and max(counts.values()) > 1.5 * min(counts.values()):
    print("\nthe arms differ in size by more than half; PC precision rises with marker "
          "count, so read a difference between them with that in mind")
```

## Cell 6 — thin and fit the PCA per arm (Batch)

One task per arm: fit the PCA with loadings on whatever that arm kept (thinned
first only if `THIN_TO` is set), then re-score the same participants through
those loadings so every arm's scores sit on `--score` coordinates.

```python
TASKS = f"{LOCAL}/pca_tasks.tsv"
with open(TASKS, "w") as f:
    f.write("--env ARM\n")
    for a in ARMS:
        f.write(f"{a}\n")

pca_cmd = f"""
    set -e
    chmod +x "$PLINK_BIN"
    P="${{TMPDIR:-/tmp}}/${{ARM}}"

    THIN=""
    if [ -n "$THIN_TO" ] && [ "$(wc -l < "$PRUNE_IN")" -gt "$THIN_TO" ]; then
      THIN="--thin-count $THIN_TO"
    fi

    "$PLINK_BIN" --pgen "$QC_PGEN" --pvar "$QC_PVAR" --psam "$QC_PSAM" \
      --extract "$PRUNE_IN" $THIN --seed {SEED} --nonfounders \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "$P"
    echo "$ARM variants: $(grep -vc '^##' "${{P}}.pvar")"

    "$PLINK_BIN" --pfile "$P" --nonfounders --freq counts \
      --pca approx {N_PCS_FIT} allele-wts \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${{OUT_DIR}}/${{ARM}}"

    W="${{OUT_DIR}}/${{ARM}}.eigenvec.allele"
    HEADER=$(head -1 "$W")
    ID=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'ID' | cut -d: -f1)
    A1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'A1' | cut -d: -f1)
    P1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC1' | cut -d: -f1)
    PK=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)

    "$PLINK_BIN" --pfile "$P" --nonfounders \
      --read-freq "${{OUT_DIR}}/${{ARM}}.acount" \
      --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
      --score-col-nums "${{P1}}-${{PK}}" \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${{OUT_DIR}}/${{ARM}}_scores"

    cp "${{P}}.pvar" "${{OUT_DIR}}/${{ARM}}_variants.pvar"
"""
open(f"{LOCAL}/pca_cmd.sh", "w").write(pca_cmd)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "covpca-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {PCA_MACHINE} --disk-size {PCA_DISK} \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input QC_PGEN="{OUT_GS}/panels/qc.pgen" \
  --input QC_PVAR="{OUT_GS}/panels/qc.pvar" \
  --input QC_PSAM="{OUT_GS}/panels/qc.psam" \
  --input PRUNE_IN="{OUT_GS}/panels/${{ARM}}.prune.in" \
  --env VCPUS={PCA_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={PCA_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/pca" \
  --tasks {TASKS} \
  --script {LOCAL}/pca_cmd.sh
"""], capture_output=True, text=True)
print(job.stdout or job.stderr)
PCA_JOB = job.stdout.strip().splitlines()[-1] if job.stdout else None
```

## Cell 7 — read each arm, project 1000G

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]
ARM_LIST = list(ARMS)


def read_scores(path, id_name):
    d = pd.read_csv(path, sep=r"\s+")
    d = d.rename(columns={("#IID" if "#IID" in d.columns else "IID"): id_name})
    d = d.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS_FIT + 1)})
    d[id_name] = d[id_name].astype(str)
    return d[[id_name] + PC]


scores = {a: read_scores(f"{OUT}/pca/{a}_scores.sscore", "person_id") for a in ARM_LIST}
eigen = {a: np.loadtxt(f"{OUT}/pca/{a}.eigenval") for a in ARM_LIST}

kg_scores = {}
for a in ARM_LIST:
    sh_ = f"""
    set -eo pipefail
    grep -v '^##' "{OUT}/pca/{a}_variants.pvar" | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > "{LOCAL}/{a}_ira.sorted"
    awk 'NR>1 {{print $2, $3, $4}}' "{KG_BFILE}.acount" | LC_ALL=C sort > "{LOCAL}/kg_ira.sorted"
    LC_ALL=C comm -12 "{LOCAL}/{a}_ira.sorted" "{LOCAL}/kg_ira.sorted" | awk '{{print $1}}' > "{LOCAL}/{a}_kg.ids"
    N=$(wc -l < "{LOCAL}/{a}_kg.ids")
    echo "{a}: $N variants shared with 1000G"
    if [ "$N" -lt 1000 ]; then exit 0; fi
    W="{OUT}/pca/{a}.eigenvec.allele"
    HEADER=$(head -1 "$W")
    ID=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'ID' | cut -d: -f1)
    A1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'A1' | cut -d: -f1)
    P1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC1' | cut -d: -f1)
    PK=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)
    plink2 --bfile "{KG_BFILE}" --extract "{LOCAL}/{a}_kg.ids" --nonfounders \
      --read-freq "{OUT}/pca/{a}.acount" \
      --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
      --score-col-nums "${{P1}}-${{PK}}" --out "{LOCAL}/kg_in_{a}"
    """
    subprocess.run(["bash", "-c", sh_], check=True)
    if os.path.isfile(f"{LOCAL}/kg_in_{a}.sscore"):
        kg_scores[a] = read_scores(f"{LOCAL}/kg_in_{a}.sscore", "sample").merge(
            pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]],
            on="sample", how="left")
    else:
        print(f"  {a}: too few shared variants to project 1000G")
```

The non-HapMap3 arm shares little with the 1000G reference by construction, so
its projection may be skipped. That is expected: it can describe structure, but
it cannot carry the anchor.

## Cell 8 — do the arms agree?

For each arm, every PC is matched to its best counterpart in the reference arm
by absolute correlation. A PC matching above about 0.95 carries nothing new; a
low match is either real extra structure or noise, which the next two cells
separate.

```python
REF_ARM = ARM_LIST[0]
K_CMP = min(10, N_PCS_FIT)

best = {}
for a in ARM_LIST:
    if a == REF_ARM:
        continue
    j = scores[REF_ARM].merge(scores[a], on="person_id", suffixes=("_ref", "_arm"))
    M = np.abs(np.corrcoef(j[[f"{p}_arm" for p in PC[:K_CMP]]].to_numpy(),
                           j[[f"{p}_ref" for p in PC[:K_CMP]]].to_numpy(),
                           rowvar=False)[:K_CMP, K_CMP:])
    best[a] = pd.Series(M.max(1), index=PC[:K_CMP])
B = pd.DataFrame(best)
B.to_csv(f"{OUT}/arm_pc_agreement.tsv", sep="\t")
print("best |r| with any PC of the reference arm")
print(B.round(3).to_string())

sep = {}
for a in ARM_LIST:
    if a not in kg_scores:
        continue
    kgs = kg_scores[a]
    anch = kgs["pop"].isin(ANCHOR_POPS)
    other = kgs["pop"].isin(set(EUR_POPS) - set(ANCHOR_POPS))
    sep[a] = pd.Series(
        [abs(kgs.loc[anch, p].mean() - kgs.loc[other, p].mean()) / scores[a][p].std()
         for p in PC[:K_CMP]], index=PC[:K_CMP])
S_ = pd.DataFrame(sep)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
for a in ARM_LIST:
    ev = eigen[a]
    axes[0].plot(range(1, len(ev) + 1), ev / ev.sum() * 100, marker="o", ms=3, label=a)
axes[0].set_xlabel("PC"); axes[0].set_ylabel("% variance"); axes[0].legend(fontsize=8)
axes[0].set_title("scree")

for a in B:
    axes[1].plot(range(1, K_CMP + 1), B[a], marker="o", ms=3, label=a)
axes[1].axhline(0.95, color="0.7", ls=":", lw=1)
axes[1].set_ylim(0, 1.02)
axes[1].set_xlabel("PC"); axes[1].set_ylabel(f"best |r| with {REF_ARM}")
axes[1].legend(fontsize=8)
axes[1].set_title("agreement with the reference arm")

for a in S_:
    axes[2].plot(range(1, K_CMP + 1), S_[a], marker="o", ms=3, label=a)
axes[2].set_xlabel("PC"); axes[2].set_ylabel("anchor separation (participant SD)")
axes[2].legend(fontsize=8)
axes[2].set_title("CEU + GBR vs other European populations")
fig.suptitle(f"{SAMPLE_SET} covariate PCA — "
             + ", ".join(f"{a} ({C_.loc[a, 'used']:,})" for a in ARM_LIST))
plt.tight_layout()
plt.savefig(f"{OUT}/arm_comparison.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 9 — LD peaks in the loadings

A PC carried by one region is LD, not structure. This flags variants whose
squared loading is far above the rest, widens each to a window, and writes them
to `ld_peaks.txt`. If anything is flagged, rerun Cells 3, 4 and 6: the exclusion
file is the long-range regions plus these, so the next fit drops them.

```python
PEAK_Q, FLANK_KB, PEAK_PCS = 0.9995, 250, min(10, N_PCS_FIT)

found = []
for a in ARM_LIST:
    L = pd.read_csv(f"{OUT}/pca/{a}.eigenvec.allele", sep=r"\s+")
    idc = "#ID" if "#ID" in L.columns else "ID"
    L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
    L["POS"] = L["POS"].astype(int)
    hit = np.zeros(len(L), dtype=bool)
    for p in PC[:PEAK_PCS]:
        v = L[p].to_numpy() ** 2
        hit |= v > np.quantile(v, PEAK_Q)
    if hit.any():
        h = L.loc[hit, ["CHROM", "POS"]].copy()
        h["arm"] = a
        found.append(h)
    print(f"{a}: {int(hit.sum())} loading outliers over PC1-{PEAK_PCS}")

if found:
    H = pd.concat(found).sort_values(["CHROM", "POS"])
    regions, k = [], 0
    for chrom, sub in H.groupby("CHROM", sort=False):
        lo = hi = None
        for pos in sub["POS"]:
            if lo is None or pos > hi + FLANK_KB * 1000:
                if lo is not None:
                    k += 1
                    regions.append(f"{chrom} {max(lo - FLANK_KB * 1000, 1)} {hi + FLANK_KB * 1000} peak{k}")
                lo = pos
            hi = pos
        k += 1
        regions.append(f"{chrom} {max(lo - FLANK_KB * 1000, 1)} {hi + FLANK_KB * 1000} peak{k}")
    open(PEAKS, "w").write("\n".join(regions) + "\n")
    print(f"\n{len(regions)} peak regions -> {PEAKS}; rerun Cells 3, 4 and 6")
    print("\n".join(regions[:10]))
else:
    print("\nno peaks; nothing to exclude")
```

## Cell 10 — is any of it technical?

Set `META` to a sample-level metadata table if one exists; the check is skipped
otherwise. R² of a PC on a batch or site factor says whether it tracks
sequencing rather than ancestry.

```python
META = None            # e.g. f"{AUX}/qc/<file>.tsv"
META_ID, META_FACTORS = "research_id", ["site_id", "sequencing_center"]

if META and os.path.isfile(META):
    m = pd.read_csv(META, sep="\t", dtype=str)
    m = m.rename(columns={META_ID: "person_id"})
    for a in ARM_LIST:
        j = scores[a].merge(m, on="person_id", how="inner")
        out = {}
        for f in META_FACTORS:
            if f not in j:
                continue
            g = j.groupby(f)
            out[f] = [float(1 - g[p].transform("var").mean() / j[p].var()) for p in PC[:5]]
        print(f"\n{a}: R² of PC1-5 on each factor")
        print(pd.DataFrame(out, index=PC[:5]).round(3).to_string())
else:
    print("no metadata table set; skipping the batch check")
```

## Cell 11 — write the covariate PCs

```python
ARM_FOR_COVARIATES = ARM_LIST[0]      # change once Cells 8-10 justify it

cov = scores[ARM_FOR_COVARIATES][["person_id"] + PC[:N_PCS_COVARIATE]].copy()
cov = cov.rename(columns={"person_id": "IID"})
COV_PATH = f"{OUT}/final_pca_pc_covariates_{SAMPLE_SET}.txt"
cov.to_csv(COV_PATH, sep="\t", index=False)

with open(f"{OUT}/covariate_pca_provenance.txt", "w") as f:
    f.write(f"panel\tunified_panel_{CDR_VERSION}, pooled --maf 0.001 at build time\n"
            f"qc\t--snps-only just-acgt --max-alleles 2 --rm-dup exclude-all "
            f"--maf {MAF_FLOOR} --geno 0.01 --hwe 1e-6 0 keep-fewhet --nonfounders\n"
            f"arms\t{ARMS}\n"
            f"prune\t--indep-pairwise {PRUNE}, {PRUNE_PASSES} passes, "
            f"long-range LD regions and detected peaks excluded (bed1)\n"
            f"thin\t{'--thin-count ' + str(THIN_TO) if THIN_TO else 'none; every pruned variant kept'}\n"
            f"variants\t{dict(C_['used'])}\n"
            f"pca\t--pca approx {N_PCS_FIT} allele-wts, fit on the participants\n"
            f"scores\t--score through the fit's own loadings, --read-freq its .acount\n"
            f"covariates\tarm={ARM_FOR_COVARIATES}, PC1-{N_PCS_COVARIATE}\n")
print(f"{len(cov):,} participants, PC1-{N_PCS_COVARIATE} -> {COV_PATH}")
```

## Cell 12 — copy this notebook to the bucket

Save the notebook first, then run this.

```python
NOTEBOOK = os.path.expanduser(f"~/notebooks/{SAMPLE_SET}_covariate_pca.ipynb")
NB_DIR = f"{B}/{SAMPLE_SET}/notebooks"
os.makedirs(NB_DIR, exist_ok=True)
assert os.path.isfile(NOTEBOOK), NOTEBOOK
dest = f"{NB_DIR}/{os.path.basename(NOTEBOOK)}"
shutil.copy(NOTEBOOK, dest)
print(f"{'OK  ' if os.path.getsize(dest) == os.path.getsize(NOTEBOOK) else 'SIZE MISMATCH'} {dest}")
```

---

Next: residualization reads `final_pca_pc_covariates_{SAMPLE_SET}.txt`.
