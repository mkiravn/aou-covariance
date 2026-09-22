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
EUR_POPS = ["CEU", "GBR", "FIN", "TSI", "IBS"]
ANCHOR_POPS = ["CEU", "GBR"]

BANDS = {"common": (0.01, 0.5), "lowfreq": (0.001, 0.01)}   # MAF floor, ceiling
TARGET_PER_BAND = 200_000
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

## Cell 4 — QC, prune and thin (Batch)

One job. QC once, then per band: restrict the frequency range, drop the
long-range LD regions, prune, and thin to the target. `bed1` is the 1-based
interval format the region file above uses (`range` is its old alias).

```python
MAF_FLOOR = min(lo for lo, _ in BANDS.values())
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

    for BAND in {" ".join(BANDS)}; do
      case "$BAND" in
        common)  FLOOR={BANDS['common'][0]};  CEIL={BANDS['common'][1]} ;;
        lowfreq) FLOOR={BANDS['lowfreq'][0]}; CEIL={BANDS['lowfreq'][1]} ;;
      esac
      "$PLINK_BIN" --pfile "$Q" --maf "$FLOOR" --max-maf "$CEIL" \
        --exclude bed1 "$LD_REGIONS" --nonfounders \
        --indep-pairwise {PRUNE} \
        --threads "$VCPUS" --memory "$MEM_MB" --out "${{Q}}_${{BAND}}"
      echo "$BAND pruned, pass 1: $(wc -l < "${{Q}}_${{BAND}}.prune.in")"

      for PASS in $(seq 2 {PRUNE_PASSES}); do
        "$PLINK_BIN" --pfile "$Q" --extract "${{Q}}_${{BAND}}.prune.in" --nonfounders \
          --indep-pairwise {PRUNE} \
          --threads "$VCPUS" --memory "$MEM_MB" --out "${{Q}}_${{BAND}}_p${{PASS}}"
        mv "${{Q}}_${{BAND}}_p${{PASS}}.prune.in" "${{Q}}_${{BAND}}.prune.in"
        echo "$BAND pruned, pass $PASS: $(wc -l < "${{Q}}_${{BAND}}.prune.in")"
      done
      NPRUNED=$(wc -l < "${{Q}}_${{BAND}}.prune.in")
      THIN=""
      if [ "$NPRUNED" -gt {TARGET_PER_BAND} ]; then THIN="--thin-count {TARGET_PER_BAND}"; fi

      "$PLINK_BIN" --pfile "$Q" --extract "${{Q}}_${{BAND}}.prune.in" \
        $THIN --seed {SEED} --nonfounders \
        --threads "$VCPUS" --memory "$MEM_MB" --make-pgen \
        --out "${{OUT_DIR}}/pca_${{BAND}}"
      echo "$BAND kept: $(grep -vc '^##' "${{OUT_DIR}}/pca_${{BAND}}.pvar")"
    done

    cat "${{OUT_DIR}}/pca_common.pvar" "${{OUT_DIR}}/pca_lowfreq.pvar" \
      | awk '$1 !~ /^#/ {{print $3}}' | sort -u > "${{TMPDIR:-/tmp}}/both.ids"
    "$PLINK_BIN" --pfile "$Q" --extract "${{TMPDIR:-/tmp}}/both.ids" --nonfounders \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "${{OUT_DIR}}/pca_both"
    echo "both kept: $(grep -vc '^##' "${{OUT_DIR}}/pca_both.pvar")"
"""
open(f"{LOCAL}/qc_cmd.sh", "w").write(cmd)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "covpca-qc-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {QC_MACHINE} --disk-size {QC_DISK} \
  --input PANEL_PGEN="{PANEL_GS}.pgen" \
  --input PANEL_PVAR="{PANEL_GS}.pvar" \
  --input PANEL_PSAM="{PANEL_GS}.psam" \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input KEEP_PATH="{KEEP_GS}" \
  --input LD_REGIONS="{LD_REGIONS_GS}" \
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

## Cell 6 — PCA per band (Batch)

One task per band. Each fits the PCA, writes loadings and allele counts, then
re-scores the same participants through those loadings, so every set of scores
sits on `--score` coordinates and is directly comparable.

```python
TASKS = f"{LOCAL}/pca_tasks.tsv"
with open(TASKS, "w") as f:
    f.write("--env BAND\t--input PGEN\t--input PVAR\t--input PSAM\n")
    for band in list(BANDS) + ["both"]:
        f.write(f"{band}\t{OUT_GS}/panels/pca_{band}.pgen\t"
                f"{OUT_GS}/panels/pca_{band}.pvar\t{OUT_GS}/panels/pca_{band}.psam\n")

pca_cmd = f"""
    set -e
    chmod +x "$PLINK_BIN"
    P="${{TMPDIR:-/tmp}}/in"
    cp "$PGEN" "${{P}}.pgen"; cp "$PVAR" "${{P}}.pvar"; cp "$PSAM" "${{P}}.psam"

    "$PLINK_BIN" --pfile "$P" --nonfounders --freq counts \
      --pca approx {N_PCS_FIT} allele-wts \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${{OUT_DIR}}/${{BAND}}"

    W="${{OUT_DIR}}/${{BAND}}.eigenvec.allele"
    HEADER=$(head -1 "$W")
    ID=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'ID' | cut -d: -f1)
    A1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'A1' | cut -d: -f1)
    P1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC1' | cut -d: -f1)
    PK=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)

    "$PLINK_BIN" --pfile "$P" --nonfounders \
      --read-freq "${{OUT_DIR}}/${{BAND}}.acount" \
      --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
      --score-col-nums "${{P1}}-${{PK}}" \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${{OUT_DIR}}/${{BAND}}_scores"
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
  --env VCPUS={PCA_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={PCA_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/pca" \
  --tasks {TASKS} \
  --script {LOCAL}/pca_cmd.sh
"""], capture_output=True, text=True)
print(job.stdout or job.stderr)
PCA_JOB = job.stdout.strip().splitlines()[-1] if job.stdout else None
```

## Cell 7 — read the three PCAs, project 1000G

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]
BAND_LIST = list(BANDS) + ["both"]


def read_scores(path, id_name):
    d = pd.read_csv(path, sep=r"\s+")
    d = d.rename(columns={("#IID" if "#IID" in d.columns else "IID"): id_name})
    d = d.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS_FIT + 1)})
    d[id_name] = d[id_name].astype(str)
    return d[[id_name] + PC]


scores = {b: read_scores(f"{OUT}/pca/{b}_scores.sscore", "person_id") for b in BAND_LIST}
eigen = {b: np.loadtxt(f"{OUT}/pca/{b}.eigenval") for b in BAND_LIST}

kg_scores = {}
for b in BAND_LIST:
    sh = f"""
    set -eo pipefail
    grep -v '^##' "{OUT}/panels/pca_{b}.pvar" | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > "{LOCAL}/{b}_ira.sorted"
    awk 'NR>1 {{print $2, $3, $4}}' "{KG_BFILE}.acount" | LC_ALL=C sort > "{LOCAL}/kg_ira.sorted"
    LC_ALL=C comm -12 "{LOCAL}/{b}_ira.sorted" "{LOCAL}/kg_ira.sorted" | awk '{{print $1}}' > "{LOCAL}/{b}_kg.ids"
    W="{OUT}/pca/{b}.eigenvec.allele"
    HEADER=$(head -1 "$W")
    ID=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'ID' | cut -d: -f1)
    A1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'A1' | cut -d: -f1)
    P1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC1' | cut -d: -f1)
    PK=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)
    plink2 --bfile "{KG_BFILE}" --extract "{LOCAL}/{b}_kg.ids" --nonfounders \
      --read-freq "{OUT}/pca/{b}.acount" \
      --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
      --score-col-nums "${{P1}}-${{PK}}" --out "{LOCAL}/kg_in_{b}"
    """
    subprocess.run(["bash", "-c", sh], check=True)
    kg_scores[b] = read_scores(f"{LOCAL}/kg_in_{b}.sscore", "sample").merge(
        pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]], on="sample", how="left")
    print(f"{b}: {sum(1 for _ in open(f'{LOCAL}/{b}_kg.ids')):,} variants for the 1000G projection")
```

## Cell 8 — do the low-frequency PCs add anything?

`anchor_sep` is how far CEU + GBR sit from the other European populations on
each PC, in participant SDs. `max_|r|` is the strongest correlation of that PC
with any common-band PC: near 1 means it is already covered.

```python
common = scores["common"]
rows = []
for b in BAND_LIST:
    s, kg = scores[b], kg_scores[b]
    for k, p in enumerate(PC, 1):
        a = kg[kg["pop"].isin(ANCHOR_POPS)][p].mean()
        o = kg[kg["pop"].isin(set(EUR_POPS) - set(ANCHOR_POPS))][p].mean()
        r = (np.abs(np.corrcoef(s[p], common[PC].to_numpy().T)[0, 1:]).max()
             if b != "common" else np.nan)
        rows.append({"band": b, "PC": p, "pct_var": eigen[b][k - 1] / eigen[b].sum() * 100,
                     "anchor_sep": abs(a - o) / s[p].std(), "max_|r| vs common": r})
D = pd.DataFrame(rows)
D.to_csv(f"{OUT}/pc_comparison.tsv", sep="\t", index=False)
print(D.pivot(index="PC", columns="band", values="anchor_sep").reindex(PC).round(2).to_string())

fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
for b in BAND_LIST:
    d = D[D["band"] == b]
    axes[0].plot(range(1, N_PCS_FIT + 1), d["pct_var"], marker="o", ms=3, label=b)
    axes[1].plot(range(1, N_PCS_FIT + 1), d["anchor_sep"], marker="o", ms=3, label=b)
    if b != "common":
        axes[2].plot(range(1, N_PCS_FIT + 1), d["max_|r| vs common"], marker="o", ms=3, label=b)
axes[0].set_ylabel("% variance"); axes[1].set_ylabel("anchor separation (SD)")
axes[2].set_ylabel("max |r| with a common-band PC"); axes[2].set_ylim(0, 1)
for ax in axes:
    ax.set_xlabel("PC"); ax.legend(fontsize=8)
fig.suptitle(f"{SAMPLE_SET} covariate PCs — common vs low-frequency")
plt.tight_layout()
plt.savefig(f"{OUT}/band_comparison.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 7b — LD peaks in the loadings

A PC carried by one region is LD, not structure. This flags variants whose
squared loading is far above the rest, widens each to a window, and writes them
to `ld_peaks.txt`. If anything is flagged, rerun Cells 3, 4 and 6: the exclusion
file is the long-range regions plus these, so the next fit drops them.

```python
PEAK_Q, FLANK_KB, PEAK_PCS = 0.9995, 250, min(10, N_PCS_FIT)

found = []
for b in BAND_LIST:
    L = pd.read_csv(f"{OUT}/pca/{b}.eigenvec.allele", sep=r"\s+")
    idc = "#ID" if "#ID" in L.columns else "ID"
    L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
    L["POS"] = L["POS"].astype(int)
    hit = np.zeros(len(L), dtype=bool)
    for p in PC[:PEAK_PCS]:
        v = L[p].to_numpy() ** 2
        hit |= v > np.quantile(v, PEAK_Q)
    if hit.any():
        h = L.loc[hit, ["CHROM", "POS"]].copy()
        h["band"] = b
        found.append(h)
    print(f"{b}: {int(hit.sum())} loading outliers over PC1-{PEAK_PCS}")

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

## Cell 9 — are the low-frequency PCs technical?

Set `META` to a sample-level metadata table if one exists; the check is skipped
otherwise. R² of each PC on a batch or site factor says whether it is tracking
sequencing rather than ancestry.

```python
META = None            # e.g. f"{AUX}/qc/<file>.tsv"
META_ID, META_FACTORS = "research_id", ["site_id", "sequencing_center"]

if META and os.path.isfile(META):
    m = pd.read_csv(META, sep="\t", dtype=str)
    m = m.rename(columns={META_ID: "person_id"})
    for b in BAND_LIST:
        j = scores[b].merge(m, on="person_id", how="inner")
        out = {}
        for f in META_FACTORS:
            if f not in j:
                continue
            g = j.groupby(f)
            out[f] = [float(1 - g[p].transform("var").mean() / j[p].var()) for p in PC[:5]]
        print(f"\n{b}: R² of PC1-5 on each factor")
        print(pd.DataFrame(out, index=PC[:5]).round(3).to_string())
else:
    print("no metadata table set; skipping the batch check")
```

## Cell 10 — write the covariate PCs

```python
BAND_FOR_COVARIATES = "common"      # change only if Cell 8 and 9 justify it

cov = scores[BAND_FOR_COVARIATES][["person_id"] + PC[:N_PCS_COVARIATE]].copy()
cov = cov.rename(columns={"person_id": "IID"})
COV_PATH = f"{OUT}/final_pca_pc_covariates_{SAMPLE_SET}.txt"
cov.to_csv(COV_PATH, sep="\t", index=False)

with open(f"{OUT}/covariate_pca_provenance.txt", "w") as f:
    f.write(f"panel\tunified_panel_{CDR_VERSION}, pooled --maf 0.001 at build time\n"
            f"qc\t--snps-only just-acgt --max-alleles 2 --rm-dup exclude-all "
            f"--maf {MAF_FLOOR} --geno 0.01 --hwe 1e-6 0 keep-fewhet --nonfounders\n"
            f"bands\t{BANDS}\n"
            f"prune\t--indep-pairwise {PRUNE}, long-range LD regions excluded (bed1)\n"
            f"thin\t--thin-count {TARGET_PER_BAND} --seed {SEED} per band\n"
            f"pca\t--pca approx {N_PCS_FIT} allele-wts, fit on the participants\n"
            f"scores\t--score through the fit's own loadings, --read-freq its .acount\n"
            f"covariates\tband={BAND_FOR_COVARIATES}, PC1-{N_PCS_COVARIATE}\n")
print(f"{len(cov):,} participants, PC1-{N_PCS_COVARIATE} -> {COV_PATH}")
```

## Cell 11 — copy this notebook to the bucket

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
