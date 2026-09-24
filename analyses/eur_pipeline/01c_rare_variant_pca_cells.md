# 1kg_eur — rare variant PCA

Explores rare variants (0.1% ≤ MAF < 1%) in the r2 EUR cohort, split into
HM3-agreeing and non-HM3, to characterise their structure and contribution to
the covariate PCA.

HWE filtering is skipped: at MAF < 1% the test is underpowered and
preferentially drops real rare heterozygotes. 1000G projection is not attempted:
most AoU-rare variants are monomorphic or near-absent in the 2,500-sample
1000G reference.

Source: `unified_panel_v9` (MAF floor 0.1%), which contains variants below the
`r1_qc` MAF cutoff. QC and HM3 splitting run as a Batch job; pruning and PCA
run locally.

---

## Cell 1 — config

```python
import os, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SAMPLE_SET  = "1kg_eur"
CDR_VERSION = "v9"
PROJECT_ID  = "wb-swift-sprout-7231"
CDR_DATASET = os.environ.get("WORKSPACE_CDR", "")   # e.g. "project.C2024Q2R9" — set manually if blank
REGION      = "us-central1"
SERVICE_ACCOUNT = "pet-27799165194323faf22e2@wb-swift-sprout-7231.iam.gserviceaccount.com"
NETWORK    = f"projects/{PROJECT_ID}/global/networks/network"
SUBNETWORK = f"projects/{PROJECT_ID}/regions/{REGION}/subnetworks/subnetwork"

WS    = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
WS_GS = "gs://cloned-shared-env-pilot-wb-swift-sprout-7231"
B, B_GS = f"{WS}/phenotypic_covariance_v9", f"{WS_GS}/phenotypic_covariance_v9"

ROUND2_KEEP    = f"{B}/1kg_eur/01_ancestry/round2/1kg_eur_keep_ids.txt"
ROUND2_KEEP_GS = f"{B_GS}/1kg_eur/01_ancestry/round2/1kg_eur_keep_ids.txt"
PANEL_GS       = f"{B_GS}/01_ancestry_filtering/unified_panel/unified_panel_{CDR_VERSION}"
PLINK2_GS      = f"{B_GS}/01_ancestry_filtering/unified_panel/bin/plink2"
KG_ACOUNT_GS   = f"{WS_GS}/1000g_reference/1kg_all_qc.acount"

KG_DIR   = f"{WS}/1000g_reference"
KG_PANEL = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

OUT    = f"{B}/{SAMPLE_SET}/01_ancestry/rare_variant_pca"
OUT_GS = f"{B_GS}/{SAMPLE_SET}/01_ancestry/rare_variant_pca"
LOGS_GS = f"{B_GS}/{SAMPLE_SET}/dsub_logs"
LOCAL  = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_rare")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOCAL, exist_ok=True)

N_PCS_FIT = 10
QC_MACHINE,  QC_MEM_MB,  QC_DISK  = "n1-highmem-32", 200_000, 500

assert os.path.isfile(ROUND2_KEEP), f"{ROUND2_KEEP} -- run round-2 gate notebook first"
print(f"round 2: {sum(1 for _ in open(ROUND2_KEEP)):,} participants")
print(OUT)
```

## Cell 2 — install dsub and plink2

```bash
%%bash
pip install --quiet --upgrade 'dsub>=0.5.3'
DSUB_DIR=$(python -c "import dsub, os; print(os.path.dirname(dsub.__file__))")
sed -i -E "s|cloud-sdk:[0-9]+\.[0-9]+\.[0-9]+-slim|cloud-sdk:581.0.0-slim|g" \
  "${DSUB_DIR}/providers/google_utils.py"
dsub --version

BIN="$HOME/bin/plink2"
if [ ! -x "$BIN" ]; then
  mkdir -p "$HOME/bin"
  wget -q -O /tmp/plink2.zip "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
  unzip -o -q /tmp/plink2.zip plink2 -d "$HOME/bin" && chmod +x "$BIN"
fi
export PATH="$HOME/bin:$PATH"
plink2 --version
```

## Cell 3 — QC and HM3 split (Batch)

Filters unified panel to r2 participants and the 0.1%–1% MAF window. No HWE
filter. Splits into HM3-agreeing and non-HM3 by ID+REF+ALT match against the
1000G acount.

```python
os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"

HIGH_LD_REGIONS_GRCH38 = """\
chr1 46761740 52761740 1
chr2 84919365 101517106 2
chr2 182427027 190427029 10
chr3 46483505 50987563 11
chr3 82368158 87868160 12
chr5 43464140 52168409 13
chr5 128636407 133636409 14
chr6 24391792 34424245 15
chr6 56788603 59453888 17
chr6 60109122 62357029 18
chr6 138637169 143137170 21
chr7 53964812 67897578 22
chr8 7105067 13105082 34
chr8 42025699 49924888 35
chr8 109918594 114918595 37
chr10 35671065 44184546 38
chr11 87127183 92127184 41
chr12 31955798 42319931 42
chr17 42159541 44159574 26
chr20 32948532 37438183 43
"""
LD_REGIONS_GS = f"{OUT_GS}/high_ld_regions_grch38.txt"
open(f"{LOCAL}/high_ld_regions_grch38.txt", "w").write(HIGH_LD_REGIONS_GRCH38)
subprocess.run(["bash", "-c",
    f"gsutil cp {LOCAL}/high_ld_regions_grch38.txt {LD_REGIONS_GS}"], check=True)

qc_cmd = """
set -eo pipefail
chmod +x "$PLINK_BIN"
Q="${TMPDIR:-/tmp}/rare_qc"

"$PLINK_BIN" --pgen "$PANEL_PGEN" --pvar "$PANEL_PVAR" --psam "$PANEL_PSAM" \
  --keep "$KEEP_PATH" --nonfounders \
  --maf 0.001 --max-maf 0.00999 --geno 0.05 --max-alleles 2 --rm-dup exclude-all \
  --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "$Q"
echo "rare_qc: $(grep -vc '^##' ${Q}.pvar) variants, $(wc -l < ${Q}.psam) samples"

"$PLINK_BIN" --pfile "$Q" --nonfounders --freq counts \
  --threads "$VCPUS" --memory "$MEM_MB" --out "$Q"

# HM3 split: ID+REF+ALT match against 1000G acount
grep -v '^##' "${Q}.pvar" | awk 'NR>1 {print $3, $4, $5}' | LC_ALL=C sort > "${Q}_ira.sorted"
awk 'NR>1 {print $2, $3, $4}' "$KG_ACOUNT" | LC_ALL=C sort > "${TMPDIR:-/tmp}/kg_ira.sorted"
LC_ALL=C comm -12 "${Q}_ira.sorted" "${TMPDIR:-/tmp}/kg_ira.sorted" \
  | awk '{print $1}' > "${OUT_DIR}/hm3_rare_agreeing.ids"
echo "HM3-agreeing rare: $(wc -l < ${OUT_DIR}/hm3_rare_agreeing.ids)"

cp "${Q}.pgen" "${Q}.pvar" "${Q}.psam" "${Q}.acount" "${OUT_DIR}/"
"""
open(f"{LOCAL}/rare_qc_cmd.sh", "w").write(qc_cmd)
subprocess.run(["bash", "-c",
    f"gsutil cp {LOCAL}/rare_qc_cmd.sh {OUT_GS}/rare_qc_cmd.sh"], check=True)
subprocess.run(["bash", "-c",
    f"gsutil cp {ROUND2_KEEP} {ROUND2_KEEP_GS}"], check=False)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "rare-qc-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {QC_MACHINE} --disk-size {QC_DISK} \
  --input PANEL_PGEN="{PANEL_GS}.pgen" \
  --input PANEL_PVAR="{PANEL_GS}.pvar" \
  --input PANEL_PSAM="{PANEL_GS}.psam" \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input KEEP_PATH="{ROUND2_KEEP_GS}" \
  --input KG_ACOUNT="{KG_ACOUNT_GS}" \
  --env VCPUS={QC_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={QC_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/qc" \
  --script {OUT_GS}/rare_qc_cmd.sh
"""], capture_output=True, text=True)
print(job.stdout or job.stderr)
QC_JOB = job.stdout.strip().splitlines()[-1] if job.stdout else None
```

## Cell 4 — watch

```python
print(subprocess.run(["bash", "-c",
    f"dstat --provider google-batch --project {PROJECT_ID} --location {REGION} "
    f"--jobs {QC_JOB} --users '*' --status '*' --full"],
    capture_output=True, text=True).stdout[-3000:])
```

## Cell 5 — pull results and count

```python
subprocess.run(["bash", "-c", f"gsutil -m cp -r {OUT_GS}/qc/ {LOCAL}/"], check=True)
```

```python
acount = pd.read_csv(f"{LOCAL}/qc/rare_qc.acount", sep=r"\s+")
acount = acount.rename(columns={"#CHROM": "CHROM"})
acount["AF"]  = acount["ALT_CTS"] / acount["OBS_CT"]
acount["MAF"] = np.minimum(acount["AF"], 1 - acount["AF"])

hm3_ids = set(open(f"{LOCAL}/qc/hm3_rare_agreeing.ids").read().split())
acount["hm3"] = acount["ID"].isin(hm3_ids)

hm3    = acount[acount["hm3"]]
nonhm3 = acount[~acount["hm3"]]
print(f"HM3-agreeing rare:  {len(hm3):>10,}")
print(f"non-HM3 rare:       {len(nonhm3):>10,}")
print(f"total rare variants:{len(acount):>10,}")

bins = np.linspace(0.001, 0.01, 50)
fig, axes = plt.subplots(1, 2, figsize=(13, 4))
for ax, yscale in zip(axes, ("linear", "log")):
    ax.hist(hm3["MAF"],    bins=bins, alpha=0.6, color="tab:blue",
            label=f"HM3-agreeing ({len(hm3):,})")
    ax.hist(nonhm3["MAF"], bins=bins, alpha=0.6, color="tab:orange",
            label=f"non-HM3 ({len(nonhm3):,})")
    ax.set_xlabel("MAF"); ax.set_ylabel("variants")
    ax.set_yscale(yscale)
    ax.legend(fontsize=8)
fig.suptitle(f"{SAMPLE_SET} — rare variant MAF distribution (0.1%–1%, r2 cohort)")
plt.tight_layout()
plt.savefig(f"{OUT}/rare_maf_distributions.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 6 — write ID lists and prune (local)

```python
for label, df in (("hm3_rare_all", hm3), ("nonhm3_rare_all", nonhm3)):
    df["ID"].to_csv(f"{LOCAL}/{label}.ids", index=False, header=False)

LD_REGIONS = f"{LOCAL}/high_ld_regions_grch38.txt"
open(LD_REGIONS, "w").write(HIGH_LD_REGIONS_GRCH38)

groups = ["hm3_rare_all", "nonhm3_rare_all"]
for label in groups:
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile {LOCAL}/qc/rare_qc --extract {LOCAL}/{label}.ids \
  --exclude bed1 {LD_REGIONS} --nonfounders \
  --indep-pairwise 1000kb 1 0.1 \
  --threads $(nproc) --out {LOCAL}/{label}
echo "{label}: $(wc -l < {LOCAL}/{label}.prune.in) variants after pruning"
"""], check=True)
```

## Cell 7 — PCA per group (local)

```python
import tempfile

os.makedirs(f"{LOCAL}/pca", exist_ok=True)
tasks = [g for g in groups if os.path.isfile(f"{LOCAL}/{g}.prune.in")]

for arm in tasks:
    out = f"{LOCAL}/pca/{arm}"
    with tempfile.TemporaryDirectory() as tmp:
        p = f"{tmp}/{arm}"
        subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile {LOCAL}/qc/rare_qc --extract {LOCAL}/{arm}.prune.in --nonfounders \
  --threads $(nproc) --make-pgen --out {p}
echo "{arm}: $(grep -vc '^##' {p}.pvar) variants going into PCA"
plink2 --pfile {p} --nonfounders --freq counts \
  --pca approx {N_PCS_FIT} allele-wts \
  --threads $(nproc) --out {out}
W={out}.eigenvec.allele
HEADER=$(head -1 "$W")
ID=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'ID'  | cut -d: -f1)
A1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'A1'  | cut -d: -f1)
P1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC1' | cut -d: -f1)
PK=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)
plink2 --pfile {p} --nonfounders \
  --read-freq {out}.acount \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${{P1}}-${{PK}}" \
  --threads $(nproc) --out {out}_scores
echo "{arm}: done"
"""], check=True)
```

## Cell 8 — read scores

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]

def read_scores(path, id_name):
    d = pd.read_csv(path, sep=r"\s+")
    d = d.rename(columns={("#IID" if "#IID" in d.columns else "IID"): id_name})
    d = d.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS_FIT + 1)})
    d[id_name] = d[id_name].astype(str)
    return d[[id_name] + PC]

scores = {}
for label in groups:
    path = f"{LOCAL}/pca/{label}_scores.sscore"
    if os.path.isfile(path):
        scores[label] = read_scores(path, "person_id")
        print(f"{label}: {len(scores[label]):,} participants")
```

## Cell 9 — PC correlation: HM3 vs non-HM3 rare

```python
K_CMP = min(5, N_PCS_FIT)
l_hm3, l_non = "hm3_rare_all", "nonhm3_rare_all"
if l_hm3 in scores and l_non in scores:
    j = scores[l_hm3].merge(scores[l_non], on="person_id", suffixes=("_hm3", "_non"))
    r = [np.corrcoef(j[f"{p}_hm3"], j[f"{p}_non"])[0, 1] for p in PC[:K_CMP]]
    print("HM3 vs non-HM3 rare |r| by PC:")
    print(pd.Series(np.abs(r), index=PC[:K_CMP]).round(3).to_string())
```

## Cell 10 — load reference coordinates

Three reference systems to colour the rare-variant PCA:
1. Common-variant PC1 from the 01b HM3 common run (same r2 participants)
2. AoU ancestry PCs (from the platform's own 1000G+HGDP training)
3. Mahalanobis distance from the round-2 CEU+GBR gate

```python
from scipy.linalg import inv as matinv

# --- 1. common-variant PCs (01b HM3 common run) ---
LOCAL_01B = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_round2")
R2_COMMON = f"{LOCAL_01B}/pca/hm3_common_scores.sscore"
N_R2_PCS  = 10   # PCs fitted in 01b

def read_scores_generic(path, id_col):
    d = pd.read_csv(path, sep=r"\s+")
    d = d.rename(columns={("#IID" if "#IID" in d.columns else "IID"): id_col})
    d = d.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_R2_PCS + 1)})
    d[id_col] = d[id_col].astype(str)
    return d

common_scores = None
if os.path.isfile(R2_COMMON):
    common_scores = read_scores_generic(R2_COMMON, "person_id")
    print(f"common-variant scores: {len(common_scores):,}")
else:
    print(f"WARNING: {R2_COMMON} not found — skip common-PC colouring")

# --- 2. AoU ancestry PCs ---
import ast

AUX = os.path.expanduser(
    "~/workspace/cdrv9/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux")
ANCESTRY_PREDS = f"{AUX}/ancestry/ancestry_preds.tsv"
_anc_raw = pd.read_csv(ANCESTRY_PREDS, sep="\t",
                        usecols=["research_id", "pca_features"])
_anc_raw = _anc_raw.rename(columns={"research_id": "person_id"})
_anc_raw["person_id"] = _anc_raw["person_id"].astype(str)
_X = np.vstack(_anc_raw["pca_features"].map(ast.literal_eval).to_numpy())
anc_pcs = [f"anc_PC{k+1}" for k in range(_X.shape[1])]
anc = _anc_raw[["person_id"]].copy()
for i, col in enumerate(anc_pcs):
    anc[col] = _X[:, i]
print(f"AoU ancestry PCs: {len(anc_pcs)} PCs for {len(anc):,} participants")

# --- 3. Mahalanobis distance from round-2 CEU+GBR gate ---
# Re-derive from round-2 participant and 1000G scores
KG_PANEL  = f"{os.path.expanduser('~/workspace/Data from All of Us Controlled Tier /shared-env-pilot')}/1000g_reference/integrated_call_samples_v3.20130502.ALL.panel"
KG_SCORES = f"{LOCAL_01B}/kg_in_r2.sscore"       # written by 02_round2 Cell 8
R2_PART_SCORES = f"{LOCAL_01B}/part_in_r2.sscore"  # participant scores from 02

ANCHOR_POPS = ["CEU", "GBR"]
K_PCS_GATE  = 4     # PCs used for the gate in 02_round2 (set K_PCS there)

mahal = None
if os.path.isfile(KG_SCORES) and os.path.isfile(R2_PART_SCORES):
    kg_meta = pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop"]]
    kg = read_scores_generic(KG_SCORES, "sample").merge(kg_meta, on="sample", how="left")
    part_r2 = read_scores_generic(R2_PART_SCORES, "person_id")
    USE = [f"PC{k}" for k in range(1, K_PCS_GATE + 1)]
    anchor = kg[kg["pop"].isin(ANCHOR_POPS)][USE].to_numpy()
    mu = anchor.mean(0)
    C  = np.cov(anchor, rowvar=False)
    Cinv = matinv(C)
    def mahal_dist(df):
        d = df[USE].to_numpy() - mu
        return np.sqrt(np.einsum("ij,jk,ik->i", d, Cinv, d))
    mahal = part_r2[["person_id"]].copy()
    mahal["d_mahal"] = mahal_dist(part_r2)
    print(f"Mahalanobis distances: {len(mahal):,} participants")
else:
    print(f"WARNING: round-2 score files not found — skip Mahalanobis colouring")
    print(f"  expected: {KG_SCORES}")
    print(f"  expected: {R2_PART_SCORES}")
```

## Cell 11 — assign EUR subpop proximity labels

Uses the common-variant PC space to label each r2 participant by their nearest
1000G EUR subpopulation centroid. These labels are then used to show "where the
reference populations sit" in the rare-variant PCA space.

```python
EUR_POPS    = ["CEU", "GBR", "FIN", "TSI", "IBS"]
POP_COLORS  = {"CEU": "#e6194b", "GBR": "#f58231", "FIN": "#3cb44b",
               "TSI": "#4363d8", "IBS": "#911eb4"}
N_COMMON_PCS = 3   # PCs to use for proximity assignment

pop_labels = None
if common_scores is not None and mahal is not None:
    kg_meta = pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop"]]
    # read_scores_generic uses N_R2_PCS; rename to avoid collision with rare PCs
    kg_common = read_scores_generic(
        f"{LOCAL_01B}/kg_in_r2.sscore", "sample").merge(kg_meta, on="sample", how="left")
    kg_common = kg_common.rename(
        columns={f"PC{k}": f"cPC{k}" for k in range(1, N_R2_PCS + 1)})

    USE_C = [f"cPC{k}" for k in range(1, N_COMMON_PCS + 1)]
    centroids = {p: kg_common.loc[kg_common["pop"] == p, USE_C].mean().to_numpy()
                 for p in EUR_POPS}

    cs_renamed = common_scores.rename(
        columns={f"PC{k}": f"cPC{k}" for k in range(1, N_R2_PCS + 1)})
    X = cs_renamed[USE_C].to_numpy()
    dists = np.stack([np.linalg.norm(X - c, axis=1) for c in centroids.values()], axis=1)
    nearest = np.array(EUR_POPS)[dists.argmin(axis=1)]
    pop_labels = cs_renamed[["person_id"]].copy()
    pop_labels["nearest_pop"] = nearest
    for pop in EUR_POPS:
        n = (nearest == pop).sum()
        print(f"  {pop}: {n:,} participants nearest")
```

## Cell 12 — batch and center metadata

Loads participant-level batch axes for detecting technical structure in the PCA:

- **data_partner_id**: EHR collection / recruitment site (from CDR `person` table)
- **sequencing batch**: from the WGS aux manifest if present

Categorical, so plotted as a separate thinned scatter figure (Cell 13) rather
than as hexbin median.

```python
from google.cloud import bigquery

batch_meta = {}   # col_name -> DataFrame[person_id, col_name]

# CDR_DATASET set in Cell 1; fill it in manually if the env var was blank:
#   CDR_DATASET = "aou-res-curation-output-prod.C2024Q2R9"

# --- EHR / collection site ---
if CDR_DATASET:
    client = bigquery.Client()
    dp = client.query(f"""
        SELECT person_id, data_partner_id
        FROM `{CDR_DATASET}.person`
    """).to_dataframe()
    dp["person_id"] = dp["person_id"].astype(str)
    # keep top-40 sites by count; collapse the rest to "other"
    top_sites = dp["data_partner_id"].value_counts().head(40).index
    dp["site"] = (dp["data_partner_id"]
                  .where(dp["data_partner_id"].isin(top_sites), other="other")
                  .astype(str))
    batch_meta["site"] = dp[["person_id", "site"]].copy()
    print(f"Collection sites loaded: {dp['site'].nunique()} categories")
else:
    print("WORKSPACE_CDR not set — skipping site metadata")

# --- sequencing batch manifest (present in some CDR versions) ---
WGS_AUX = os.path.join(WS, "cdrv9/vwb-aou-datasets-controlled-v9/v9"
                        "/wgs/short_read/snpindel/aux")
for fname in ["sequencing_batch_manifest.tsv", "sample_batch_info.tsv",
              "wgs_sample_list.tsv"]:
    fp = os.path.join(WGS_AUX, fname)
    if os.path.exists(fp):
        bm = pd.read_csv(fp, sep="\t")
        bm.columns = bm.columns.str.lower()
        bm["person_id"] = bm["person_id"].astype(str)
        bcols = [c for c in bm.columns if "batch" in c or "center" in c or "site" in c]
        for bc in bcols[:2]:   # at most 2 extra axes per manifest
            batch_meta[bc] = bm[["person_id", bc]].copy()
            print(f"  {fname}: {bc} — {bm[bc].nunique()} categories")
        break

print("Batch axes available:", list(batch_meta.keys()))
```

## Cell 13 — coloured plots (continuous reference axes)

One figure per arm per PC pair. Columns = reference variables (common PC1-3, AoU
ancestry PC1-3, Mahalanobis distance); overlaid scatter shows participants
coloured by nearest EUR subpop.

```python
PC_PAIRS = [(PC[0], PC[1]), (PC[2], PC[3])]

def hexbin_c(ax, x, y, c, cmap, label):
    hb = ax.hexbin(x, y, C=c, reduce_C_function=np.median,
                   gridsize=70, cmap=cmap, mincnt=1, linewidths=0)
    plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04, label=label)

REF_VARS = []
if common_scores is not None:
    cs = common_scores.rename(columns={f"PC{k}": f"cPC{k}" for k in range(1, N_R2_PCS+1)})
    for k in range(1, N_COMMON_PCS + 1):
        REF_VARS.append((cs[["person_id", f"cPC{k}"]], f"cPC{k}",
                         f"common PC{k}", "RdYlBu_r"))
for k in range(1, N_COMMON_PCS + 1):
    REF_VARS.append((anc[["person_id", f"anc_PC{k}"]], f"anc_PC{k}",
                     f"AoU ancestry PC{k}", "PiYG"))
if mahal is not None:
    REF_VARS.append((mahal, "d_mahal", "CEU+GBR Mahal. dist.", "YlOrRd"))

rng = np.random.default_rng(42)

for arm in groups:
    if arm not in scores:
        continue
    s = scores[arm]

    for a, b in PC_PAIRS:
        ncols = len(REF_VARS)
        fig, axes = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.5))
        fig.suptitle(f"{arm.replace('_', ' ')} — {a} vs {b}", y=1.01)

        for ax, (ref_df, col, label, cmap) in zip(np.atleast_1d(axes), REF_VARS):
            m = s.merge(ref_df, on="person_id", how="inner")
            hexbin_c(ax, m[a], m[b], m[col], cmap, label)

            # overlay nearest-pop scatter (thinned)
            if pop_labels is not None:
                m2 = s.merge(pop_labels, on="person_id", how="inner")
                for pop, color in POP_COLORS.items():
                    sub = m2[m2["nearest_pop"] == pop]
                    idx = rng.choice(len(sub), size=min(500, len(sub)), replace=False)
                    ax.scatter(sub.iloc[idx][a], sub.iloc[idx][b],
                               s=4, color=color, alpha=0.6, label=pop,
                               linewidths=0, zorder=3, rasterized=True)

            ax.set_xlabel(a); ax.set_ylabel(b); ax.set_title(label, fontsize=9)

        np.atleast_1d(axes)[0].legend(fontsize=7, markerscale=2,
                                       loc="upper left", framealpha=0.6)
        plt.tight_layout()
        plt.savefig(f"{OUT}/rare_coloured_{arm}_{a}_{b}.png",
                    dpi=130, bbox_inches="tight")
        plt.show()
```

## Cell 14 — batch / center scatter plots

One scatter figure per arm per PC pair, one panel per batch axis. Points are
thinned to ≤50,000 to keep rendering fast; colors are drawn from a combined
tab20 × tab20b palette (up to 60 categories). If `batch_meta` is empty this
cell is a no-op.

```python
import itertools, matplotlib.colors as mcolors

def _cat_palette(n):
    """Return up to n distinct colors cycling through tab20/tab20b/tab20c."""
    pools = [plt.cm.tab20.colors, plt.cm.tab20b.colors, plt.cm.tab20c.colors]
    return list(itertools.islice(itertools.chain.from_iterable(pools), n))

MAX_SCATTER = 50_000

for arm in groups:
    if arm not in scores:
        continue
    s = scores[arm]

    for col_name, ref_df in batch_meta.items():
        m = s.merge(ref_df, on="person_id", how="inner")
        if m.empty:
            continue

        cats = m[col_name].astype(str)
        uniq = cats.value_counts().index.tolist()   # ordered by frequency
        palette = _cat_palette(len(uniq))
        cat_color = {c: palette[i] for i, c in enumerate(uniq)}

        # thin to at most MAX_SCATTER points (stratified by category)
        n_each = max(1, MAX_SCATTER // len(uniq))
        rows = []
        for cat in uniq:
            sub = m[cats == cat]
            rows.append(sub.sample(min(len(sub), n_each), random_state=42))
        thin = pd.concat(rows, ignore_index=True).sample(frac=1, random_state=42)

        for a, b in PC_PAIRS:
            fig, ax = plt.subplots(figsize=(6, 5))
            for cat in uniq:
                sub = thin[thin[col_name].astype(str) == cat]
                ax.scatter(sub[a], sub[b], s=3, alpha=0.25,
                           color=cat_color[cat], label=cat,
                           linewidths=0, rasterized=True)
            ax.set_xlabel(a); ax.set_ylabel(b)
            ax.set_title(f"{arm.replace('_', ' ')} — coloured by {col_name}", fontsize=9)
            # legend only if few enough categories to read
            if len(uniq) <= 20:
                ax.legend(fontsize=6, markerscale=3, ncol=2,
                          loc="upper left", framealpha=0.6)
            plt.tight_layout()
            plt.savefig(f"{OUT}/batch_{col_name}_{arm}_{a}_{b}.png",
                        dpi=130, bbox_inches="tight")
            plt.show()
```

## Cell 15 — plain projection plots

```python
POP_COLORS = {"CEU": "#e6194b", "GBR": "#f58231", "FIN": "#3cb44b",
              "TSI": "#4363d8", "IBS": "#911eb4"}
PC_PAIRS = [(PC[0], PC[1]), (PC[2], PC[3])]

for label in groups:
    if label not in scores:
        continue
    s = scores[label]
    fig, axes = plt.subplots(1, len(PC_PAIRS), figsize=(5 * len(PC_PAIRS), 4.5))
    for ax, (a, b) in zip(axes, PC_PAIRS):
        ax.hexbin(s[a], s[b], gridsize=80, cmap="Greys", mincnt=1,
                  bins="log", linewidths=0, zorder=1)
        ax.set_xlabel(a); ax.set_ylabel(b)
    fig.suptitle(f"{label.replace('_', ' ')} — participant cloud (no 1000G projection)")
    plt.tight_layout()
    plt.savefig(f"{OUT}/projection_{label}.png", dpi=130, bbox_inches="tight")
    plt.show()
```

## Cell 16 — loadings by genomic position

```python
for label in groups:
    lpath = f"{LOCAL}/pca/{label}.eigenvec.allele"
    if not os.path.isfile(lpath):
        continue
    L = pd.read_csv(lpath, sep=r"\s+")
    idc = "#ID" if "#ID" in L.columns else "ID"
    L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
    L["CHROM"] = L["CHROM"].str.replace("chr", "", regex=False)
    L = L[L["CHROM"].isin([str(c) for c in range(1, 23)])].copy()
    L["CHROM"] = pd.Categorical(L["CHROM"], categories=[str(c) for c in range(1, 23)], ordered=True)
    L["POS"] = L["POS"].astype(int)
    L = L.sort_values(["CHROM", "POS"])

    offset, centres = 0, {}
    cum = np.empty(len(L))
    for c, idx in L.groupby("CHROM", observed=True).indices.items():
        p = L["POS"].to_numpy()[idx]
        cum[idx] = p + offset
        centres[c] = offset + p.max() / 2
        offset += p.max()
    L["CUM"] = cum

    n_show = min(4, N_PCS_FIT)
    fig, axs = plt.subplots(n_show, 1, figsize=(14, 2.4 * n_show), sharex=True)
    for ax, pc in zip(np.atleast_1d(axs), PC[:n_show]):
        for k, c in enumerate(centres):
            s = L[L["CHROM"] == c]
            ax.scatter(s["CUM"], s[pc] ** 2, s=2, alpha=0.5,
                       color=["tab:blue", "tab:orange"][k % 2], rasterized=True)
        ax.set_ylabel(f"{pc} loading²")
    np.atleast_1d(axs)[-1].set_xticks(list(centres.values()))
    np.atleast_1d(axs)[-1].set_xticklabels(list(centres), fontsize=7)
    np.atleast_1d(axs)[-1].set_xlabel("chromosome")
    fig.suptitle(f"{label.replace('_', ' ')} — PCA loadings")
    plt.tight_layout()
    plt.savefig(f"{OUT}/loadings_{label}.png", dpi=130, bbox_inches="tight")
    plt.show()
```
