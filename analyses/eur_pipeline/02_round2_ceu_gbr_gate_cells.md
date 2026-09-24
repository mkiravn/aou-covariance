# 1kg_eur — round 2: CEU + GBR gate

Refits PCA on the round-1 participants, projects 1000G into that space, and
keeps those closest to the CEU + GBR centroid. Two alternative gates are built
for comparison: the participants' own centroid (as eur_D2 did), and the Kemper
direction, PCA fit on the 1000G Europeans with participants projected in.

The gate is a multivariate normal fitted to the anchor population: mean and
covariance from the projected CEU + GBR samples, Mahalanobis distance under it,
and a radius tuned to the coverage level. The anchor supplies the orientation
and relative PC weighting; projection shrinks its scores toward the origin but
that only affects the scale, which the coverage quantile absorbs.

Compute: QC+prune on n1-highmem-16; PCA on n1-highmem-32. The PCA fit is the
long step. Scoring runs on this VM after the jobs finish.

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

ROUND1_KEEP    = f"{WS}/phenotypic_covariance_v9/analyses/ancestry_filtering/keep/EUR_99pct_keep_ids.txt"
ROUND1_KEEP_GS = f"{WS_GS}/phenotypic_covariance_v9/analyses/ancestry_filtering/keep/EUR_99pct_keep_ids.txt"
PANEL_GS       = f"{B_GS}/01_ancestry_filtering/unified_panel/unified_panel_{CDR_VERSION}"
PLINK2_GS      = f"{B_GS}/01_ancestry_filtering/unified_panel/bin/plink2"
KG_BFILE_GS    = f"{WS_GS}/1000g_reference/1kg_all_qc"
KG_ACOUNT_GS   = f"{WS_GS}/1000g_reference/1kg_all_qc.acount"

KG_DIR   = f"{WS}/1000g_reference"
KG_BFILE = f"{KG_DIR}/1kg_all_qc"
KG_PANEL = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

OUT    = f"{B}/{SAMPLE_SET}/01_ancestry/round2"
OUT_GS = f"{B_GS}/{SAMPLE_SET}/01_ancestry/round2"
LOGS_GS = f"{B_GS}/{SAMPLE_SET}/dsub_logs"
LOCAL  = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_round2")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOCAL, exist_ok=True)

ANCHOR_POPS = ["CEU", "GBR"]
EUR_POPS    = ["CEU", "GBR", "FIN", "TSI", "IBS"]
K_PCS       = 4                 # set after the scree in Cell 6
N_PCS_FIT   = 20
ANCHOR_PCT  = 0.90              # fraction of CEU+GBR founders inside the gate

QC_MACHINE,  QC_MEM_MB,  QC_DISK  = "n1-highmem-16",  96_000, 500
PCA_MACHINE, PCA_MEM_MB, PCA_DISK = "n1-highmem-32", 200_000, 500

assert os.path.isfile(ROUND1_KEEP), f"{ROUND1_KEEP} -- run the ancestry survey first"
print(f"round 1: {sum(1 for _ in open(ROUND1_KEEP)):,} participants")
assert os.path.isfile(f"{KG_BFILE}.bed"), KG_BFILE
assert os.path.isfile(KG_PANEL), KG_PANEL
print(OUT)
```

## Cell 2 — dsub

```bash
%%bash
set -e
pip install --quiet --upgrade 'dsub>=0.5.3'
DSUB_DIR=$(python -c "import dsub, os; print(os.path.dirname(dsub.__file__))")
sed -i -E "s|cloud-sdk:[0-9]+\.[0-9]+\.[0-9]+-slim|cloud-sdk:581.0.0-slim|g" \
  "${DSUB_DIR}/providers/google_utils.py"
dsub --version
grep -o "cloud-sdk:[^'\"]*" "${DSUB_DIR}/providers/google_utils.py" | sort -u
```

## Cell 3 — QC, HM3 filtering, and LD pruning (Batch)

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
LD_REGIONS_GS = f"{OUT_GS}/high_ld_regions_grch38.txt"
open(f"{LOCAL}/high_ld_regions_grch38.txt", "w").write(HIGH_LD_REGIONS_GRCH38)

qc_cmd = """
    set -eo pipefail
    chmod +x "$PLINK_BIN"
    Q="${TMPDIR:-/tmp}/r1_qc"

    "$PLINK_BIN" --pgen "$PANEL_PGEN" --pvar "$PANEL_PVAR" --psam "$PANEL_PSAM" \
      --keep "$KEEP_PATH" --nonfounders \
      --maf 0.01 --hwe 1e-6 0 keep-fewhet --geno 0.05 --max-alleles 2 --rm-dup exclude-all \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "$Q"

    grep -v '^##' "${Q}.pvar" | awk 'NR>1 {print $3, $4, $5}' | LC_ALL=C sort > "${Q}_ira.sorted"
    awk 'NR>1 {print $2, $3, $4}' "$KG_ACOUNT" | LC_ALL=C sort > "${Q}_kg.sorted"
    LC_ALL=C comm -12 "${Q}_ira.sorted" "${Q}_kg.sorted" | awk '{print $1}' > "${OUT_DIR}/hm3_agreeing.ids"
    echo "HM3-agreeing: $(wc -l < "${OUT_DIR}/hm3_agreeing.ids")"

    "$PLINK_BIN" --pfile "$Q" --extract "${OUT_DIR}/hm3_agreeing.ids" \
      --exclude bed1 "$LD_REGIONS" --nonfounders --indep-pairwise 1000kb 1 0.05 \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${TMPDIR:-/tmp}/prune"

    "$PLINK_BIN" --pfile "$Q" --extract "${TMPDIR:-/tmp}/prune.prune.in" \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "${OUT_DIR}/pca_input"

    cp "${TMPDIR:-/tmp}/prune.prune.in" "${OUT_DIR}/prune.prune.in"
    echo "PCA variants: $(wc -l < "${OUT_DIR}/prune.prune.in")"
"""
open(f"{LOCAL}/qc_cmd.sh", "w").write(qc_cmd)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "r2-qc-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {QC_MACHINE} --disk-size {QC_DISK} \
  --input PANEL_PGEN="{PANEL_GS}.pgen" \
  --input PANEL_PVAR="{PANEL_GS}.pvar" \
  --input PANEL_PSAM="{PANEL_GS}.psam" \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input KEEP_PATH="{ROUND1_KEEP_GS}" \
  --input LD_REGIONS="{LD_REGIONS_GS}" \
  --input KG_ACOUNT="{KG_ACOUNT_GS}" \
  --env VCPUS={QC_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={QC_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/panels" \
  --script {LOCAL}/qc_cmd.sh
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

## Cell 5 — PCA (Batch)

```python
pca_cmd = f"""
    set -eo pipefail
    chmod +x "$PLINK_BIN"
    P="${{TMPDIR:-/tmp}}/pca_input"

    "$PLINK_BIN" --pgen "$QC_PGEN" --pvar "$QC_PVAR" --psam "$QC_PSAM" \
      --extract "$PRUNE_IN" --nonfounders \
      --threads "$VCPUS" --memory "$MEM_MB" --make-pgen --out "$P"

    "$PLINK_BIN" --pfile "$P" --nonfounders --freq counts \
      --pca approx {N_PCS_FIT} allele-wts \
      --threads "$VCPUS" --memory "$MEM_MB" --out "${{OUT_DIR}}/round2_pca"

    echo "done"
"""
open(f"{LOCAL}/pca_cmd.sh", "w").write(pca_cmd)

job = subprocess.run(["bash", "-c", f"""
set -e
dsub --provider google-batch --project {PROJECT_ID} --regions {REGION} \
  --logging {LOGS_GS} --service-account {SERVICE_ACCOUNT} \
  --network {NETWORK} --subnetwork {SUBNETWORK} --use-private-address \
  --image "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim" \
  --name "r2-pca-{SAMPLE_SET.replace('_', '-')}" \
  --machine-type {PCA_MACHINE} --disk-size {PCA_DISK} \
  --input PLINK_BIN="{PLINK2_GS}" \
  --input QC_PGEN="{OUT_GS}/panels/pca_input.pgen" \
  --input QC_PVAR="{OUT_GS}/panels/pca_input.pvar" \
  --input QC_PSAM="{OUT_GS}/panels/pca_input.psam" \
  --input PRUNE_IN="{OUT_GS}/panels/prune.prune.in" \
  --env VCPUS={PCA_MACHINE.rsplit('-', 1)[-1]} --env MEM_MB={PCA_MEM_MB} \
  --output-recursive OUT_DIR="{OUT_GS}/pca" \
  --script {LOCAL}/pca_cmd.sh
"""], capture_output=True, text=True)
print(job.stdout or job.stderr)
PCA_JOB = job.stdout.strip().splitlines()[-1] if job.stdout else None
```

## Cell 6 — watch

```python
print(subprocess.run(["bash", "-c",
    f"dstat --provider google-batch --project {PROJECT_ID} --location {REGION} "
    f"--jobs {PCA_JOB} --users '*' --status '*' --full"],
    capture_output=True, text=True).stdout[-3000:])
```

## Cell 7 — plink2 on this VM (for scoring)

```bash
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
"$HOME/bin/plink2" --version
```

```python
os.environ.update(dict(
    LOCAL=LOCAL, OUT_GS=OUT_GS,
    KG_BFILE=KG_BFILE, KG_PANEL=KG_PANEL,
    N_PCS_FIT=str(N_PCS_FIT),
))
os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"
```

## Cell 8 — score everyone through the same loadings, pick K

plink's `--pca` eigenvectors and a `--score` projection do not land on the
same coordinates, so participants are re-scored through their own loadings
rather than read from the eigenvec. Both sets then sit on identical axes.

```bash
%%bash
set -eo pipefail
THREADS=$(nproc)
W="$LOCAL/round2_pca.eigenvec.allele"
FREQ="$LOCAL/round2_pca.acount"

# pull results from the bucket
gsutil -m cp \
  "$OUT_GS/pca/round2_pca.eigenval" \
  "$OUT_GS/pca/round2_pca.eigenvec.allele" \
  "$OUT_GS/pca/round2_pca.acount" \
  "$OUT_GS/panels/pca_input.pvar" \
  "$OUT_GS/panels/prune.prune.in" \
  "$LOCAL/"

HEADER=$(head -1 "$W")
ID=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'ID'  | cut -d: -f1)
A1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'A1'  | cut -d: -f1)
P1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC1' | cut -d: -f1)
PK=$(echo "$HEADER" | tr '\t' '\n' | grep -nx "PC${N_PCS_FIT}" | cut -d: -f1)

# find variants shared with 1000G for projection
grep -v '^##' "$LOCAL/pca_input.pvar" | awk 'NR>1 {print $3, $4, $5}' | LC_ALL=C sort > "$LOCAL/pca_ira.sorted"
awk 'NR>1 {print $2, $3, $4}' "$KG_BFILE.acount" | LC_ALL=C sort > "$LOCAL/kg_ira.sorted"
LC_ALL=C comm -12 "$LOCAL/pca_ira.sorted" "$LOCAL/kg_ira.sorted" | awk '{print $1}' > "$LOCAL/kg_project.ids"
echo "variants for the 1000G projection: $(wc -l < "$LOCAL/kg_project.ids")"

# score 1000G into the participant PCA space
plink2 --bfile "$KG_BFILE" --extract "$LOCAL/kg_project.ids" \
  --read-freq "$FREQ" \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${P1}-${PK}" --threads "$THREADS" --out "$LOCAL/kg_in_r2"

# re-score participants through their own loadings
plink2 --pgen "$LOCAL/pca_input.pgen" --pvar "$LOCAL/pca_input.pvar" --psam "$LOCAL/pca_input.psam" \
  --read-freq "$FREQ" \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${P1}-${PK}" --threads "$THREADS" --out "$LOCAL/part_in_r2"
```

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]


def read_scores(path, id_name):
    d = pd.read_csv(path, sep=r"\s+")
    d = d.rename(columns={("#IID" if "#IID" in d.columns else "IID"): id_name})
    d = d.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS_FIT + 1)})
    d[id_name] = d[id_name].astype(str)
    return d[[id_name] + PC]


kg = read_scores(f"{LOCAL}/kg_in_r2.sscore", "sample").merge(
    pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]], on="sample", how="left")
part = read_scores(f"{LOCAL}/part_in_r2.sscore", "person_id")

ev = np.loadtxt(f"{LOCAL}/round2_pca.eigenval")
sep = [np.abs(kg[kg["pop"].isin(ANCHOR_POPS)][p].mean()
              - kg[kg["pop"].isin(set(EUR_POPS) - set(ANCHOR_POPS))][p].mean()) / part[p].std()
       for p in PC]
print(pd.DataFrame({"PC": PC, "pct_var": (ev / ev.sum() * 100).round(2),
                    "anchor_vs_other_EUR": np.round(sep, 2)}).to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(range(1, len(ev) + 1), ev / ev.sum() * 100, color="royalblue")
ax.set_xlabel("PC"); ax.set_ylabel("% variance")
plt.tight_layout(); plt.savefig(f"{OUT}/scree.png", dpi=130, bbox_inches="tight"); plt.show()
```

`anchor_vs_other_EUR` is how far the anchor sits from the other European
populations on each PC, in participant SDs. Use the PCs where it is large.

## Cell 8b — loadings

Loading² by genomic position. A PC carried by one region is an inversion, an
LD block or a mapping artefact, not structure.

```python
L = pd.read_csv(f"{LOCAL}/round2_pca.eigenvec.allele", sep=r"\s+")
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

n_show = min(K_PCS, 4)
fig, axes = plt.subplots(n_show, 1, figsize=(14, 2.4 * n_show), sharex=True)
for ax, p in zip(np.atleast_1d(axes), PC[:n_show]):
    for k, c in enumerate(centres):
        s = L[L["CHROM"] == c]
        ax.scatter(s["CUM"], s[p] ** 2, s=2, alpha=0.5,
                   color=["tab:blue", "tab:orange"][k % 2], rasterized=True)
    ax.set_ylabel(f"{p} loading²")
np.atleast_1d(axes)[-1].set_xticks(list(centres.values()))
np.atleast_1d(axes)[-1].set_xticklabels(list(centres), fontsize=7)
np.atleast_1d(axes)[-1].set_xlabel("chromosome")
fig.suptitle(f"{SAMPLE_SET} round 2 — PCA loadings")
plt.tight_layout()
plt.savefig(f"{OUT}/loadings.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 9 — the gate

```python
USE = PC[:K_PCS]


def mvn(frame):
    A = frame[USE].to_numpy()
    assert len(A) > 2 * K_PCS, f"too few samples for a {K_PCS}x{K_PCS} covariance"
    return A.mean(0), np.cov(A, rowvar=False)


def distance(frame, mu, C):
    d = frame[USE].to_numpy() - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, np.linalg.inv(C), d))


anchor = kg[kg["pop"].isin(ANCHOR_POPS)]
mu, C = mvn(anchor)
print(f"anchor: {len(anchor)} samples; SD per PC {np.sqrt(np.diag(C)).round(4)}")

d_anchor_self = distance(anchor, mu, C)
THRESHOLD = float(np.quantile(d_anchor_self, ANCHOR_PCT))
print(f"threshold: {ANCHOR_PCT:.0%} of anchor distances = {THRESHOLD:.4f}")

d_part = distance(part, mu, C)
keep = d_part <= THRESHOLD
KEEP_PATH = f"{OUT}/{SAMPLE_SET}_keep_ids.txt"
part.loc[keep, "person_id"].to_csv(KEEP_PATH, index=False, header=False)

with open(f"{OUT}/round2_provenance.txt", "w") as f:
    f.write(f"pca\tfit on the round-1 set, plink2 --pca approx {N_PCS_FIT} allele-wts\n"
            f"variants\t{sum(1 for _ in open(f'{LOCAL}/prune.prune.in'))} pruned HM3 sites\n"
            f"anchor\t1000G {'+'.join(ANCHOR_POPS)} projected into that space\n"
            f"pcs\t1-{K_PCS}\n"
            f"metric\tMahalanobis under the anchor's own covariance\n"
            f"threshold\t{ANCHOR_PCT:.0%} quantile of anchor distances = {THRESHOLD:.6f}\n"
            f"kept\t{int(keep.sum())}\n")
print(f"{int(keep.sum()):,} kept -> {KEEP_PATH}")
```

## Cell 10 — plots

```python
from matplotlib.patches import Ellipse

PAIRS = [(0, 1), (2, 3)]
rng = np.random.default_rng(0)
shown = rng.choice(np.flatnonzero(keep), size=min(50_000, int(keep.sum())), replace=False)
cols = dict(zip(EUR_POPS, plt.cm.tab10.colors))


def gate_ellipse(ax, centre, i, j, radius, cov=None, **kw):
    M = (C if cov is None else cov)[np.ix_([i, j], [i, j])]
    vals, vecs = np.linalg.eigh(M)
    ax.add_patch(Ellipse((centre[i], centre[j]), 2 * radius * np.sqrt(vals[-1]),
                         2 * radius * np.sqrt(vals[0]),
                         angle=np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1])),
                         fill=False, **kw))


fig, axes = plt.subplots(1, len(PAIRS), figsize=(15, 6.5))
for ax, (i, j) in zip(axes, PAIRS):
    a, b = PC[i], PC[j]
    ax.scatter(part[a], part[b], s=1, alpha=0.1, color="0.8", rasterized=True, zorder=0)
    ax.scatter(part[a].to_numpy()[shown], part[b].to_numpy()[shown], s=2, alpha=0.3,
               color="tab:blue", rasterized=True, zorder=1)
    for pop in EUR_POPS:
        s = kg[kg["pop"] == pop]
        ax.scatter(s[a], s[b], s=30 if pop in ANCHOR_POPS else 16, marker="x",
                   color=cols[pop], label=pop, zorder=2)
    if i < K_PCS and j < K_PCS:
        gate_ellipse(ax, mu, i, j, THRESHOLD, edgecolor="tab:blue", lw=1.4, zorder=3)
    ax.set_xlabel(a); ax.set_ylabel(b)
axes[0].legend(fontsize=8, markerscale=1.4)
fig.suptitle(f"{SAMPLE_SET} round 2 — kept (blue) among the round-1 set (grey), 1000G Europeans (crosses)")
plt.tight_layout()
plt.savefig(f"{OUT}/round2_pcs.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 11 — alternative gates

Same PC space and the same coverage level, centred on the participants instead
of the anchor; and the Kemper direction, PCA fit on the 1000G Europeans with
participants projected in.

```python
# a. participants' own centroid and covariance, same coverage level
mu_self, C_self = mvn(part)
d_self = distance(part, mu_self, C_self)
keep_self = d_self <= np.quantile(d_self, ANCHOR_PCT)

# b. Kemper direction: fit on the 1000G Europeans, project participants
eur_ids = pd.read_csv(KG_PANEL, sep=r"\s+")
eur_ids = eur_ids.loc[eur_ids["super_pop"].eq("EUR"), "sample"]
eur_ids.to_csv(f"{LOCAL}/kg_eur.keep", index=False, header=False)
```

```bash
%%bash
set -eo pipefail
THREADS=$(nproc)

plink2 --bfile "$KG_BFILE" --keep "$LOCAL/kg_eur.keep" \
  --extract "$LOCAL/kg_project.ids" --nonfounders --freq counts \
  --pca "$N_PCS_FIT" allele-wts --threads "$THREADS" --out "$LOCAL/kgeur_pca"

W="$LOCAL/kgeur_pca.eigenvec.allele"
HEADER=$(head -1 "$W")
ID=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'ID'  | cut -d: -f1)
A1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'A1'  | cut -d: -f1)
P1=$(echo "$HEADER" | tr '\t' '\n' | grep -nx 'PC1' | cut -d: -f1)
PK=$(echo "$HEADER" | tr '\t' '\n' | grep -nx "PC${N_PCS_FIT}" | cut -d: -f1)

plink2 --pgen "$LOCAL/pca_input.pgen" --pvar "$LOCAL/pca_input.pvar" --psam "$LOCAL/pca_input.psam" \
  --extract "$LOCAL/kg_project.ids" \
  --read-freq "$LOCAL/kgeur_pca.acount" \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${P1}-${PK}" --threads "$THREADS" --out "$LOCAL/part_in_kgeur"
```

```python
part_k = read_scores(f"{LOCAL}/part_in_kgeur.sscore", "person_id")
kg_k = pd.read_csv(f"{LOCAL}/kgeur_pca.eigenvec", sep=r"\s+")
kg_k = kg_k.rename(columns={("#IID" if "#IID" in kg_k.columns else "IID"): "sample"})
kg_k = kg_k.merge(pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop"]], on="sample", how="left")

kg_k_anchor = kg_k[kg_k["pop"].isin(ANCHOR_POPS)]
mu_k, C_k = mvn(kg_k_anchor)
d_k = distance(part_k, mu_k, C_k)
thresh_k = np.quantile(distance(kg_k_anchor, mu_k, C_k), ANCHOR_PCT)
keep_kemper = pd.Series(d_k <= thresh_k, index=part_k["person_id"]).reindex(
    part["person_id"]).fillna(False).to_numpy()

sets = {"anchor (this round)": set(part.loc[keep, "person_id"]),
        "participant centroid": set(part.loc[keep_self, "person_id"]),
        "Kemper direction": set(part.loc[keep_kemper, "person_id"])}
d2 = f"{B}/01_ancestry_filtering/r2_eur_gbr_ceu_projection/aou_D2_keep_ids.txt"
if os.path.isfile(d2):
    s = pd.read_csv(d2, sep=r"\s+", header=None, dtype=str).iloc[:, -1]
    sets["eur_D2"] = set(s[s.str.isdigit()])

names = list(sets)
J = pd.DataFrame([[len(sets[a] & sets[b]) / len(sets[a] | sets[b]) for b in names] for a in names],
                 index=names, columns=names)
print("Jaccard overlap")
print(J.round(3).to_string())
J.to_csv(f"{OUT}/gate_overlap.tsv", sep="\t")

only_anchor = keep & ~keep_self
only_self = keep_self & ~keep
fig, axes = plt.subplots(1, len(PAIRS), figsize=(15, 6.5))
for ax, (i, j) in zip(axes, PAIRS):
    a, b = PC[i], PC[j]
    ax.scatter(part[a], part[b], s=1, alpha=0.1, color="0.85", rasterized=True)
    for m, c, lab in ((only_self, "tab:orange", "participant centroid only"),
                      (only_anchor, "tab:blue", "anchor only")):
        ax.scatter(part[a][m], part[b][m], s=2, alpha=0.4, color=c, label=lab, rasterized=True)
    for pop in ANCHOR_POPS:
        s = kg[kg["pop"] == pop]
        ax.scatter(s[a], s[b], s=30, marker="x", color="k", zorder=3)
    if i < K_PCS and j < K_PCS:
        gate_ellipse(ax, mu, i, j, THRESHOLD, edgecolor="tab:blue", lw=1.4, zorder=4)
        gate_ellipse(ax, mu_self, i, j, np.quantile(d_self, ANCHOR_PCT), cov=C_self,
                     edgecolor="tab:orange", lw=1.4, ls="--", zorder=4)
    ax.set_xlabel(a); ax.set_ylabel(b)
axes[0].legend(fontsize=8, markerscale=4)
fig.suptitle("who the anchor changes, against the participants' own centroid")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_difference.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 12 — copy this notebook to the bucket

Save the notebook first, then run this.

```python
NOTEBOOK = os.path.expanduser(f"~/notebooks/{SAMPLE_SET}_round2.ipynb")
NB_DIR = f"{B}/{SAMPLE_SET}/notebooks"
os.makedirs(NB_DIR, exist_ok=True)
assert os.path.isfile(NOTEBOOK), NOTEBOOK
dest = f"{NB_DIR}/{os.path.basename(NOTEBOOK)}"
shutil.copy(NOTEBOOK, dest)
print(f"{'OK  ' if os.path.getsize(dest) == os.path.getsize(NOTEBOOK) else 'SIZE MISMATCH'} {dest}")
```

---

Next: `02a_covariate_pca_cells.md` — covariate PCs on the round-2 keep list.
