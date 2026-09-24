# 1kg_eur — variant set exploration

Compares HM3 and non-HM3 variants from the unified panel at matched frequency
ranges, to inform the choice of variant set for the round-2 ancestry PCA.

Within each frequency range (rare <5%, common >=5%), variants are sampled to
produce a flat MAF histogram so that differences in PCA structure reflect the
variant set rather than its frequency spectrum.

QC outputs are reused from the round-2 local run. Pruning runs locally.
PCA per group runs as Batch tasks (potentially >1 h across all groups).

---

## Cell 1 — config

```python
import os, shutil, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SAMPLE_SET  = "1kg_eur"
CDR_VERSION = "v9"

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
B  = f"{WS}/phenotypic_covariance_v9"

ROUND1_KEEP = f"{B}/analyses/ancestry_filtering/keep/EUR_99pct_keep_ids.txt"
ROUND2_KEEP = f"{B}/1kg_eur/01_ancestry/round2/1kg_eur_keep_ids.txt"

KG_DIR   = f"{WS}/1000g_reference"
KG_BFILE = f"{KG_DIR}/1kg_all_qc"
KG_PANEL = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

OUT   = f"{B}/{SAMPLE_SET}/01_ancestry/variant_exploration"
LOCAL = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_round2")
os.makedirs(OUT, exist_ok=True)

ANCHOR_POPS = ["CEU", "GBR"]
EUR_POPS    = ["CEU", "GBR", "FIN", "TSI", "IBS"]
N_PCS_FIT   = 10

N_RARE   = 50_000
N_COMMON = 100_000
N_BINS   = 20

assert os.path.isfile(ROUND1_KEEP), f"{ROUND1_KEEP} -- run the ancestry survey first"
assert os.path.isfile(ROUND2_KEEP), f"{ROUND2_KEEP} -- run the round-2 gate notebook first"
assert os.path.isfile(f"{KG_BFILE}.bed"), KG_BFILE
print(OUT)
```


## Cell 2 — install plink2

```bash
%%bash
BIN="$HOME/bin/plink2"
if [ ! -x "$BIN" ]; then
  mkdir -p "$HOME/bin"
  wget -q -O /tmp/plink2.zip "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260504.zip"
  unzip -o -q /tmp/plink2.zip plink2 -d "$HOME/bin" && chmod +x "$BIN"
fi
export PATH="$HOME/bin:$PATH"
plink2 --version
```

## Cell 3 — check QC outputs and create r2_qc

Reuses `r1_qc.*` and `hm3_agreeing.ids` from the round-2 gate notebook, then
filters to round-2 participants and reapplies site QC (MAF, HWE, missingness)
in the r2 cohort to produce `r2_qc`.

```python
for f in ("r1_qc.pvar", "r1_qc.acount", "r1_qc.pgen", "r1_qc.psam", "hm3_agreeing.ids"):
    assert os.path.isfile(f"{LOCAL}/{f}"), f"missing: {LOCAL}/{f}"
    print(f"{f}: {os.path.getsize(f'{LOCAL}/{f}'):,} bytes")

os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"

r2_qc = f"{LOCAL}/r2_qc"
if not os.path.isfile(f"{r2_qc}.psam"):
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile {LOCAL}/r1_qc --keep "{ROUND2_KEEP}" --nonfounders \
  --maf 0.01 --hwe 1e-6 0 keep-fewhet --geno 0.05 \
  --threads $(nproc) --make-pgen --out {r2_qc}
echo "r2_qc: $(wc -l < {r2_qc}.psam) samples, $(grep -vc '^##' {r2_qc}.pvar) variants"
"""], check=True)
else:
    n = sum(1 for _ in open(f"{r2_qc}.psam")) - 1
    print(f"r2_qc already exists: {n:,} samples")

if not os.path.isfile(f"{r2_qc}.acount"):
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile {r2_qc} --nonfounders --freq counts \
  --threads $(nproc) --out {r2_qc}
echo "r2_qc.acount: $(wc -l < {r2_qc}.acount) variants"
"""], check=True)
```

## Cell 4 — MAF distributions

```python
pvar = pd.read_csv(f"{LOCAL}/r2_qc.pvar", sep=r"\s+", comment="#",
                   names=["CHROM", "POS", "ID", "REF", "ALT"])
acount = pd.read_csv(f"{LOCAL}/r2_qc.acount", sep=r"\s+")
acount = acount.rename(columns={"#CHROM": "CHROM"})
acount["AF"]  = acount["ALT_CTS"] / acount["OBS_CT"]
acount["MAF"] = np.minimum(acount["AF"], 1 - acount["AF"])

hm3_ids = set(open(f"{LOCAL}/hm3_agreeing.ids").read().split())
acount["hm3"] = acount["ID"].isin(hm3_ids)

hm3    = acount[acount["hm3"]]
nonhm3 = acount[~acount["hm3"]]
print(f"HM3: {len(hm3):,}   non-HM3: {len(nonhm3):,}   total: {len(acount):,}")

bins = np.linspace(0, 0.5, 101)
fig, axes = plt.subplots(1, 2, figsize=(13, 4))
for ax, yscale in zip(axes, ("linear", "log")):
    ax.hist(hm3["MAF"],    bins=bins, alpha=0.6, color="tab:blue",   label=f"HM3 ({len(hm3):,})")
    ax.hist(nonhm3["MAF"], bins=bins, alpha=0.6, color="tab:orange", label=f"non-HM3 ({len(nonhm3):,})")
    ax.set_xlabel("MAF"); ax.set_ylabel("variants")
    ax.set_yscale(yscale)
    ax.legend(fontsize=8)
fig.suptitle(f"{SAMPLE_SET} — HM3 vs non-HM3 MAF distribution (r2 participants)")
plt.tight_layout()
plt.savefig(f"{OUT}/maf_distributions.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5 — flat-MAF sampling

Sample each group to a flat MAF histogram so the PCA comparison is not driven
by frequency spectrum differences. Reports mean/median MAF and plots the
sampled distributions to confirm flatness.

```python
def flat_sample(variants, maf_lo, maf_hi, n_total, n_bins, seed):
    """Sample variants to produce a flat MAF histogram."""
    v = variants[(variants["MAF"] >= maf_lo) & (variants["MAF"] < maf_hi)].copy()
    edges = np.linspace(maf_lo, maf_hi, n_bins + 1)
    v["bin"] = pd.cut(v["MAF"], bins=edges, labels=False, include_lowest=True)
    n_per_bin = min(v["bin"].value_counts().min(), n_total // n_bins)
    rng = np.random.default_rng(seed)
    return (v.groupby("bin", group_keys=False)
             .apply(lambda g: g.sample(min(len(g), n_per_bin),
                                       random_state=rng.integers(1e9))))


groups = {
    "hm3_rare":      flat_sample(hm3,    0.01, 0.05, N_RARE,   N_BINS, seed=0),
    "hm3_common":    flat_sample(hm3,    0.05, 0.50, N_COMMON, N_BINS, seed=1),
    "nonhm3_rare":   flat_sample(nonhm3, 0.01, 0.05, N_RARE,   N_BINS, seed=2),
    "nonhm3_common": flat_sample(nonhm3, 0.05, 0.50, N_COMMON, N_BINS, seed=3),
}

for k, s in groups.items():
    print(f"{k}: {len(s):,}  mean MAF {s['MAF'].mean():.3f}  median {s['MAF'].median():.3f}")

fig, axes = plt.subplots(2, 2, figsize=(12, 7))
for ax, (k, s) in zip(axes.flat, groups.items()):
    lo, hi = (0.01, 0.05) if "rare" in k else (0.05, 0.50)
    ax.hist(s["MAF"], bins=N_BINS, range=(lo, hi), color="tab:blue", edgecolor="none")
    ax.set_title(k); ax.set_xlabel("MAF"); ax.set_ylabel("variants")
plt.tight_layout()
plt.savefig(f"{OUT}/sampled_maf_distributions.png", dpi=130, bbox_inches="tight")
plt.show()

for k, s in groups.items():
    s["ID"].to_csv(f"{LOCAL}/{k}.ids", index=False, header=False)
```

## Cell 6 — prune per (bin, source)

Fast enough to run locally.

```python
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
# 1 Mb buffer added around each canonical high-LD region; small singleton
# intervals from the original list are subsumed by their neighbours.
LD_REGIONS = f"{LOCAL}/high_ld_regions_grch38.txt"
open(LD_REGIONS, "w").write(HIGH_LD_REGIONS_GRCH38)

os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"

for label in groups:
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile {LOCAL}/r2_qc --extract {LOCAL}/{label}.ids \
  --exclude bed1 {LD_REGIONS} --nonfounders \
  --indep-pairwise 1000kb 1 0.05 \
  --threads $(nproc) --out {LOCAL}/{label}
echo "{label}: $(wc -l < {LOCAL}/{label}.prune.in) variants after pruning"
"""], check=True)
```

<!-- r2 run (r²=0.05, 1Mb exclusion buffer, n=223,209):
     hm3_common   97,381 → 22,723  (74,658 pruned, 77%)
     nonhm3_common 95,542 → 21,404  (74,138 pruned, 78%)
     hm3_rare     48,567 → 20,278  (28,289 pruned, 59%)
     nonhm3_rare  47,894 → 29,847  (18,047 pruned, 38%)
     HM3 common is most aggressively pruned — it was designed to tag LD blocks,
     so r²=0.05 is expensive. For covariate PCA (02a), consider r²=0.1 on HM3
     common only to recover ~10–20K more variants for PC stability. -->

## Cell 7 — PCA per (bin, source) (local)

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
plink2 --pfile {LOCAL}/r2_qc --extract {LOCAL}/{arm}.prune.in --nonfounders \
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

## Cell 8 — read participant scores

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]
kg_meta = pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]]


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

## Cell 9 — project 1000G into each space

```python
kg_scores = {}
for label in [l for l in scores if l.startswith("hm3_")]:
    wpath = f"{LOCAL}/pca/{label}.eigenvec.allele"
    apath = f"{LOCAL}/pca/{label}.acount"
    if not os.path.isfile(wpath):
        continue
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
grep -v '^##' {LOCAL}/r1_qc.pvar | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > {LOCAL}/{label}_ira.sorted
awk 'NR>1 {{print $2, $3, $4}}' "{KG_BFILE}.acount" | LC_ALL=C sort > {LOCAL}/kg_ira.sorted
LC_ALL=C comm -12 {LOCAL}/{label}_ira.sorted {LOCAL}/kg_ira.sorted | awk '{{print $1}}' > {LOCAL}/{label}_kg.ids
N=$(wc -l < {LOCAL}/{label}_kg.ids)
echo "{label}: $N variants shared with 1000G"
if [ "$N" -lt 100 ]; then exit 0; fi
W={wpath}
HEADER=$(head -1 "$W")
ID=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'ID'  | cut -d: -f1)
A1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'A1'  | cut -d: -f1)
P1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC1' | cut -d: -f1)
PK=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)
plink2 --bfile "{KG_BFILE}" --extract {LOCAL}/{label}_kg.ids --nonfounders \
  --read-freq {apath} \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${{P1}}-${{PK}}" --out {LOCAL}/kg_{label}
"""], check=True)
    kpath = f"{LOCAL}/kg_{label}.sscore"
    if os.path.isfile(kpath):
        kg_scores[label] = read_scores(kpath, "sample").merge(kg_meta, on="sample", how="left")
        print(f"{label}: {len(kg_scores[label]):,} 1000G samples projected")
```

## Cell 10 — do HM3 and non-HM3 agree?

```python
K_CMP = min(5, N_PCS_FIT)

for freq in ("rare", "common"):
    l_hm3, l_non = f"hm3_{freq}", f"nonhm3_{freq}"
    if l_hm3 not in scores or l_non not in scores:
        continue
    j = scores[l_hm3].merge(scores[l_non], on="person_id", suffixes=("_hm3", "_non"))
    r = [np.corrcoef(j[f"{p}_hm3"], j[f"{p}_non"])[0, 1] for p in PC[:K_CMP]]
    print(f"\n{freq}  |r| between HM3 and non-HM3 PCs:")
    print(pd.Series(np.abs(r), index=PC[:K_CMP]).round(3).to_string())
```

## Cell 11 — projection plots

Participants as a grey density layer; 1000G EUR populations coloured by pop.
One figure per arm (4 total), PC pairs 1-2 and 3-4 side by side.

```python
POP_COLORS = {"CEU": "#e6194b", "GBR": "#f58231", "FIN": "#3cb44b",
              "TSI": "#4363d8", "IBS": "#911eb4"}

PC_PAIRS = [(PC[0], PC[1]), (PC[2], PC[3])]

for label in groups:
    if label not in scores:
        continue
    s   = scores[label]
    kgs = kg_scores.get(label)

    fig, axes = plt.subplots(1, len(PC_PAIRS), figsize=(5 * len(PC_PAIRS), 4.5))
    for ax, (a, b) in zip(axes, PC_PAIRS):
        ax.hexbin(s[a], s[b], gridsize=80, cmap="Greys", mincnt=1,
                  bins="log", linewidths=0, zorder=1)
        if kgs is not None:
            for pop, color in POP_COLORS.items():
                sub = kgs[kgs["pop"] == pop]
                if sub.empty:
                    continue
                ax.scatter(sub[a], sub[b], s=18, color=color, label=pop,
                           edgecolors="white", linewidths=0.3, zorder=3)
        ax.set_xlabel(a); ax.set_ylabel(b)
    if kgs is not None:
        axes[0].legend(fontsize=8, markerscale=1.4, framealpha=0.7)
    suffix = "+ 1000G EUR" if kgs is not None else "(no 1000G projection)"
    fig.suptitle(f"{label.replace('_', ' ')} — participant cloud {suffix}")
    plt.tight_layout()
    plt.savefig(f"{OUT}/projection_{label}.png", dpi=130, bbox_inches="tight")
    plt.show()
```

## Cell 12 — loadings by genomic position

A PC carried by one region is an inversion or artefact. Check whether non-HM3
variants concentrate on specific chromosomes.

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

    n_show = min(K_CMP, 4)
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
