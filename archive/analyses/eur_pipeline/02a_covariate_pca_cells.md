# 1kg_eur — covariate PCA

Produces the 20 covariate PCs for the round-2 EUR sample (n ≈ 223 K) following
Kemper et al. 2021: HapMap3 common variants, LD-pruned, PCA fit on the analysis
sample, scores exported.

Design decisions from variant-set exploration (01b):
- HM3 common only — non-HM3 arms are not projectable and add little
- r²=0.1 (not 0.05) — at r²=0.05 HM3 common loses 77% of variants; r²=0.1
  yields ~35–40K while still removing most LD
- QC re-run from `r1_qc` here (not reusing 01b's r2_qc): same filters but
  self-contained for reproducibility
- PCA runs locally; the approx algorithm handles 223K samples

Depends on: `r1_qc` pfile + round-2 keep list (from the round-2 gate notebook)

---

## Cell 1 — config

```python
import os, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SAMPLE_SET  = "1kg_eur"
CDR_VERSION = "v9"
CDR_DATASET = os.environ.get("WORKSPACE_CDR", "")   # fill in manually if blank

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
B  = f"{WS}/phenotypic_covariance_v9"

# r1_qc comes from the round-2 gate notebook (same scratch dir)
SCRATCH_GATE = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_round2")
R1_QC  = f"{SCRATCH_GATE}/r1_qc"
KEEP   = f"{B}/{SAMPLE_SET}/01_ancestry/round2/{SAMPLE_SET}_keep_ids.txt"

KG_DIR    = f"{WS}/1000g_reference"
KG_BFILE  = f"{KG_DIR}/1kg_all_qc"
KG_ACOUNT = f"{KG_DIR}/1kg_all_qc.acount"
KG_PANEL  = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"
EUR_POPS    = ["CEU", "GBR", "FIN", "TSI", "IBS"]
ANCHOR_POPS = ["CEU", "GBR"]

OUT   = f"{B}/{SAMPLE_SET}/01_ancestry/covariate_pca"
LOCAL = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_covpca")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOCAL, exist_ok=True)

PRUNE           = "1000kb 1 0.1"
N_PCS_FIT       = 20
N_PCS_COVARIATE = 20

os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"

assert os.path.isfile(f"{R1_QC}.pvar"), f"{R1_QC} — run the round-2 gate notebook first"
assert os.path.isfile(KEEP),            f"{KEEP} — run the round-2 gate notebook first"
print(f"r1_qc:  {R1_QC}")
print(f"keep:   {sum(1 for _ in open(KEEP)):,} participants")
print(f"output: {OUT}")
```

## Cell 2 — install plink2

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

## Cell 3 — site QC in round-2 cohort and HM3 ID list

Filters `r1_qc` to the round-2 keep list, reapplies MAF / HWE / missingness in
the r2 cohort, then identifies HM3-agreeing variants by ID+REF+ALT match against
the 1000G allele count file.

```python
r2_qc = f"{LOCAL}/r2_qc"
hm3_ids = f"{LOCAL}/hm3.ids"

if not os.path.isfile(f"{r2_qc}.psam"):
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
plink2 --pfile "{R1_QC}" --keep "{KEEP}" --nonfounders \
  --maf 0.01 --hwe 1e-6 0 keep-fewhet --geno 0.05 \
  --threads $(nproc) --make-pgen --out "{r2_qc}"
echo "r2_qc: $(wc -l < "{r2_qc}.psam") samples, $(grep -vc '^##' "{r2_qc}.pvar") variants"
"""], check=True)
else:
    print(f"r2_qc exists: {sum(1 for _ in open(r2_qc + '.psam')) - 1:,} samples")

if not os.path.isfile(hm3_ids):
    subprocess.run(["bash", "-c", f"""
set -eo pipefail
grep -v '^##' "{r2_qc}.pvar" | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > "{LOCAL}/r2_ira.sorted"
awk 'NR>1 {{print $2, $3, $4}}' "{KG_ACOUNT}" | LC_ALL=C sort > "{LOCAL}/kg_ira.sorted"
LC_ALL=C comm -12 "{LOCAL}/r2_ira.sorted" "{LOCAL}/kg_ira.sorted" | awk '{{print $1}}' > "{hm3_ids}"
echo "HM3-agreeing: $(wc -l < "{hm3_ids}") variants"
"""], check=True)
else:
    print(f"hm3.ids exists: {sum(1 for _ in open(hm3_ids)):,} variants")
```

## Cell 4 — prune HM3 common variants (r²=0.1)

```python
HIGH_LD_REGIONS_GRCH38 = """\
chr1 47761740 51761740 1
chr2 85919365 100517106 2
chr2 182427027 189427029 3
chr3 47483505 49987563 4
chr3 83368158 86868160 5
chr5 44464140 51168409 6
chr5 129636407 132636409 7
chr6 25391792 33424245 8
chr6 57788603 58453888 9
chr6 61109122 61357029 10
chr6 139637169 142137170 11
chr7 54964812 66897578 12
chr8 8105067 12105082 13
chr8 43025699 48924888 14
chr8 110918594 113918595 15
chr10 36671065 43184546 16
chr11 88127183 91127184 17
chr12 32955798 41319931 18
chr20 33948532 36438183 19
"""
LD_REGIONS = f"{LOCAL}/high_ld_regions.txt"
open(LD_REGIONS, "w").write(HIGH_LD_REGIONS_GRCH38)

prune_out = f"{LOCAL}/hm3_common_covpca"
subprocess.run(["bash", "-c", f"""
set -e
plink2 --pfile "{r2_qc}" --nonfounders \
  --extract "{hm3_ids}" \
  --maf 0.01 \
  --exclude bed1 "{LD_REGIONS}" \
  --indep-pairwise {PRUNE} \
  --threads $(nproc) --out "{prune_out}"
echo "HM3 common after pruning (r²=0.1): $(wc -l < "{prune_out}.prune.in") variants"
"""], check=True)
```

## Cell 5 — fit PCA (local, approx)

Uses `--pca approx` with allele weights so participants can be re-scored through
the loadings later (and so 1000G can be projected).

```python
pca_out = f"{LOCAL}/covpca"
subprocess.run(["bash", "-c", f"""
set -e
plink2 --pfile "{r2_qc}" --nonfounders \
  --extract "{prune_out}.prune.in" \
  --freq counts \
  --pca approx {N_PCS_FIT} allele-wts \
  --threads $(nproc) --out "{pca_out}"
echo "PCA done; eigenvalues: {pca_out}.eigenval"
"""], check=True)
```

## Cell 6 — score participants through loadings

Re-scoring through the fit's own loadings (not `--pca` direct scores) means 1000G
can be projected using the same command.

```python
PC = [f"PC{k}" for k in range(1, N_PCS_FIT + 1)]

def read_scores(path, id_col):
    d = pd.read_csv(path, sep=r"\s+")
    col = "#IID" if "#IID" in d.columns else "IID"
    d = d.rename(columns={col: id_col,
                           **{f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS_FIT + 1)}})
    d[id_col] = d[id_col].astype(str)
    return d[[id_col] + PC]

W   = f"{pca_out}.eigenvec.allele"
_h  = open(W).readline().split()
_id = _h.index("#ID") + 1 if "#ID" in _h else _h.index("ID") + 1
_a1 = _h.index("A1") + 1
_p1 = _h.index("PC1") + 1
_pk = _h.index(f"PC{N_PCS_FIT}") + 1

score_out = f"{LOCAL}/part_covpca"
subprocess.run(["bash", "-c", f"""
set -e
plink2 --pfile "{r2_qc}" --nonfounders \
  --extract "{prune_out}.prune.in" \
  --read-freq "{pca_out}.acount" \
  --score "{W}" {_id} {_a1} header-read no-mean-imputation variance-standardize \
  --score-col-nums {_p1}-{_pk} \
  --threads $(nproc) --out "{score_out}"
echo "scored: $(wc -l < "{score_out}.sscore") lines"
"""], check=True)

part_scores = read_scores(f"{score_out}.sscore", "person_id")
eigen       = np.loadtxt(f"{pca_out}.eigenval")
print(f"{len(part_scores):,} participants, {N_PCS_FIT} PCs")
print(f"variance explained: PC1 {eigen[0]/eigen.sum()*100:.2f}%  "
      f"PC2 {eigen[1]/eigen.sum()*100:.2f}%  "
      f"PC20 {eigen[19]/eigen.sum()*100:.2f}%")
```

## Cell 7 — project 1000G EUR

```python
kg_out = f"{LOCAL}/kg_covpca"
subprocess.run(["bash", "-c", f"""
set -e
# variants shared with 1000G (ID+REF+ALT match)
awk 'NR>1 {{print $2, $3, $4}}' "{pca_out}.acount" \
  | LC_ALL=C sort > "{LOCAL}/covpca_ira.sorted"
awk 'NR>1 {{print $2, $3, $4}}' "{KG_BFILE}.acount" \
  | LC_ALL=C sort > "{LOCAL}/kg_ira.sorted"
LC_ALL=C comm -12 "{LOCAL}/covpca_ira.sorted" "{LOCAL}/kg_ira.sorted" \
  | awk '{{print $1}}' > "{LOCAL}/covpca_kg.ids"
N=$(wc -l < "{LOCAL}/covpca_kg.ids")
echo "shared with 1000G: $N variants"

plink2 --bfile "{KG_BFILE}" --nonfounders \
  --extract "{LOCAL}/covpca_kg.ids" \
  --read-freq "{pca_out}.acount" \
  --score "{W}" {_id} {_a1} header-read no-mean-imputation variance-standardize \
  --score-col-nums {_p1}-{_pk} \
  --out "{kg_out}"
echo "projected: $(grep -c . "{kg_out}.sscore") 1000G samples"
"""], check=True)

kg_scores = read_scores(f"{kg_out}.sscore", "sample").merge(
    pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]],
    on="sample", how="left")
print(kg_scores.groupby("pop").size().to_string())
```

## Cell 8 — scree + EUR population separation

```python
POP_COLORS = {"CEU": "#e6194b", "GBR": "#f58231", "FIN": "#3cb44b",
              "TSI": "#4363d8", "IBS": "#911eb4"}
EUR_COLOR  = "#aaaaaa"

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# scree
axes[0].plot(range(1, N_PCS_FIT + 1), eigen / eigen.sum() * 100,
             marker="o", ms=4, color="steelblue")
axes[0].set_xlabel("PC"); axes[0].set_ylabel("% variance explained")
axes[0].set_title("scree")

# PC1 vs PC2
eur_kg = kg_scores[kg_scores["super_pop"] == "EUR"]
axes[1].hexbin(part_scores["PC1"], part_scores["PC2"],
               gridsize=80, cmap="Greys", mincnt=1, bins="log",
               linewidths=0, zorder=1)
for pop, col in POP_COLORS.items():
    sub = eur_kg[eur_kg["pop"] == pop]
    axes[1].scatter(sub["PC1"], sub["PC2"], s=10, color=col,
                    label=pop, alpha=0.8, linewidths=0, zorder=3)
axes[1].set_xlabel("PC1"); axes[1].set_ylabel("PC2")
axes[1].set_title("PC1 vs PC2 — 1000G EUR overlaid")
axes[1].legend(fontsize=8, markerscale=1.5)

# EUR subpop separation: mean diff between anchors and others in participant SD
sep = []
anch  = eur_kg["pop"].isin(ANCHOR_POPS)
other = eur_kg["pop"].isin(set(EUR_POPS) - set(ANCHOR_POPS))
for p in PC:
    d = abs(eur_kg.loc[anch, p].mean() - eur_kg.loc[other, p].mean())
    sep.append(d / part_scores[p].std())
axes[2].bar(range(1, N_PCS_FIT + 1), sep, color="steelblue", alpha=0.8)
axes[2].set_xlabel("PC"); axes[2].set_ylabel("anchor separation (participant SD)")
axes[2].set_title(f"CEU+GBR vs FIN+TSI+IBS")

plt.suptitle(f"{SAMPLE_SET} covariate PCA — HM3 common, r²=0.1, n={len(part_scores):,}", y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/covpca_overview.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 9 — PC pairs coloured by 1000G EUR subpopulation

```python
PC_PAIRS = [(PC[i], PC[i+1]) for i in range(0, min(8, N_PCS_FIT - 1), 2)]
eur_kg   = kg_scores[kg_scores["super_pop"] == "EUR"]

ncols = len(PC_PAIRS)
fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 4.8))
for ax, (a, b) in zip(np.atleast_1d(axes), PC_PAIRS):
    ax.hexbin(part_scores[a], part_scores[b],
              gridsize=70, cmap="Greys", mincnt=1, bins="log",
              linewidths=0, zorder=1)
    for pop, col in POP_COLORS.items():
        sub = eur_kg[eur_kg["pop"] == pop]
        ax.scatter(sub[a], sub[b], s=12, color=col, label=pop,
                   alpha=0.85, linewidths=0, zorder=3)
    ax.set_xlabel(a); ax.set_ylabel(b)
np.atleast_1d(axes)[0].legend(fontsize=7, markerscale=1.5)
plt.suptitle(f"{SAMPLE_SET} — PC pairs 1–8 with 1000G EUR", y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/covpca_pc_pairs.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 10 — LD peak check

A PC dominated by a single region is LD, not ancestry structure. Flag outlier
loadings, widen each to a ±250 kb window, and re-add to the exclusion list.
If anything is flagged: add those regions to Cell 1's LD_REGIONS, rerun Cells
3–8.

```python
PEAK_Q, FLANK_KB, PEAK_PCS = 0.9995, 250, min(10, N_PCS_FIT)

L = pd.read_csv(f"{pca_out}.eigenvec.allele", sep=r"\s+")
idc = "#ID" if "#ID" in L.columns else "ID"
L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
L["POS"] = L["POS"].astype(int)
hit = np.zeros(len(L), dtype=bool)
for p in PC[:PEAK_PCS]:
    v = L[p].to_numpy() ** 2
    hit |= v > np.quantile(v, PEAK_Q)
print(f"{hit.sum()} loading outliers across PC1-{PEAK_PCS}")

if hit.any():
    H = L.loc[hit, ["CHROM", "POS"]].sort_values(["CHROM", "POS"])
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
    print(f"\n{len(regions)} peak regions — add to LD_REGIONS and rerun Cells 3–8:")
    print("\n".join(regions))
else:
    print("no peaks; PCs look clean")
```

## Cell 11 — loadings by genomic position

```python
n_show = min(4, N_PCS_FIT)
L = pd.read_csv(f"{pca_out}.eigenvec.allele", sep=r"\s+")
idc = "#ID" if "#ID" in L.columns else "ID"
L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
L["POS"] = L["POS"].astype(int)
L2 = L.copy()
L2 = L2[L2["CHROM"].isin([str(c) for c in range(1, 23)])].copy()
L2["CHROM"] = pd.Categorical(L2["CHROM"], categories=[str(c) for c in range(1, 23)], ordered=True)
L2 = L2.sort_values(["CHROM", "POS"])

offset, centres = 0, {}
cum = np.empty(len(L2))
for c, idx in L2.groupby("CHROM", observed=True).indices.items():
    p = L2["POS"].to_numpy()[idx]
    cum[idx] = p + offset
    centres[c] = offset + p.max() / 2
    offset += p.max()
L2["CUM"] = cum

fig, axs = plt.subplots(n_show, 1, figsize=(14, 2.4 * n_show), sharex=True)
for ax, pc in zip(np.atleast_1d(axs), PC[:n_show]):
    for k, c in enumerate(centres):
        s = L2[L2["CHROM"] == c]
        ax.scatter(s["CUM"], s[pc] ** 2, s=2, alpha=0.5,
                   color=["tab:blue", "tab:orange"][k % 2], rasterized=True)
    ax.set_ylabel(f"{pc} loading²")
np.atleast_1d(axs)[-1].set_xticks(list(centres.values()))
np.atleast_1d(axs)[-1].set_xticklabels(list(centres), fontsize=7)
np.atleast_1d(axs)[-1].set_xlabel("chromosome")
fig.suptitle(f"{SAMPLE_SET} covariate PCA — loadings")
plt.tight_layout()
plt.savefig(f"{OUT}/covpca_loadings.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 12 — batch effect check

R² of each PC on metadata factors — requires a metadata table with `person_id`
and batch/site columns. Set `META` to a file path when one becomes available
(e.g. from BigQuery `person.data_partner_id`).

```python
META = None         # e.g. path to a TSV with person_id, data_partner_id, ...
META_FACTORS = ["data_partner_id", "site_id", "sequencing_center"]

if META and os.path.isfile(META):
    m = pd.read_csv(META, sep="\t", dtype=str).rename(columns={"research_id": "person_id"})
    j = part_scores.merge(m, on="person_id", how="inner")
    rows = {}
    for f in META_FACTORS:
        if f not in j:
            continue
        rows[f] = [float(1 - j.groupby(f)[p].transform("var").mean() / j[p].var())
                   for p in PC[:10]]
    if rows:
        print("R² of PC1–10 on batch factors:")
        print(pd.DataFrame(rows, index=PC[:10]).round(3).to_string())
else:
    print("no metadata table available — skipping batch check")
```

<!-- Run notes (2026-09-24, n=223,209):
     r2_qc: 1,183,682 variants after MAF/HWE/geno in r2 cohort
     HM3 common after pruning (r²=0.1, 1000kb, exclusion zones): 64,379 variants
     Variance explained: PC1 16.90%, PC2 13.95%, ..., PC20 3.38%
     Informative PCs: 1–5 show signal; 6–20 noise
     Loadings plots: no obvious LD peaks in PC1–4; loadings broadly distributed
     Recommendation: use PC1–5 (or 1–10 conservatively) as GWAS covariates -->

## Cell 13 — write covariate PCs

```python
cov = part_scores[["person_id"] + PC[:N_PCS_COVARIATE]].copy()
cov = cov.rename(columns={"person_id": "IID"})
COV_PATH = f"{OUT}/final_pca_pc_covariates_{SAMPLE_SET}.txt"
cov.to_csv(COV_PATH, sep="\t", index=False)

with open(f"{OUT}/covariate_pca_provenance.txt", "w") as fh:
    fh.write(
        f"input\tr2_qc pfile (01b: MAF/HWE/geno in r2 cohort)\n"
        f"variants\tHM3 common (maf≥0.01, HM3 ID+REF+ALT match)\n"
        f"prune\t--indep-pairwise {PRUNE}, long-range LD regions excluded (bed1)\n"
        f"n_variants_pruned\t{sum(1 for _ in open(prune_out + '.prune.in')):,}\n"
        f"pca\t--pca approx {N_PCS_FIT} allele-wts, fit on r2 participants\n"
        f"scores\t--score through fit loadings, --read-freq its .acount\n"
        f"covariates\tPC1-{N_PCS_COVARIATE}\n"
        f"n_participants\t{len(cov):,}\n"
    )
print(f"{len(cov):,} participants, PC1-{N_PCS_COVARIATE} → {COV_PATH}")
```

---

Next: residualization reads `final_pca_pc_covariates_{SAMPLE_SET}.txt`.
