# 1kg_eur — round 2: CEU + GBR gate

Refits PCA on the round-1 participants, projects 1000G into that space, and
keeps those closest to the CEU + GBR centroid. Two alternative gates are built
for comparison: the participants' own centroid (as eur_D2 did), and the Kemper
direction, PCA fit on the 1000G Europeans with participants projected in.

The gate is a multivariate normal model of the anchor population, assembled
from the part of each estimate that is trustworthy: **centre** and **shape**
from the projected CEU + GBR samples, **scale** from the participants, **size**
from a target. Projection shrinks each PC by its own factor, which biases the
anchor's variances but cancels in its correlations, so the correlation matrix is
taken from the anchor and the per-PC variances from the participants.

Participants are uncorrelated on their own PCs by construction, so any tilt in
the gate comes from the anchor's shape, not theirs. Cell 7 checks both.

Compute: 16 vCPU, ~100 GB RAM, ~500 GB disk. The PCA fit is the long step.

---

## Cell 1 — config

```python
import os, shlex, shutil, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SAMPLE_SET = "1kg_eur"
WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
B  = f"{WS}/phenotypic_covariance_v9"
CDR_VERSION = "v9"

ROUND1_KEEP = f"{B}/{SAMPLE_SET}/01_ancestry/round1/{SAMPLE_SET}_round1_keep_ids.txt"
PANEL       = f"{B}/01_ancestry_filtering/unified_panel/unified_panel_{CDR_VERSION}"
KG_DIR      = f"{WS}/1000g_reference"
KG_BFILE    = f"{KG_DIR}/1kg_all_qc"          # HM3-restricted, QC'd 1000G
KG_PANEL    = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

OUT = f"{B}/{SAMPLE_SET}/01_ancestry/round2"
LOCAL = os.path.expanduser(f"~/scratch_{SAMPLE_SET}_round2")
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOCAL, exist_ok=True)

ANCHOR_POPS = ["CEU", "GBR"]
EUR_POPS    = ["CEU", "GBR", "FIN", "TSI", "IBS"]
K_PCS       = 4                 # set after the scree in Cell 6
N_PCS_FIT   = 20
TARGET_N    = 221_992           # eur_D2's size
N_THREADS   = os.cpu_count()

assert os.path.isfile(ROUND1_KEEP), f"{ROUND1_KEEP} -- run round 1 first"
print(f"round 1: {sum(1 for _ in open(ROUND1_KEEP)):,} participants")
assert os.path.isfile(f"{KG_BFILE}.bed"), KG_BFILE
assert os.path.isfile(KG_PANEL), KG_PANEL


def sh(cmd):
    subprocess.run(["bash", "-c", "set -eo pipefail\n" + cmd], check=True)


def q(p):
    return shlex.quote(str(p))


print(OUT)
```

## Cell 2 — plink2

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
nproc; free -h; df -h "$HOME"
```

```python
os.environ["PATH"] = f"{os.path.expanduser('~/bin')}:{os.environ['PATH']}"
```

## Cell 3 — stage the panel locally

gcsfuse reads of large files have corrupted before; copy and check sizes.

```python
for ext in ("pgen", "pvar", "psam"):
    src, dst = f"{PANEL}.{ext}", f"{LOCAL}/panel.{ext}"
    if not os.path.isfile(dst) or os.path.getsize(dst) != os.path.getsize(src):
        shutil.copy(src, dst)
    assert os.path.getsize(dst) == os.path.getsize(src), dst
    print(f"{ext}: {os.path.getsize(dst):,}")
```

## Cell 4 — variants for the PCA

Round-1 participants only; QC within them, keep the HM3 variants that agree
with 1000G on ID+REF+ALT, drop long-range LD regions, then prune.

`--hwe 1e-6 0 keep-fewhet`: a fixed threshold (the `0` is the sample-size
term, which plink2 requires at this n), removing only heterozygote excess.
Heterozygote deficit here is population structure, not error.

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
LD_REGIONS = f"{LOCAL}/high_ld_regions_grch38.txt"
open(LD_REGIONS, "w").write(HIGH_LD_REGIONS_GRCH38)

sh(f"""
plink2 --pfile {q(f'{LOCAL}/panel')} --keep {q(ROUND1_KEEP)} --nonfounders \
  --maf 0.01 --hwe 1e-6 0 keep-fewhet --geno 0.05 --max-alleles 2 --rm-dup exclude-all \
  --threads {N_THREADS} --make-pgen --out {q(f'{LOCAL}/r1_qc')}

grep -v '^##' {q(f'{LOCAL}/r1_qc.pvar')} | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > {q(f'{LOCAL}/r1_ira.sorted')}
awk 'NR>1 {{print $2, $3, $4}}' {q(f'{KG_BFILE}.acount')} | LC_ALL=C sort > {q(f'{LOCAL}/kg_ira.sorted')}
LC_ALL=C comm -12 {q(f'{LOCAL}/r1_ira.sorted')} {q(f'{LOCAL}/kg_ira.sorted')} | awk '{{print $1}}' > {q(f'{LOCAL}/hm3_agreeing.ids')}
echo "HM3-agreeing: $(wc -l < {q(f'{LOCAL}/hm3_agreeing.ids')})"

plink2 --pfile {q(f'{LOCAL}/r1_qc')} --extract {q(f'{LOCAL}/hm3_agreeing.ids')} \
  --exclude range {q(LD_REGIONS)} --nonfounders --indep-pairwise 1000 50 0.1 \
  --threads {N_THREADS} --out {q(f'{LOCAL}/prune')}

plink2 --pfile {q(f'{LOCAL}/r1_qc')} --extract {q(f'{LOCAL}/prune.prune.in')} \
  --threads {N_THREADS} --make-pgen --out {q(f'{LOCAL}/pca_input')}
""")
print("PCA variants:", sum(1 for _ in open(f"{LOCAL}/prune.prune.in")))
```

## Cell 5 — fit the PCA on the participants

```python
sh(f"""
plink2 --pfile {q(f'{LOCAL}/pca_input')} --nonfounders --freq counts \
  --pca approx {N_PCS_FIT} allele-wts --threads {N_THREADS} \
  --out {q(f'{LOCAL}/round2_pca')}
""")
sh(f"cp {q(f'{LOCAL}/round2_pca.eigenval')} {q(f'{LOCAL}/round2_pca.eigenvec.allele')} "
   f"{q(f'{LOCAL}/round2_pca.log')} {q(OUT)}/")
print(open(f"{LOCAL}/round2_pca.eigenval").read())
```

## Cell 6 — score everyone through the same loadings, pick K

plink's `--pca` eigenvectors and a `--score` projection do not land on the
same coordinates, so participants are re-scored through their own loadings
rather than read from the eigenvec. Both sets then sit on identical axes.

```python
def score_cmd(bfile_flag, prefix, weights, freq, out, extract=None):
    return f"""
W={q(weights)}
HEADER=$(head -1 "$W")
ID=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'ID' | cut -d: -f1)
A1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'A1' | cut -d: -f1)
P1=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC1' | cut -d: -f1)
PK=$(echo "$HEADER" | tr '\\t' '\\n' | grep -nx 'PC{N_PCS_FIT}' | cut -d: -f1)
plink2 {bfile_flag} {q(prefix)} {'--extract ' + q(extract) if extract else ''} \
  --read-freq {q(freq)} \
  --score "$W" "$ID" "$A1" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${{P1}}-${{PK}}" --threads {N_THREADS} --out {q(out)}
"""


sh(f"""
grep -v '^##' {q(f'{LOCAL}/pca_input.pvar')} | awk 'NR>1 {{print $3, $4, $5}}' | LC_ALL=C sort > {q(f'{LOCAL}/pca_ira.sorted')}
LC_ALL=C comm -12 {q(f'{LOCAL}/pca_ira.sorted')} {q(f'{LOCAL}/kg_ira.sorted')} | awk '{{print $1}}' > {q(f'{LOCAL}/kg_project.ids')}
echo "variants for the 1000G projection: $(wc -l < {q(f'{LOCAL}/kg_project.ids')})"
""")
sh(score_cmd("--bfile", KG_BFILE, f"{LOCAL}/round2_pca.eigenvec.allele",
             f"{LOCAL}/round2_pca.acount", f"{LOCAL}/kg_in_r2", f"{LOCAL}/kg_project.ids"))
sh(score_cmd("--pfile", f"{LOCAL}/pca_input", f"{LOCAL}/round2_pca.eigenvec.allele",
             f"{LOCAL}/round2_pca.acount", f"{LOCAL}/part_in_r2"))

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

## Cell 6b — loadings

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

## Cell 7 — the gate

Centroid from the projected anchor populations, scale from the participants.

```python
SHRINK = 0.2            # anchor correlations toward independence; a few hundred samples

USE = PC[:K_PCS]
anchor = kg[kg["pop"].isin(ANCHOR_POPS)]
mu = anchor[USE].mean().to_numpy()
sd = part[USE].std(ddof=1).to_numpy()

off_part = np.abs(np.corrcoef(part[USE].to_numpy(), rowvar=False) - np.eye(K_PCS)).max()
Ra = np.corrcoef(anchor[USE].to_numpy(), rowvar=False)
off_anchor = np.abs(Ra - np.eye(K_PCS)).max()
print(f"largest off-diagonal correlation: participants {off_part:.3f}, anchor {off_anchor:.3f} "
      f"({len(anchor)} anchor samples)")
assert off_part < 0.05, "participant PCs are correlated; the PCA fit sample is not what we think"

Ra = (1 - SHRINK) * Ra + SHRINK * np.eye(K_PCS)
C = np.outer(sd, sd) * Ra            # anchor shape, participant scale
Cinv = np.linalg.inv(C)


def distance(frame, centre):
    d = frame[USE].to_numpy() - centre
    return np.sqrt(np.einsum("ij,jk,ik->i", d, Cinv, d))


d_anchor = distance(part, mu)
d_sorted = np.sort(d_anchor)

scan = pd.DataFrame({"target_n": np.linspace(0.8 * TARGET_N, 1.2 * TARGET_N, 5).astype(int)})
scan["threshold"] = d_sorted[scan["target_n"] - 1]
print(scan.round(4).to_string(index=False))

THRESHOLD = float(d_sorted[TARGET_N - 1])
keep = d_anchor <= THRESHOLD
KEEP_PATH = f"{OUT}/{SAMPLE_SET}_keep_ids.txt"
part.loc[keep, "person_id"].to_csv(KEEP_PATH, index=False, header=False)

with open(f"{OUT}/round2_provenance.txt", "w") as f:
    f.write(f"pca\tfit on the round-1 set, plink2 --pca approx {N_PCS_FIT} allele-wts\n"
            f"variants\t{sum(1 for _ in open(f'{LOCAL}/prune.prune.in'))} pruned HM3 sites\n"
            f"anchor\t1000G {'+'.join(ANCHOR_POPS)} projected into that space\n"
            f"pcs\t1-{K_PCS}\n"
            f"metric\tMahalanobis: anchor correlations (shrunk {SHRINK:g}), participant SDs\n"
            f"threshold\t{THRESHOLD:.6f}\nkept\t{int(keep.sum())}\n")
print(f"{int(keep.sum()):,} kept -> {KEEP_PATH}")
```

## Cell 8 — plots

```python
from matplotlib.patches import Ellipse

PAIRS = [(0, 1), (2, 3)]
rng = np.random.default_rng(0)
shown = rng.choice(np.flatnonzero(keep), size=min(50_000, int(keep.sum())), replace=False)
cols = dict(zip(EUR_POPS, plt.cm.tab10.colors))


def gate_ellipse(ax, centre, i, j, radius, **kw):
    """Gate boundary on PCs i, j: the ellipsoid's shadow, from the 2x2 marginal
    covariance, so it tilts with the anchor's correlation structure."""
    M = C[np.ix_([i, j], [i, j])]
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

## Cell 9 — alternative gates

Same PC space and the same size, centred on the participants instead of the
anchor; and the Kemper direction, PCA fit on the 1000G Europeans with
participants projected in.

```python
# a. participants' own centroid, as eur_D2 did
d_self = distance(part, part[USE].mean().to_numpy())
keep_self = d_self <= np.sort(d_self)[TARGET_N - 1]

# b. Kemper direction: fit on the 1000G Europeans, project participants
eur_ids = pd.read_csv(KG_PANEL, sep=r"\s+")
eur_ids = eur_ids.loc[eur_ids["super_pop"].eq("EUR"), "sample"]
eur_ids.to_csv(f"{LOCAL}/kg_eur.keep", index=False, header=False)

sh(f"""
plink2 --bfile {q(KG_BFILE)} --keep {q(f'{LOCAL}/kg_eur.keep')} \
  --extract {q(f'{LOCAL}/kg_project.ids')} --nonfounders --freq counts \
  --pca {N_PCS_FIT} allele-wts --threads {N_THREADS} --out {q(f'{LOCAL}/kgeur_pca')}
""")
sh(score_cmd("--pfile", f"{LOCAL}/pca_input", f"{LOCAL}/kgeur_pca.eigenvec.allele",
             f"{LOCAL}/kgeur_pca.acount", f"{LOCAL}/part_in_kgeur", f"{LOCAL}/kg_project.ids"))

part_k = read_scores(f"{LOCAL}/part_in_kgeur.sscore", "person_id")
kg_k = pd.read_csv(f"{LOCAL}/kgeur_pca.eigenvec", sep=r"\s+")
kg_k = kg_k.rename(columns={("#IID" if "#IID" in kg_k.columns else "IID"): "sample"})
kg_k = kg_k.merge(pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop"]], on="sample", how="left")

mu_k = kg_k[kg_k["pop"].isin(ANCHOR_POPS)][USE].mean().to_numpy()
sd_k = part_k[USE].std(ddof=1).to_numpy()
Rk = np.corrcoef(kg_k[kg_k["pop"].isin(ANCHOR_POPS)][USE].to_numpy(), rowvar=False)
Ck_inv = np.linalg.inv(np.outer(sd_k, sd_k) * ((1 - SHRINK) * Rk + SHRINK * np.eye(K_PCS)))
dk = part_k[USE].to_numpy() - mu_k
d_k = np.sqrt(np.einsum("ij,jk,ik->i", dk, Ck_inv, dk))
keep_kemper = pd.Series(d_k <= np.sort(d_k)[TARGET_N - 1], index=part_k["person_id"]).reindex(
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
        gate_ellipse(ax, part[USE].mean().to_numpy(), i, j, np.sort(d_self)[TARGET_N - 1],
                     edgecolor="tab:orange", lw=1.4, ls="--", zorder=4)
    ax.set_xlabel(a); ax.set_ylabel(b)
axes[0].legend(fontsize=8, markerscale=4)
fig.suptitle("who the anchor changes, against the participants' own centroid")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_difference.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 10 — copy this notebook to the bucket

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

Next: `04_final_pca.ipynb` style refit on this keep list for the covariate PCs,
then the GRM panel QC.
