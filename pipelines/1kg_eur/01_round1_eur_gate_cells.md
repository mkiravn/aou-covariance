# 1kg_eur — round 1: 1000G European gate

Keeps AoU participants near the centroid of the five 1000G European
populations in AoU's ancestry PCs. Same gate as the ancestry survey's EUR
group; this writes it into the sample set's own folder.

**The distance.** Mahalanobis distance to the European reference centroid over
PCs 1–5, under the reference founders' own covariance: a multivariate normal
model of that cloud, so the gate follows its shape and orientation rather than
a sphere. The covariance is shrunk toward a sphere, since a 5×5 covariance from
a few hundred founders is otherwise noisy enough to stretch the gate along a
direction that is only sampling error. Width 1× is the distance that holds 99%
of the reference founders themselves, and wider gates scale that radius.

**No AoU label enters the selection.** From `ancestry_preds.tsv` this notebook
reads only `research_id` and `pca_features`. The centroid, the SDs and the
threshold all come from 1000G founders, labelled by 1000G's own sample file.
AoU's `ancestry_pred` is read once, in the last cell, after the keep list is
written, purely to compare. The PC space is AoU's, but it was fit on the
reference panel's genotypes; their labels only entered the classifier we
ignore.

The centroid uses founders only.

Compute: a small VM (4 vCPU, 16 GB).

---

## Cell 1 — config

```python
import os, ast, shutil, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

AUX = os.path.expanduser("~/workspace/cdrv9/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux")
ANCESTRY_PREDS = f"{AUX}/ancestry/ancestry_preds.tsv"
TRAINING_PCA   = f"{AUX}/ancestry/training_pca.tsv"

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
B = f"{WS}/phenotypic_covariance_v9"
SAMPLE_SET = "1kg_eur"
OUT = f"{B}/{SAMPLE_SET}/01_ancestry/round1"
os.makedirs(OUT, exist_ok=True)

KG_LABELS_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
                 "20130606_g1k_3202_samples_ped_population.txt")
KG_LABELS = f"{WS}/reference_labels/{os.path.basename(KG_LABELS_URL)}"

N_PCS     = 5
EUR_POPS  = ["CEU", "GBR", "FIN", "TSI", "IBS"]
WIDTH     = 1.0        # × the distance holding 99% of the reference founders

EUR_LOOSE_KEEP = f"{B}/01_ancestry_filtering/ancestry_pca_filter/final_pca/eur_loose/final_keep_ids_eur_loose_p50.txt"
EUR_D2_KEEP    = f"{B}/01_ancestry_filtering/r2_eur_gbr_ceu_projection/aou_D2_keep_ids.txt"

if not os.path.isfile(KG_LABELS):
    os.makedirs(os.path.dirname(KG_LABELS), exist_ok=True)
    tmp = f"/tmp/{os.path.basename(KG_LABELS)}"
    subprocess.run(["curl", "-fsSL", "--retry", "3", "-o", tmp, KG_LABELS_URL], check=True)
    shutil.copy(tmp, KG_LABELS)
for p in (ANCESTRY_PREDS, TRAINING_PCA):
    assert os.path.isfile(p), p
print(OUT)
```

## Cell 2 — load

```python
def parse_pcs(col, n=16):
    arr = np.vstack(col.map(ast.literal_eval).to_numpy())
    assert arr.shape[1] == n, arr.shape
    return arr

def show(n):
    return "≤20" if 0 < n <= 20 else f"{n:,}"

aou = pd.read_csv(ANCESTRY_PREDS, sep="\t", usecols=["research_id", "pca_features"])
X = parse_pcs(aou.pop("pca_features"))
aou = aou.rename(columns={"research_id": "person_id"})
aou["person_id"] = aou["person_id"].astype(str)

kg = pd.read_csv(KG_LABELS, sep=r"\s+", dtype=str)
kg = pd.DataFrame({"s": kg["SampleID"], "pop": kg["Population"],
                   "founder": (kg["FatherID"] == "0") & (kg["MotherID"] == "0")})

ref = pd.read_csv(TRAINING_PCA, sep="\t", usecols=["s", "scores"])
R = parse_pcs(ref.pop("scores"))
ref = ref.merge(kg, on="s", how="left", validate="one_to_one")
assert len(ref) == len(R)

is_eur = (ref["pop"].isin(EUR_POPS) & ref["founder"].eq(True)).to_numpy()
print(f"{len(aou):,} AoU participants, {len(ref):,} reference samples")
print(ref[ref["pop"].isin(EUR_POPS)].groupby("pop")
      .agg(n=("s", "size"), founders=("founder", "sum")).to_string())
assert is_eur.sum() > 2 * N_PCS
```

## Cell 3 — gate and keep list

`d99` is the distance holding 99% of the reference founders; the gate keeps
participants within `WIDTH` × that. Inputs here are the PC matrix and the
1000G labels only.

```python
Rg = R[is_eur, :N_PCS]
mu = Rg.mean(0)
try:
    from sklearn.covariance import LedoitWolf
    C = LedoitWolf().fit(Rg).covariance_
except ImportError:                      # same shrinkage target, fixed weight
    S_ = np.cov(Rg, rowvar=False)
    C = 0.9 * S_ + 0.1 * np.trace(S_) / S_.shape[0] * np.eye(S_.shape[0])
Cinv = np.linalg.inv(C)
dist = lambda Z: np.sqrt(np.einsum("ij,jk,ik->i", Z[:, :N_PCS] - mu, Cinv, Z[:, :N_PCS] - mu))
d_aou, d_ref = dist(X), dist(Rg)
d99 = np.quantile(d_ref, 0.99)

keep = d_aou <= WIDTH * d99
KEEP_PATH = f"{OUT}/{SAMPLE_SET}_round1_keep_ids.txt"
aou.loc[keep, "person_id"].to_csv(KEEP_PATH, index=False, header=False)

with open(f"{OUT}/round1_provenance.txt", "w") as f:
    f.write(f"reference\t1000G {', '.join(EUR_POPS)} founders in training_pca.tsv (n={int(is_eur.sum())})\n"
            f"labels\t{KG_LABELS_URL}\n"
            f"pcs\t1-{N_PCS} of ancestry_preds.tsv pca_features\n"
            f"distance\tMahalanobis over PCs 1-{N_PCS}, shrunk reference covariance\n"
            f"selection\tPC coordinates and 1000G labels only; no AoU label\n"
            f"width\t{WIDTH:g}x the distance holding 99% of the reference founders\n"
            f"threshold\t{WIDTH * d99:.6f}\nkept\t{int(keep.sum())}\n")
print(f"{show(int(keep.sum()))} kept -> {KEEP_PATH}")

mult = np.linspace(0.25, 3, 60)
ds = np.sort(d_aou / d99)
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(mult, np.maximum(np.searchsorted(ds, mult, side="right"), 1), color="k")
ax.axvline(WIDTH, color="tab:blue", ls="--", lw=1, label=f"{WIDTH:g}×")
ax.set_yscale("log")
ax.set_xlabel("gate width (× distance holding 99% of the reference)")
ax.set_ylabel("AoU participants inside")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/kept_vs_width.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4 — where the kept set sits

```python
from matplotlib.patches import Ellipse

rng = np.random.default_rng(0)
shown = rng.choice(np.flatnonzero(keep), size=min(30_000, int(keep.sum())), replace=False)
cols = dict(zip(EUR_POPS, plt.cm.tab10.colors))

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
for ax, (i, j) in zip(axes, [(0, 1), (2, 3)]):
    ax.hist2d(X[:, i], X[:, j], bins=250, cmap="Greys", norm=LogNorm(), zorder=0)
    for pop in EUR_POPS:
        m = (ref["pop"] == pop).to_numpy()
        ax.scatter(R[m, i], R[m, j], s=16, marker="x", color=cols[pop], label=pop, zorder=1)
    ax.scatter(X[shown, i], X[shown, j], s=2, alpha=0.3, color="tab:blue", rasterized=True, zorder=2)
    r = WIDTH * d99
    M = C[np.ix_([i, j], [i, j])]              # the ellipsoid's shadow on these two PCs
    vals, vecs = np.linalg.eigh(M)
    ax.add_patch(Ellipse((mu[i], mu[j]), 2 * r * np.sqrt(vals[-1]), 2 * r * np.sqrt(vals[0]),
                         angle=np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1])),
                         fill=False, edgecolor="tab:blue", lw=1.4, zorder=3))
    pts = np.r_[X[shown][:, [i, j]], R[is_eur][:, [i, j]]]
    pad = 0.2 * np.ptp(pts, axis=0)
    ax.set_xlim(pts[:, 0].min() - pad[0], pts[:, 0].max() + pad[0])
    ax.set_ylim(pts[:, 1].min() - pad[1], pts[:, 1].max() + pad[1])
    ax.set_xlabel(f"PC{i + 1}"); ax.set_ylabel(f"PC{j + 1}")
axes[0].legend(fontsize=8, markerscale=1.4)
fig.suptitle(f"{SAMPLE_SET} round 1 — kept (blue), 1000G Europeans (crosses), gate at {WIDTH:g}×")
plt.tight_layout()
plt.savefig(f"{OUT}/round1_pcs.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5 — comparison with AoU labels and earlier sets

The only cell that touches an AoU label, and it runs after the keep list is
written. Delete it and the sample set is unchanged.

```python
aou_labels = pd.read_csv(ANCESTRY_PREDS, sep="\t", usecols=["ancestry_pred"])["ancestry_pred"]
print(pd.crosstab(aou_labels.to_numpy(), np.where(keep, "kept", "not kept"))
      .apply(lambda c: c.map(show)).to_string())


def read_ids(path):
    s = pd.read_csv(path, sep=r"\s+", header=None, dtype=str).iloc[:, -1]
    return set(s[s.str.isdigit()])


kept_ids = set(aou.loc[keep, "person_id"])
for name, path in (("eur_loose", EUR_LOOSE_KEEP), ("eur_D2", EUR_D2_KEEP)):
    if not os.path.isfile(path):
        print(f"\n{name}: not found at {path}")
        continue
    other = read_ids(path)
    both = len(kept_ids & other)
    print(f"\n{name}: {show(len(other))}   in both {show(both)}   "
          f"only {SAMPLE_SET} {show(len(kept_ids - other))}   only {name} {show(len(other - kept_ids))}   "
          f"Jaccard {both / len(kept_ids | other):.3f}")
```

## Cell 6 — copy this notebook to the bucket

Save the notebook first, then run this.

```python
NOTEBOOK = os.path.expanduser(f"~/notebooks/{SAMPLE_SET}_round1.ipynb")
NB_DIR = f"{B}/{SAMPLE_SET}/notebooks"
os.makedirs(NB_DIR, exist_ok=True)
assert os.path.isfile(NOTEBOOK), NOTEBOOK
dest = f"{NB_DIR}/{os.path.basename(NOTEBOOK)}"
shutil.copy(NOTEBOOK, dest)
print(f"{'OK  ' if os.path.getsize(dest) == os.path.getsize(NOTEBOOK) else 'SIZE MISMATCH'} {dest}")
```
