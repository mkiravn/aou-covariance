# Ancestry survey

For each reference group of 1000 Genomes or HGDP populations:

1. gate AoU participants around the group's centroid in AoU's ancestry PCs, at several widths, and save the lists
2. test whether the group is a cluster or a continuum
3. plot each gate in PC space
4. measure overlap between gates
5. count related pairs inside each gate, by kinship class

Population labels come from the projects' own sample files. Kinship comes from
AoU's `samples_relatedness.tsv` (Hail PC-Relate, pairs above 0.1 only).
Counts of 1–20 are shown as `≤20`.

Compute: a small VM (4 vCPU, 16 GB).

---

## Cell 1 — config

Paths, gate populations and widths. Downloads the reference label files once.

```python
import os, ast, shutil, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

AUX = os.path.expanduser("~/workspace/cdrv9/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux")
ANCESTRY_PREDS = f"{AUX}/ancestry/ancestry_preds.tsv"
TRAINING_PCA   = f"{AUX}/ancestry/training_pca.tsv"
REL_PATH       = f"{AUX}/relatedness/samples_relatedness.tsv"
REL_COLS       = ("i.s", "j.s", "kin")     # participant 1, participant 2, kinship coefficient

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
OUT = f"{WS}/phenotypic_covariance_v9/analyses/ancestry_survey"
KEEP_DIR = f"{OUT}/keep"
os.makedirs(KEEP_DIR, exist_ok=True)

LABEL_DIR = f"{WS}/reference_labels"
LABEL_FILES = {
    "1kg":  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
            "20130606_g1k_3202_samples_ped_population.txt",
    "hgdp": "https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/"
            "hgdp_wgs.20190516.metadata.txt",
}

N_PCS = 5
WIDTHS = [0.5, 1.0, 2.0]                     # multiples of the distance holding 99% of the reference
W_NARROW, W_WIDE = min(WIDTHS), max(WIDTHS)
KIN_CUTS = (0.177, 0.354)                    # KING first-degree bounds, kinship scale

GATES = {
    "EUR": ["CEU", "GBR", "FIN", "TSI", "IBS"],                    # 1000G
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN"],                    # 1000G, without ASW, ACB
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],                    # 1000G
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],                    # 1000G
    "MID": ["Druze", "Bedouin", "Palestinian", "Mozabite"],        # HGDP
    "NAT": ["Maya", "Pima", "Karitiana", "Surui", "Colombian"],    # HGDP, Indigenous American
    "AMR": ["MXL", "PUR", "CLM", "PEL"],                           # 1000G, admixed
}
GCOL = dict(zip(GATES, ["tab:blue", "tab:orange", "tab:green", "tab:purple",
                        "tab:brown", "tab:olive", "tab:red"]))

os.makedirs(LABEL_DIR, exist_ok=True)
for url in LABEL_FILES.values():
    path = f"{LABEL_DIR}/{os.path.basename(url)}"
    if not os.path.isfile(path):
        tmp = f"/tmp/{os.path.basename(url)}"
        subprocess.run(["curl", "-fsSL", "--retry", "3", "-o", tmp, url], check=True)
        shutil.copy(tmp, path)
for p in (ANCESTRY_PREDS, TRAINING_PCA):
    assert os.path.isfile(p), p
print(OUT)
```

## Cell 2 — load

AoU and reference PCs, with the original population labels. Stops if a gate
population is missing. Only founders define the gates.

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

kg = pd.read_csv(f"{LABEL_DIR}/{os.path.basename(LABEL_FILES['1kg'])}", sep=r"\s+", dtype=str)
hg = pd.read_csv(f"{LABEL_DIR}/{os.path.basename(LABEL_FILES['hgdp'])}", sep="\t", dtype=str)
labels = pd.concat([
    pd.DataFrame({"s": kg["SampleID"], "pop": kg["Population"], "region": kg["Superpopulation"],
                  "project": "1000G", "founder": (kg["FatherID"] == "0") & (kg["MotherID"] == "0")}),
    pd.DataFrame({"s": hg["sample"], "pop": hg["population"], "region": hg["region"],
                  "project": "HGDP", "founder": True}),
    # HGDP genomes sequenced by the SGDP appear under their SGDP ID, the library suffix
    pd.DataFrame({"s": hg["library"].str.split(".", n=1).str[1], "pop": hg["population"],
                  "region": hg["region"], "project": "HGDP", "founder": True}),
]).dropna(subset=["s"]).drop_duplicates("s")

ref = pd.read_csv(TRAINING_PCA, sep="\t", usecols=["s", "scores"])
R = parse_pcs(ref.pop("scores"))
ref = ref.merge(labels, on="s", how="left", validate="one_to_one")
assert len(ref) == len(R)

unlabelled = ref.loc[ref["pop"].isna(), "s"]
print(f"{len(aou):,} AoU participants, {len(ref):,} reference samples, "
      f"{len(unlabelled):,} without a 1000G or HGDP label"
      + (f" (e.g. {unlabelled.head(3).tolist()})" if len(unlabelled) else ""))
print(ref.groupby(["project", "region", "pop"]).agg(n=("s", "size"), founders=("founder", "sum"))
      .to_string())
missing = [p for ps in GATES.values() for p in ps if p not in set(ref["pop"])]
assert not missing, f"not among the reference labels: {missing}"


def gate_rows(pops):
    """Founder reference samples of these populations, as a boolean mask over R."""
    return (ref["pop"].isin(pops) & ref["founder"].eq(True)).to_numpy()
```

## Cell 3 — gates

Mahalanobis distance to each group's centroid under the group's own covariance:
a multivariate normal fitted to its founders, so the gate follows the shape of
the cloud rather than a sphere. Saves the participant list of every gate at
every width.

- `growth`: rise in count from 1× to the widest gate; near 1 for a cluster, large for a continuum
- `spread`: dispersion of the selected participants relative to the reference

```python
def mahalanobis(Z, mu, C):
    d = Z[:, :N_PCS] - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, np.linalg.inv(C), d))


dist, d99, kept, centre = {}, {}, {}, {}
for g, pops in GATES.items():
    Rg = R[gate_rows(pops), :N_PCS]
    assert len(Rg) > 2 * N_PCS, f"{g}: too few founders for a {N_PCS}x{N_PCS} covariance"
    mu, C = Rg.mean(0), np.cov(Rg, rowvar=False)
    centre[g] = (mu, C)
    dist[g] = mahalanobis(X, mu, C)
    d99[g] = np.quantile(mahalanobis(Rg, mu, C), 0.99)
    for w in WIDTHS:
        kept[(g, w)] = dist[g] <= w * d99[g]
        aou.loc[kept[(g, w)], "person_id"].to_csv(
            f"{KEEP_DIR}/{g}_{w:g}x_keep_ids.txt", index=False, header=False)

rows = []
for g, pops in GATES.items():
    Rg = R[gate_rows(pops), :N_PCS]
    mu, C = centre[g]
    k1 = kept[(g, 1.0)]
    spread = (X[k1, :N_PCS].std(0, ddof=1) / np.sqrt(np.diag(C))).max() if k1.sum() > 1 else np.nan
    ref_d = mahalanobis(Rg, mu, C)
    rows.append({"gate": g, "ref_n": len(Rg),
                 "ref_in_narrow": (ref_d <= W_NARROW * d99[g]).mean(),
                 **{f"n_{w:g}x": int(kept[(g, w)].sum()) for w in WIDTHS},
                 "growth": kept[(g, W_WIDE)].sum() / max(k1.sum(), 1), "spread": spread})
T = pd.DataFrame(rows)
print(T.assign(**{c: T[c].map(show) for c in T if c.startswith("n_")}).round(2).to_string(index=False))
print(f"\nkeep lists -> {KEEP_DIR}")

mult = np.linspace(0.25, 3, 60)
fig, ax = plt.subplots(figsize=(8, 5))
for g in GATES:
    ds = np.sort(dist[g] / d99[g])
    ax.plot(mult, np.maximum(np.searchsorted(ds, mult, side="right"), 1), color=GCOL[g], label=g)
ax.axvline(1, color="grey", ls=":", lw=1)
ax.set_yscale("log")
ax.set_xlabel("gate width (× distance holding 99% of the reference)")
ax.set_ylabel("AoU participants inside")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/kept_vs_width.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4 — all gates in PC space

Gates at `W_PLOT` as dots, their boundaries as ellipses, reference founders as
crosses.

```python
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

W_PLOT = 1.0
rng = np.random.default_rng(0)


def gate_ellipse(ax, g, i, j, width, **kw):
    """Outline of the gate on PCs i, j: the shadow of the ellipsoid, which is set
    by the 2x2 marginal covariance, so a full metric tilts it."""
    mu, C = centre[g]
    r = width * d99[g]
    M = C[np.ix_([i, j], [i, j])]
    vals, vecs = np.linalg.eigh(M)
    angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
    ax.add_patch(Ellipse((mu[i], mu[j]), 2 * r * np.sqrt(vals[-1]), 2 * r * np.sqrt(vals[0]),
                         angle=angle, fill=False, **kw))


fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for ax, (i, j) in zip(axes, [(0, 1), (2, 3)]):
    ax.hist2d(X[:, i], X[:, j], bins=300, cmap="Greys", norm=LogNorm(), zorder=0)
    for g in GATES:
        r = gate_rows(GATES[g])
        ax.scatter(R[r, i], R[r, j], s=14, marker="x", color=GCOL[g], lw=0.8, alpha=0.6, zorder=1)
    for g in GATES:
        m = kept[(g, W_PLOT)]
        idx = rng.choice(np.flatnonzero(m), size=min(5_000, m.sum()), replace=False) if m.any() else []
        ax.scatter(X[idx, i], X[idx, j], s=3, alpha=0.5, color=GCOL[g], lw=0, rasterized=True, zorder=2)
        gate_ellipse(ax, g, i, j, W_PLOT, edgecolor=GCOL[g], lw=1.4, zorder=3)
    ax.set_xlabel(f"PC{i + 1}"); ax.set_ylabel(f"PC{j + 1}")
axes[0].legend(handles=[Line2D([], [], lw=0, marker="o", ms=6, color=GCOL[g], label=g) for g in GATES]
                       + [Line2D([], [], lw=0, marker="x", ms=6, color="0.4", label="reference founders")],
               fontsize=8)
fig.suptitle(f"AoU (grey), reference founders (crosses), participants inside each gate at "
             f"{W_PLOT:g}× (dots) and its boundary (ellipse)")
plt.tight_layout()
plt.savefig(f"{OUT}/gates_pcs_{W_PLOT:g}x.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4b — each gate

Dark dots: narrow gate. Light dots: added by the widest gate. Ellipses: narrow
(solid), 1× (dotted), widest (dashed).

```python
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D

ZOOM = True              # False: full PC range
PC_PAIRS = [(0, 1), (2, 3)]

for g, pops in GATES.items():
    in_narrow = kept[(g, W_NARROW)]
    in_wide = kept[(g, W_WIDE)] & ~in_narrow
    r = gate_rows(pops)
    pop_col = dict(zip(pops, plt.cm.Dark2.colors))
    light = 0.35 * np.array(to_rgb(GCOL[g])) + 0.65

    fig, axes = plt.subplots(1, len(PC_PAIRS), figsize=(15, 6.5))
    for ax, (i, j) in zip(axes, PC_PAIRS):
        if ZOOM:
            mu, C = centre[g]
            half = W_WIDE * d99[g] * np.sqrt(np.diag(C[np.ix_([i, j], [i, j])]))
            pts = np.r_[X[in_narrow | in_wide][:, [i, j]], R[r][:, [i, j]],
                        [mu[[i, j]] - half, mu[[i, j]] + half]]
            lo, hi = pts.min(0), pts.max(0)
            pad = 0.15 * np.maximum(hi - lo, 1e-9)
            lo, hi = lo - pad, hi + pad
        else:
            lo, hi = X[:, [i, j]].min(0), X[:, [i, j]].max(0)
        window = ((X[:, i] >= lo[0]) & (X[:, i] <= hi[0]) & (X[:, j] >= lo[1]) & (X[:, j] <= hi[1]))
        if window.any():
            ax.hist2d(X[window, i], X[window, j], bins=200, cmap="Greys", norm=LogNorm(),
                      range=[[lo[0], hi[0]], [lo[1], hi[1]]])
        for p in pops:
            m = gate_rows([p])
            ax.scatter(R[m, i], R[m, j], s=28, marker="x", color=pop_col[p], lw=1.2, zorder=1)
        for m, col, a in ((in_wide, light, 0.5), (in_narrow, GCOL[g], 0.6)):
            idx = np.flatnonzero(m)
            idx = rng.choice(idx, size=min(20_000, len(idx)), replace=False)
            ax.scatter(X[idx, i], X[idx, j], s=3, color=col, alpha=a, lw=0, rasterized=True, zorder=2)
        gate_ellipse(ax, g, i, j, W_NARROW, edgecolor=GCOL[g], lw=1.6, zorder=5)
        gate_ellipse(ax, g, i, j, 1.0, edgecolor=GCOL[g], lw=1.0, ls=":", zorder=5)
        gate_ellipse(ax, g, i, j, W_WIDE, edgecolor=GCOL[g], lw=1.2, ls="--", zorder=5)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
        ax.set_xlabel(f"PC{i + 1}"); ax.set_ylabel(f"PC{j + 1}")

    handles = ([Line2D([], [], lw=0, marker="o", ms=6, color=GCOL[g],
                       label=f"inside {W_NARROW:g}× ({show(int(in_narrow.sum()))})"),
                Line2D([], [], lw=0, marker="o", ms=6, color=light,
                       label=f"added up to {W_WIDE:g}× ({show(int(in_wide.sum()))})")]
               + [Line2D([], [], lw=0, marker="x", ms=7, mew=1.5, color=pop_col[p],
                         label=f"{p} ({int(gate_rows([p]).sum())})") for p in pops]
               + [Line2D([], [], color=GCOL[g], lw=1.6, label=f"gate at {W_NARROW:g}×"),
                  Line2D([], [], color=GCOL[g], lw=1.0, ls=":", label="gate at 1×"),
                  Line2D([], [], color=GCOL[g], lw=1.2, ls="--", label=f"gate at {W_WIDE:g}×")])
    axes[0].legend(handles=handles, fontsize=8, loc="best", framealpha=0.9)
    fig.suptitle(f"{g} gate", color=GCOL[g])
    plt.tight_layout()
    plt.savefig(f"{OUT}/gate_{g}_pcs{'_zoom' if ZOOM else ''}.png", dpi=130, bbox_inches="tight")
    plt.show()
```

## Cell 4c — selection vs reference

Offset of each reference population, and of the selected participants, from
the gate centroid, in reference SDs. `inside` checks the ellipses (should be 1).

```python
W_REPORT = 1.0
PCS = [f"PC{k + 1}" for k in range(N_PCS)]
for g, pops in GATES.items():
    mu, C = centre[g]
    sd = np.sqrt(np.diag(C))
    m = kept[(g, W_REPORT)]
    rows = {f"{p} (founders)": (R[gate_rows([p]), :N_PCS].mean(0) - mu) / sd for p in pops}
    if m.sum() > 20:
        rows[f"AoU inside {W_REPORT:g}×"] = (X[m, :N_PCS].mean(0) - mu) / sd
    inside = {}
    if m.any():
        for i, j in PC_PAIRS:
            y = X[m][:, [i, j]] - mu[[i, j]]
            Minv = np.linalg.inv(C[np.ix_([i, j], [i, j])])
            d2 = np.einsum("ij,jk,ik->i", y, Minv, y)
            inside[f"PC{i + 1}/PC{j + 1}"] = round(float((d2 <= (W_REPORT * d99[g]) ** 2).mean()), 3)
    print(f"\n{g}  n inside {W_REPORT:g}× = {show(int(m.sum()))}   inside: {inside}")
    print(pd.DataFrame(rows, index=PCS).T.round(2).to_string())
```

## Cell 4e — why the drawn ellipse is not the gate

The gate uses PCs 1–5; the ellipse is its outline on two of them. Participants
inside the ellipse are coloured by what the remaining PCs add to their distance,
given these two, in units of the gate limit. Above 1 they are excluded whatever
this panel shows.

```python
G = "NAT"                 # gate to inspect
I, J = 0, 1               # PCs on the axes
W_SHOW = W_WIDE

# Under the MVN, the squared distance splits exactly into what these two PCs
# show (the marginal term, which is what the drawn ellipse bounds) and what the
# other PCs add given them (the conditional term).
mu, C = centre[G]
d = X[:, :N_PCS] - mu
M = C[np.ix_([I, J], [I, J])]
d_plane = np.sqrt(np.einsum("ij,jk,ik->i", d[:, [I, J]], np.linalg.inv(M), d[:, [I, J]]))
d_other = np.sqrt(np.maximum(dist[G] ** 2 - d_plane ** 2, 0))
lim = W_SHOW * d99[G]
inside_ellipse = d_plane <= lim
passes = inside_ellipse & (dist[G] <= lim)

print(f"{G}: inside the {W_SHOW:g}× ellipse on PC{I + 1}/PC{J + 1}: {show(int(inside_ellipse.sum()))}; "
      f"of those, in the gate: {show(int(passes.sum()))}")
excluded = inside_ellipse & ~passes
if excluded.sum() > 20:
    sd = np.sqrt(np.diag(C))
    print("mean |PC - centroid| in reference SDs among those excluded:",
          dict(zip([f"PC{k + 1}" for k in range(N_PCS)],
                   np.abs(d[excluded] / sd).mean(0).round(2))))

fig, ax = plt.subplots(figsize=(8.5, 7.5))
ax.hist2d(X[:, I], X[:, J], bins=300, cmap="Greys", norm=LogNorm(), zorder=0)
sc = ax.scatter(X[inside_ellipse, I], X[inside_ellipse, J], c=d_other[inside_ellipse] / lim,
                s=3, cmap="RdYlBu_r", vmin=0, vmax=2, lw=0, rasterized=True, zorder=2)
for w, ls in ((W_NARROW, "-"), (1.0, ":"), (W_SHOW, "--")):
    gate_ellipse(ax, G, I, J, w, edgecolor=GCOL[G], lw=1.3, ls=ls, zorder=3)
pts = X[inside_ellipse][:, [I, J]]
pad = 0.1 * np.ptp(pts, axis=0)
ax.set_xlim(pts[:, 0].min() - pad[0], pts[:, 0].max() + pad[0])
ax.set_ylim(pts[:, 1].min() - pad[1], pts[:, 1].max() + pad[1])
ax.set_xlabel(f"PC{I + 1}"); ax.set_ylabel(f"PC{J + 1}")
fig.colorbar(sc, ax=ax, label="what the other PCs add, given these two (× gate limit); >1 excluded")
ax.set_title(f"{G}: participants inside the {W_SHOW:g}× ellipse on these two PCs")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_{G}_other_pcs.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4d — overlap between gates

Share of each row gate's participants also inside each column gate.

```python
W_OVERLAP = 1.0
names = list(GATES)
sets = {g: kept[(g, W_OVERLAP)] for g in names}
C = np.array([[int((sets[a] & sets[b]).sum()) for b in names] for a in names])
n = np.diag(C)
share = np.divide(C, n[:, None], out=np.full(C.shape, np.nan), where=n[:, None] > 0)
small = (C > 0) & (C <= 20)

counts = pd.DataFrame(C, index=names, columns=names)
counts.to_csv(f"{OUT}/gate_overlap_counts_{W_OVERLAP:g}x.tsv", sep="\t")
print(f"participants in both gates at {W_OVERLAP:g}× (row ∩ column)")
print(counts.apply(lambda c: c.map(show)).to_string())

fig, ax = plt.subplots(figsize=(7.5, 6.5))
im = ax.imshow(np.where(small, np.nan, share), vmin=0, vmax=1, cmap="Blues")
for a in range(len(names)):
    for b in range(len(names)):
        if a == b:
            txt = f"n={show(int(n[a]))}"
        elif small[a, b]:
            txt = "≤20"
        elif np.isfinite(share[a, b]):
            txt = f"{share[a, b]:.2f}"
        else:
            txt = ""
        ax.text(b, a, txt, ha="center", va="center", fontsize=8,
                color="white" if np.isfinite(share[a, b]) and share[a, b] > 0.6 else "black")
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names)
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_color(GCOL[lab.get_text()])
ax.set_xlabel("also inside gate …")
ax.set_ylabel("participants of gate …")
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="share of row gate")
ax.set_title(f"overlap between gates at {W_OVERLAP:g}×")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_overlap_{W_OVERLAP:g}x.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5 — load kinship

AoU's kinship table and the participant set of every gate.

```python
if not os.path.isfile(REL_PATH):
    rel_dir = os.path.dirname(REL_PATH)
    print(f"no {os.path.basename(REL_PATH)}; set REL_PATH and REL_COLS from what {rel_dir} holds:")
    for f in sorted(os.listdir(rel_dir)):
        p = f"{rel_dir}/{f}"
        cols = (pd.read_csv(p, sep="\t", nrows=0).columns.tolist()
                if os.path.isfile(p) and f.endswith((".tsv", ".txt", ".tsv.gz")) else "")
        print(f"  {f}  {os.path.getsize(p) / 1e6:.1f} MB  {cols}")
else:
    header = pd.read_csv(REL_PATH, sep="\t", nrows=0).columns.tolist()
    assert all(c in header for c in REL_COLS), f"set REL_COLS from {header}"
    i_col, j_col, k_col = REL_COLS
    rel = pd.read_csv(REL_PATH, sep="\t", usecols=list(REL_COLS), dtype={i_col: str, j_col: str})
    kin = rel[k_col].to_numpy()
    print(f"{show(len(rel))} pairs, kinship {kin.min():.3f} to {kin.max():.3f}")

    members = {(g, f"{w:g}×"): set(aou.loc[kept[(g, w)], "person_id"])
               for g in GATES for w in WIDTHS}
```

## Cell 5b — related pairs by kinship class

Pairs inside each gate: second degree, first degree, duplicates or MZ twins.

```python
KIN_CLASSES = ["< 0.177", "0.177–0.354", "≥ 0.354"]
kin_class = np.select([kin < KIN_CUTS[0], kin < KIN_CUTS[1]], KIN_CLASSES[:2], KIN_CLASSES[2])

rows = []
for (g, w), ids in members.items():
    both = (rel[i_col].isin(ids) & rel[j_col].isin(ids)).to_numpy()
    by_class = pd.Series(kin_class[both]).value_counts().reindex(KIN_CLASSES, fill_value=0)
    rows.append({"group": g, "width": w, "participants": len(ids), **by_class.to_dict()})
K = pd.DataFrame(rows).set_index(["group", "width"])
K.to_csv(f"{OUT}/kinship_classes.tsv", sep="\t")

display(K.apply(lambda c: c.map(show)).set_axis(
    pd.MultiIndex.from_tuples([("", "participants")] + [("pairs by kinship", c) for c in KIN_CLASSES]),
    axis=1))
```

## Cell 5c — related pairs, plotted

```python
KIN_TITLES = {"< 0.177": "kinship 0.1–0.177\n(2nd degree)",
              "0.177–0.354": "kinship 0.177–0.354\n(1st degree)",
              "≥ 0.354": "kinship ≥ 0.354\n(duplicates, MZ twins)"}
groups = list(GATES)
widths = [f"{w:g}×" for w in WIDTHS]
alphas = dict(zip(widths, np.linspace(1.0, 0.35, len(widths))))
bw = 0.8 / len(widths)

fig, axes = plt.subplots(1, len(KIN_CLASSES), figsize=(17, 5), sharey=True)
for ax, c in zip(axes, KIN_CLASSES):
    for x, g in enumerate(groups):
        for k, w in enumerate(widths):
            v = int(K.loc[(g, w), c])
            xpos = x - 0.4 + bw * (k + 0.5)
            if v > 20:
                ax.bar(xpos, v, width=bw, color=GCOL[g], alpha=alphas[w], edgecolor="0.25", lw=0.4)
            elif v > 0:
                ax.text(xpos, 22, "≤20", rotation=90, ha="center", va="bottom", fontsize=7, color="0.35")
    ax.set_yscale("log")
    ax.set_ylim(bottom=20)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    for lab in ax.get_xticklabels():
        lab.set_color(GCOL[lab.get_text()])
    ax.set_title(KIN_TITLES[c], fontsize=10)
    ax.grid(axis="y", color="0.9", lw=0.6)
    ax.set_axisbelow(True)
axes[0].set_ylabel("pairs with both members in the gate")
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="0.3", alpha=alphas[w], label=f"gate at {w}")
                    for w in widths],
           loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=len(widths), frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUT}/kinship_classes.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5d — kinship histograms

All pairs, and pairs inside each gate. Dotted lines mark degree cutoffs.

```python
W_REL = 1.0
CUTS = {"2nd": KIN_CUTS[0] / 2, "1st": KIN_CUTS[0], "dup/MZ": KIN_CUTS[1]}
edges = np.linspace(kin.min(), kin.max(), 90)


def counts_shown(v):
    c = np.histogram(v, edges)[0].astype(float)
    c[c <= 20] = np.nan
    return c


fig, axes = plt.subplots(1, 2, figsize=(16, 5.5), sharex=True)
axes[0].stairs(counts_shown(kin), edges, fill=True, color="0.6")
axes[0].set_title(f"all pairs in the table ({show(len(kin))})")

for g in GATES:
    ids = members[(g, f"{W_REL:g}×")]
    both = (rel[i_col].isin(ids) & rel[j_col].isin(ids)).to_numpy()
    if both.sum() > 20:
        axes[1].stairs(counts_shown(kin[both]), edges, color=GCOL[g], lw=1.4,
                       label=f"{g} ({show(int(both.sum()))})")
axes[1].set_title(f"pairs with both members inside the gate at {W_REL:g}×")
axes[1].legend(fontsize=8)

for ax in axes:
    ax.set_yscale("log")
    ax.set_xlabel("kinship coefficient (PC-Relate)")
    for lab, x in CUTS.items():
        if edges[0] <= x <= edges[-1]:
            ax.axvline(x, color="0.4", ls=":", lw=1)
            ax.text(x, 1.0, f" {lab}", transform=ax.get_xaxis_transform(),
                    ha="left", va="top", fontsize=8, color="0.3")
axes[0].set_ylabel("pairs")
plt.tight_layout()
plt.savefig(f"{OUT}/kinship_histograms.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 6 — summary

Tables hold raw counts and stay in the workbench.

```python
T.to_csv(f"{OUT}/gate_summary.tsv", sep="\t", index=False)
print(sorted(os.listdir(OUT)))
```

## Cell 7 — copy this notebook to the bucket

Save the notebook first, then run this. Set `NOTEBOOK` to its path.

```python
NOTEBOOK = os.path.expanduser("~/notebooks/ancestry_survey.ipynb")
NB_DIR = f"{OUT}/notebooks"
os.makedirs(NB_DIR, exist_ok=True)
assert os.path.isfile(NOTEBOOK), NOTEBOOK

dest = f"{NB_DIR}/{os.path.basename(NOTEBOOK)}"
shutil.copy(NOTEBOOK, dest)
ok = os.path.getsize(dest) == os.path.getsize(NOTEBOOK)
print(f"{'OK  ' if ok else 'SIZE MISMATCH'} {dest}")
```
