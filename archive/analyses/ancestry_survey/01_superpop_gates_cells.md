# Ancestry survey

For each reference group of 1000 Genomes or HGDP populations:

1. fit a multivariate normal to its founders in AoU's ancestry PCs
2. gate AoU participants by Mahalanobis distance, at several coverage levels, and save the lists
3. plot each gate in PC space
4. measure overlap between gates, and assign participants where a gate is unambiguous
5. count related pairs inside each gate, by kinship class

A width is the **fraction of the reference population inside the gate**: 0.8 is
the ellipsoid holding 80% of that group's own founders. Any number of widths
can be listed; every cell follows `WIDTHS`.

Population labels come from the projects' own sample files. Kinship comes from
AoU's `samples_relatedness.tsv` (Hail PC-Relate, pairs above 0.1 only).
Counts of 1–20 are shown as `≤20`.

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
REL_PATH       = f"{AUX}/relatedness/samples_relatedness.tsv"
REL_COLS       = ("i.s", "j.s", "kin")     # participant 1, participant 2, kinship coefficient

WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
OUT = f"{WS}/phenotypic_covariance_v9/analyses/ancestry_filtering"
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
WIDTHS = [0.8, 0.9, 0.99]                    # fraction of the reference population inside the gate
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
pct = lambda w: f"{w * 100:g}%"

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

AoU and reference PCs, with the original population labels. Only founders
define the gates.

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

An MVN per group; each width's radius is the quantile of the founders' own
distances, so the gate contains that fraction of them by construction.

```python
def mahalanobis(Z, mu, C):
    d = Z[:, :N_PCS] - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, np.linalg.inv(C), d))


dist, radius, kept, centre = {}, {}, {}, {}
for g, pops in GATES.items():
    Rg = R[gate_rows(pops), :N_PCS]
    assert len(Rg) > 2 * N_PCS, f"{g}: too few founders for a {N_PCS}x{N_PCS} covariance"
    mu, C = Rg.mean(0), np.cov(Rg, rowvar=False)
    centre[g] = (mu, C)
    dist[g] = mahalanobis(X, mu, C)
    ref_d = mahalanobis(Rg, mu, C)
    for w in WIDTHS:
        radius[(g, w)] = float(np.quantile(ref_d, w))
        kept[(g, w)] = dist[g] <= radius[(g, w)]
        aou.loc[kept[(g, w)], "person_id"].to_csv(
            f"{KEEP_DIR}/{g}_{w * 100:g}pct_keep_ids.txt", index=False, header=False)

rows = []
for g, pops in GATES.items():
    Rg = R[gate_rows(pops), :N_PCS]
    mu, C = centre[g]
    k = kept[(g, WIDTHS[0])]
    rows.append({"gate": g, "ref_n": len(Rg),
                 **{f"n_{w * 100:g}pct": int(kept[(g, w)].sum()) for w in WIDTHS},
                 "growth": kept[(g, WIDTHS[-1])].sum() / max(kept[(g, WIDTHS[0])].sum(), 1),
                 "spread": (X[k, :N_PCS].std(0, ddof=1) / np.sqrt(np.diag(C))).max()
                 if k.sum() > 1 else np.nan})
T = pd.DataFrame(rows)
print(T.assign(**{c: T[c].map(show) for c in T if c.startswith("n_")}).round(2).to_string(index=False))
print(f"\nkeep lists -> {KEEP_DIR}")

cover = np.linspace(0.05, 0.999, 60)
fig, ax = plt.subplots(figsize=(8, 5))
for g in GATES:
    mu, C = centre[g]
    radii = np.quantile(mahalanobis(R[gate_rows(GATES[g]), :N_PCS], mu, C), cover)
    ax.plot(cover * 100, [np.sum(dist[g] <= r) for r in radii], color=GCOL[g], label=g)
for w in WIDTHS:
    ax.axvline(w * 100, color="0.85", lw=0.8, zorder=0)
ax.set_yscale("log")
ax.set_xlabel("reference population inside the gate (%)")
ax.set_ylabel("AoU participants inside")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/kept_vs_reference_coverage.png", dpi=130, bbox_inches="tight")
plt.show()
```

`growth` is the widest gate's count over the narrowest: near 1 for a cluster,
large for a continuum. `spread` compares the selected participants' dispersion
with the reference's, at the narrowest width.

## Cell 3b — which reference populations land in each gate

Rows are reference populations, columns are gates, cells are the percentage of
that population's founders falling inside. A gate holds its own populations by
construction, at the width's coverage; the informative part is everything off
that diagonal, which is where gates reach into groups they don't represent.
Reference samples are public, so counts are unsuppressed.

```python
W_CHECK = WIDTHS[0]
EXTRA_POPS = ["ASW", "ACB"]          # in neither gate, shown for context

pops_shown = [p for ps in GATES.values() for p in ps] + [
    p for p in EXTRA_POPS if p in set(ref["pop"])]
rows_ = []
for p in pops_shown:
    m = gate_rows([p])
    if not m.any():
        continue
    row = {"population": p, "founders": int(m.sum()),
           "in GATES": next((g for g, ps in GATES.items() if p in ps), "—")}
    for g in GATES:
        mu, C = centre[g]
        row[g] = 100 * (mahalanobis(R[m], mu, C) <= radius[(g, W_CHECK)]).mean()
    rows_.append(row)
P = pd.DataFrame(rows_).set_index("population")
P.to_csv(f"{OUT}/reference_in_gates_{W_CHECK * 100:g}pct.tsv", sep="\t")
print(P.round(1).to_string())

fig, ax = plt.subplots(figsize=(1.1 * len(GATES) + 3, 0.28 * len(P) + 2))
im = ax.imshow(P[list(GATES)].to_numpy(), vmin=0, vmax=100, cmap="Blues", aspect="auto")
for i in range(len(P)):
    for j, g in enumerate(GATES):
        v = P.iloc[i][g]
        if v >= 1:
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7,
                    color="white" if v > 60 else "black")
ax.set_xticks(range(len(GATES)))
ax.set_xticklabels(list(GATES))
for lab in ax.get_xticklabels():
    lab.set_color(GCOL[lab.get_text()])
ax.set_yticks(range(len(P)))
ax.set_yticklabels([f"{p}  ({P.loc[p, 'in GATES']})" for p in P.index], fontsize=8)
ax.tick_params(length=0)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="% of the population's founders inside")
ax.set_title(f"reference populations inside each gate at {pct(W_CHECK)}")
plt.tight_layout()
plt.savefig(f"{OUT}/reference_in_gates_{W_CHECK * 100:g}pct.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4 — all gates in PC space

```python
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

W_PLOT = WIDTHS[0]
LS = ["-", "--", ":", "-."]
rng = np.random.default_rng(0)


def gate_ellipse(ax, g, i, j, w, **kw):
    """The ellipsoid's shadow on PCs i, j: the 2x2 marginal covariance, so it
    tilts with the group."""
    mu, C = centre[g]
    M = C[np.ix_([i, j], [i, j])]
    vals, vecs = np.linalg.eigh(M)
    r = radius[(g, w)]
    ax.add_patch(Ellipse((mu[i], mu[j]), 2 * r * np.sqrt(vals[-1]), 2 * r * np.sqrt(vals[0]),
                         angle=np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1])), fill=False, **kw))


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
fig.suptitle(f"AoU (grey), reference founders (crosses), participants inside each gate at {pct(W_PLOT)} "
             "(dots) and its boundary (ellipse)")
plt.tight_layout()
plt.savefig(f"{OUT}/gates_pcs_{W_PLOT * 100:g}pct.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 4b — each gate

One figure per gate, with every width drawn. Dots are shaded by the narrowest
gate that contains them.

```python
from matplotlib.colors import to_rgb

ZOOM = True              # False: full PC range
PC_PAIRS = [(0, 1), (2, 3)]


def shade(base, k, n):
    """Full colour for the narrowest gate, lighter for each wider one."""
    t = 0.25 + 0.75 * (1 - k / max(n - 1, 1))
    return tuple(t * np.array(to_rgb(base)) + (1 - t))


for g, pops in GATES.items():
    pop_col = dict(zip(pops, plt.cm.Dark2.colors))
    bands, prev = [], np.zeros(len(X), dtype=bool)
    for k, w in enumerate(WIDTHS):
        m = kept[(g, w)] & ~prev
        bands.append((m, shade(GCOL[g], k, len(WIDTHS)), f"{pct(w)} ({show(int(m.sum()))})"))
        prev = prev | kept[(g, w)]

    fig, axes = plt.subplots(1, len(PC_PAIRS), figsize=(15, 6.5))
    for ax, (i, j) in zip(axes, PC_PAIRS):
        mu, C = centre[g]
        if ZOOM:
            half = radius[(g, WIDTHS[-1])] * np.sqrt(np.diag(C[np.ix_([i, j], [i, j])]))
            pts = np.r_[X[prev][:, [i, j]], R[gate_rows(pops)][:, [i, j]],
                        [mu[[i, j]] - half, mu[[i, j]] + half]]
            lo, hi = pts.min(0), pts.max(0)
            pad = 0.15 * np.maximum(hi - lo, 1e-9)
            lo, hi = lo - pad, hi + pad
        else:
            lo, hi = X[:, [i, j]].min(0), X[:, [i, j]].max(0)
        window = (X[:, i] >= lo[0]) & (X[:, i] <= hi[0]) & (X[:, j] >= lo[1]) & (X[:, j] <= hi[1])
        if window.any():
            ax.hist2d(X[window, i], X[window, j], bins=200, cmap="Greys", norm=LogNorm(),
                      range=[[lo[0], hi[0]], [lo[1], hi[1]]])
        for p in pops:
            m = gate_rows([p])
            ax.scatter(R[m, i], R[m, j], s=28, marker="x", color=pop_col[p], lw=1.2, zorder=1)
        for m, col, _ in reversed(bands):
            idx = np.flatnonzero(m)
            if len(idx):
                idx = rng.choice(idx, size=min(20_000, len(idx)), replace=False)
                ax.scatter(X[idx, i], X[idx, j], s=3, color=col, alpha=0.6, lw=0,
                           rasterized=True, zorder=2)
        for k, w in enumerate(WIDTHS):
            gate_ellipse(ax, g, i, j, w, edgecolor=GCOL[g], lw=1.3, ls=LS[k % len(LS)], zorder=5)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
        ax.set_xlabel(f"PC{i + 1}"); ax.set_ylabel(f"PC{j + 1}")

    handles = ([Line2D([], [], lw=0, marker="o", ms=6, color=col, label=lab)
                for _, col, lab in bands]
               + [Line2D([], [], lw=0, marker="x", ms=7, mew=1.5, color=pop_col[p],
                         label=f"{p} ({int(gate_rows([p]).sum())})") for p in pops]
               + [Line2D([], [], color=GCOL[g], lw=1.3, ls=LS[k % len(LS)], label=f"gate at {pct(w)}")
                  for k, w in enumerate(WIDTHS)])
    axes[0].legend(handles=handles, fontsize=8, loc="best", framealpha=0.9)
    fig.suptitle(f"{g} gate", color=GCOL[g])
    plt.tight_layout()
    plt.savefig(f"{OUT}/gate_{g}_pcs{'_zoom' if ZOOM else ''}.png", dpi=130, bbox_inches="tight")
    plt.show()
```

## Cell 4c — selection vs reference

Offsets from the gate centroid in reference SDs. `inside` should be 1.

```python
W_REPORT = WIDTHS[0]
PCS = [f"PC{k + 1}" for k in range(N_PCS)]
for g, pops in GATES.items():
    mu, C = centre[g]
    sd = np.sqrt(np.diag(C))
    m = kept[(g, W_REPORT)]
    rows_ = {f"{p} (founders)": (R[gate_rows([p]), :N_PCS].mean(0) - mu) / sd for p in pops}
    if m.sum() > 20:
        rows_[f"AoU inside {pct(W_REPORT)}"] = (X[m, :N_PCS].mean(0) - mu) / sd
    inside = {}
    if m.any():
        for i, j in PC_PAIRS:
            y = X[m][:, [i, j]] - mu[[i, j]]
            d2 = np.einsum("ij,jk,ik->i", y, np.linalg.inv(C[np.ix_([i, j], [i, j])]), y)
            inside[f"PC{i + 1}/PC{j + 1}"] = round(float((d2 <= radius[(g, W_REPORT)] ** 2).mean()), 3)
    print(f"\n{g}  n inside {pct(W_REPORT)} = {show(int(m.sum()))}   inside: {inside}")
    print(pd.DataFrame(rows_, index=PCS).T.round(2).to_string())
```

## Cell 4d — why the ellipse is not the gate

The gate uses PCs 1–5; the ellipse is its outline on two of them. Colour is
what the remaining PCs add, given these two: above 1 the participant is out.

```python
G = "NAT"
I, J = 0, 1
W_SHOW = WIDTHS[-1]

mu, C = centre[G]
d = X[:, :N_PCS] - mu
M = C[np.ix_([I, J], [I, J])]
d_plane = np.sqrt(np.einsum("ij,jk,ik->i", d[:, [I, J]], np.linalg.inv(M), d[:, [I, J]]))
d_other = np.sqrt(np.maximum(dist[G] ** 2 - d_plane ** 2, 0))
lim = radius[(G, W_SHOW)]
in_ellipse = d_plane <= lim
passes = in_ellipse & (dist[G] <= lim)

print(f"{G}: inside the {pct(W_SHOW)} ellipse on PC{I + 1}/PC{J + 1}: {show(int(in_ellipse.sum()))}; "
      f"of those, in the gate: {show(int(passes.sum()))}")
excluded = in_ellipse & ~passes
if excluded.sum() > 20:
    print("mean |PC - centroid| in reference SDs among those excluded:",
          dict(zip(PCS, np.abs(d[excluded] / np.sqrt(np.diag(C))).mean(0).round(2))))

fig, ax = plt.subplots(figsize=(8.5, 7.5))
ax.hist2d(X[:, I], X[:, J], bins=300, cmap="Greys", norm=LogNorm(), zorder=0)
sc = ax.scatter(X[in_ellipse, I], X[in_ellipse, J], c=d_other[in_ellipse] / lim,
                s=3, cmap="RdYlBu_r", vmin=0, vmax=2, lw=0, rasterized=True, zorder=2)
for k, w in enumerate(WIDTHS):
    gate_ellipse(ax, G, I, J, w, edgecolor=GCOL[G], lw=1.3, ls=LS[k % len(LS)], zorder=3)
pts = X[in_ellipse][:, [I, J]]
pad = 0.1 * np.ptp(pts, axis=0)
ax.set_xlim(pts[:, 0].min() - pad[0], pts[:, 0].max() + pad[0])
ax.set_ylim(pts[:, 1].min() - pad[1], pts[:, 1].max() + pad[1])
ax.set_xlabel(f"PC{I + 1}"); ax.set_ylabel(f"PC{J + 1}")
fig.colorbar(sc, ax=ax, label="what the other PCs add, given these two (× gate radius); >1 excluded")
ax.set_title(f"{G}: participants inside the {pct(W_SHOW)} ellipse on these two PCs")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_{G}_other_pcs.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5 — assignment at one width

A participant is assigned only when exactly one gate contains them.

```python
W = WIDTHS[0]
LABELS = list(GATES) + ["multiple", "unassigned"]

masks = np.column_stack([kept[(g, W)] for g in GATES])
n_in = masks.sum(1)
assigned = np.where(n_in > 1, "multiple", "unassigned").astype(object)
for k, g in enumerate(GATES):
    assigned[masks[:, k] & (n_in == 1)] = g
counts = pd.Series(assigned).value_counts().reindex(LABELS, fill_value=0)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.barh(range(len(counts)), np.maximum(counts.values, 1),
        color=[GCOL.get(l, "0.6") for l in counts.index])
ax.set_yticks(range(len(counts)))
ax.set_yticklabels(counts.index)
ax.invert_yaxis()
ax.set_xscale("log")
ax.set_xlabel("AoU participants")
for i, n in enumerate(counts.values):
    ax.text(max(n, 1) * 1.15, i, f"{show(int(n))}  ({n / len(assigned) * 100:.1f}%)",
            va="center", fontsize=8, color="0.3")
ax.set_title(f"assignment at {pct(W)}: a gate is claimed only when no other gate also contains them")
plt.tight_layout()
plt.savefig(f"{OUT}/assignment_{W * 100:g}pct.png", dpi=130, bbox_inches="tight")
plt.show()

print(pd.DataFrame({"assignment": counts.index, "n": counts.map(show).values,
                    "pct": (counts.values / len(assigned) * 100).round(2)}).to_string(index=False))
```

## Cell 5b — assignment across widths

How the composition moves as the gates widen: more people covered, more of them
claimed by two gates at once.

```python
comp, per_gate = [], []
for w in WIDTHS:
    masks = np.column_stack([kept[(g, w)] for g in GATES])
    n_in = masks.sum(1)
    a = np.where(n_in > 1, "multiple", "unassigned").astype(object)
    for k, g in enumerate(GATES):
        a[masks[:, k] & (n_in == 1)] = g
    comp.append({"width": pct(w), "assigned": int((n_in == 1).sum()),
                 "multiple": int((n_in > 1).sum()), "unassigned": int((n_in == 0).sum())})
    per_gate.append(pd.Series(a).value_counts().reindex(list(GATES), fill_value=0).rename(pct(w)))

A = pd.DataFrame(comp).set_index("width")
G_ = pd.concat(per_gate, axis=1)
A.to_csv(f"{OUT}/assignment_by_width.tsv", sep="\t")
G_.to_csv(f"{OUT}/assignment_by_width_per_gate.tsv", sep="\t")
print(A.apply(lambda c: c.map(show)).to_string())
print("\nuniquely assigned to each gate")
print(G_.apply(lambda c: c.map(show)).to_string())

fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
bottom = np.zeros(len(A))
for col, c in (("assigned", "#2a78d6"), ("multiple", "#eb6834"), ("unassigned", "0.8")):
    axes[0].bar(A.index, A[col], bottom=bottom, color=c, label=col)
    bottom += A[col].to_numpy()
axes[0].set_ylabel("AoU participants")
axes[0].set_xlabel("reference population inside the gate")
axes[0].legend(fontsize=8)
axes[0].set_title("composition")

for g in GATES:
    axes[1].plot(G_.columns, np.maximum(G_.loc[g], 1), marker="o", ms=4, color=GCOL[g], label=g)
axes[1].set_yscale("log")
axes[1].set_ylabel("uniquely assigned")
axes[1].set_xlabel("reference population inside the gate")
axes[1].legend(fontsize=8, ncol=2)
axes[1].set_title("per gate")
plt.tight_layout()
plt.savefig(f"{OUT}/assignment_by_width.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 5c — overlap between gates, at every width

Row gate A, column gate B: the share of A's participants also inside B.

```python
names = list(GATES)
fig, axes = plt.subplots(1, len(WIDTHS), figsize=(5.2 * len(WIDTHS), 5), squeeze=False)
for ax, w in zip(axes[0], WIDTHS):
    sets = {g: kept[(g, w)] for g in names}
    Cn = np.array([[int((sets[a] & sets[b]).sum()) for b in names] for a in names])
    n = np.diag(Cn)
    share = np.divide(Cn, n[:, None], out=np.full(Cn.shape, np.nan), where=n[:, None] > 0)
    small = (Cn > 0) & (Cn <= 20)
    pd.DataFrame(Cn, index=names, columns=names).to_csv(
        f"{OUT}/gate_overlap_counts_{w * 100:g}pct.tsv", sep="\t")

    im = ax.imshow(np.where(small, np.nan, share), vmin=0, vmax=1, cmap="Blues")
    for a in range(len(names)):
        for b in range(len(names)):
            txt = (f"n={show(int(n[a]))}" if a == b else
                   "≤20" if small[a, b] else
                   f"{share[a, b]:.2f}" if np.isfinite(share[a, b]) else "")
            ax.text(b, a, txt, ha="center", va="center", fontsize=7,
                    color="white" if np.isfinite(share[a, b]) and share[a, b] > 0.6 else "black")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=8)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_color(GCOL[lab.get_text()])
    ax.set_title(f"{pct(w)} of the reference inside")
axes[0][0].set_ylabel("participants of gate …")
fig.suptitle("overlap between gates: share of the row gate also inside the column gate")
plt.tight_layout()
plt.savefig(f"{OUT}/gate_overlap.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 6 — load kinship

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

    members = {(g, w): set(aou.loc[kept[(g, w)], "person_id"]) for g in GATES for w in WIDTHS}
```

## Cell 7 — related pairs by kinship class

Second degree, first degree, duplicates or MZ twins.

A big gate holds many pairs simply by being big: both members have to be inside,
so counts grow with the square of the gate's share of the cohort. `expected_1st`
is that null — the first-degree pairs a gate this size would hold if relatives
were spread evenly across the cohort — and `ratio` is observed over expected. A
gate that is merely large sits near 1; a gate enriched for relatives sits above.

```python
KIN_CLASSES = ["< 0.177", "0.177–0.354", "≥ 0.354"]
kin_class = np.select([kin < KIN_CUTS[0], kin < KIN_CUTS[1]], KIN_CLASSES[:2], KIN_CLASSES[2])
n_first = int((kin_class == KIN_CLASSES[1]).sum())

rows_ = []
for (g, w), ids in members.items():
    both = (rel[i_col].isin(ids) & rel[j_col].isin(ids)).to_numpy()
    by_class = pd.Series(kin_class[both]).value_counts().reindex(KIN_CLASSES, fill_value=0)
    share = len(ids) / len(aou)
    expected = n_first * share ** 2
    rows_.append({"gate": g, "width": pct(w), "participants": len(ids),
                  **by_class.to_dict(), "share_of_cohort": share,
                  "expected_1st": expected,
                  "ratio": by_class[KIN_CLASSES[1]] / expected if expected > 0 else np.nan})
K = pd.DataFrame(rows_).set_index(["gate", "width"])
K.to_csv(f"{OUT}/kinship_classes.tsv", sep="\t")

out = K.copy()
for c in KIN_CLASSES + ["participants"]:
    out[c] = out[c].map(show)
out["expected_1st"] = out["expected_1st"].round(0).astype(int).map(show)
out["share_of_cohort"] = (out["share_of_cohort"] * 100).round(1)
out["ratio"] = out["ratio"].round(2)
display(out)
```

**Check the ID spaces match** before reading anything into these counts: the
relatedness table's IDs must be the same `person_id` the PCs use. If they were
not, every count would be zero rather than wrong, so a table of zeros means a
join problem, not an absence of relatives.

## Cell 8 — related pairs, plotted

```python
KIN_TITLES = {"< 0.177": "kinship 0.1–0.177\n(2nd degree)",
              "0.177–0.354": "kinship 0.177–0.354\n(1st degree)",
              "≥ 0.354": "kinship ≥ 0.354\n(duplicates, MZ twins)"}
alphas = dict(zip(WIDTHS, np.linspace(1.0, 0.35, len(WIDTHS))))
bw = 0.8 / len(WIDTHS)

fig, axes = plt.subplots(1, len(KIN_CLASSES), figsize=(17, 5), sharey=True)
for ax, c in zip(axes, KIN_CLASSES):
    for x, g in enumerate(GATES):
        for k, w in enumerate(WIDTHS):
            v = int(K.loc[(g, pct(w)), c])
            xpos = x - 0.4 + bw * (k + 0.5)
            if v > 20:
                ax.bar(xpos, v, width=bw, color=GCOL[g], alpha=alphas[w], edgecolor="0.25", lw=0.4)
            elif v > 0:
                ax.text(xpos, 22, "≤20", rotation=90, ha="center", va="bottom", fontsize=7, color="0.35")
    ax.set_yscale("log")
    ax.set_ylim(bottom=20)
    ax.set_xticks(range(len(GATES)))
    ax.set_xticklabels(list(GATES))
    for lab in ax.get_xticklabels():
        lab.set_color(GCOL[lab.get_text()])
    ax.set_title(KIN_TITLES[c], fontsize=10)
    ax.grid(axis="y", color="0.9", lw=0.6)
    ax.set_axisbelow(True)
axes[0].set_ylabel("pairs with both members in the gate")
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="0.3", alpha=alphas[w], label=pct(w))
                    for w in WIDTHS],
           loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=len(WIDTHS), frameon=False, fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUT}/kinship_classes.png", dpi=130, bbox_inches="tight")
plt.show()
```

## Cell 9 — kinship histograms

All pairs, and pairs inside each gate at `W_REL`.

```python
W_REL = WIDTHS[0]
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
    ids = members[(g, W_REL)]
    both = (rel[i_col].isin(ids) & rel[j_col].isin(ids)).to_numpy()
    if both.sum() > 20:
        axes[1].stairs(counts_shown(kin[both]), edges, color=GCOL[g], lw=1.4,
                       label=f"{g} ({show(int(both.sum()))})")
axes[1].set_title(f"pairs with both members inside the gate at {pct(W_REL)}")
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

## Cell 10 — summary

Tables hold raw counts and stay in the workbench.

```python
T.to_csv(f"{OUT}/gate_summary.tsv", sep="\t", index=False)
pd.DataFrame([{"gate": g, "width": pct(w), "radius": radius[(g, w)]}
              for g in GATES for w in WIDTHS]).to_csv(f"{OUT}/gate_radii.tsv", sep="\t", index=False)
print(sorted(os.listdir(OUT)))
```

## Cell 11 — copy this notebook to the bucket

Save the notebook first, then run this.

```python
NOTEBOOK = os.path.expanduser("~/notebooks/ancestry_survey.ipynb")
NB_DIR = f"{OUT}/notebooks"
os.makedirs(NB_DIR, exist_ok=True)
assert os.path.isfile(NOTEBOOK), NOTEBOOK
dest = f"{NB_DIR}/{os.path.basename(NOTEBOOK)}"
shutil.copy(NOTEBOOK, dest)
print(f"{'OK  ' if os.path.getsize(dest) == os.path.getsize(NOTEBOOK) else 'SIZE MISMATCH'} {dest}")
```
