"""Standard diagnostic figures.

Each returns the saved path. Participants are always drawn as a hexbin or a
rasterised subsample -- never one marker per person, which produces a 40 MB PNG
and hides the density anyway.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from .plink import pc_names, read_eigenval
from .refs import ANCHOR_POPS, EUR_POPS, POP_COLORS


def scree(pca_prefix, out_png, n_pcs=20, title=""):
    _, pct = read_eigenval(pca_prefix)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(1, len(pct) + 1), pct, marker="o", ms=4, color="steelblue")
    ax.set_xlabel("PC")
    ax.set_ylabel("% variance explained")
    ax.set_title(title or "scree")
    plt.tight_layout()
    plt.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.show()
    return out_png


def anchor_separation(kg_scores, part_scores, n_pcs=20,
                      anchor=ANCHOR_POPS, eur=EUR_POPS):
    """How far the anchor sits from the other Europeans on each PC.

    In participant SDs, so it is comparable across PCs. Pick K where this is
    large -- a PC on which the anchor is indistinguishable from the rest of
    Europe contributes noise to the gate, not selectivity.
    """
    eur_kg = kg_scores[kg_scores["pop"].isin(eur)]
    a = eur_kg["pop"].isin(anchor)
    o = eur_kg["pop"].isin(set(eur) - set(anchor))
    rows = []
    for p in pc_names(n_pcs):
        d = abs(eur_kg.loc[a, p].mean() - eur_kg.loc[o, p].mean())
        rows.append(d / part_scores[p].std())
    return np.array(rows)


def pc_pairs(part, kg, out_png, pairs=((0, 1), (2, 3)), gate=None,
             pops=EUR_POPS, title="", subsample=50_000, mask=None):
    """Participants as density, reference populations as crosses, gate as ellipse."""
    PC = pc_names(max(max(p) for p in pairs) + 1)
    rng = np.random.default_rng(0)
    idx = np.flatnonzero(mask) if mask is not None else np.arange(len(part))
    shown = rng.choice(idx, size=min(subsample, len(idx)), replace=False)
    cols = dict(zip(pops, plt.cm.tab10.colors))

    fig, axes = plt.subplots(1, len(pairs), figsize=(7.5 * len(pairs), 6.5))
    for ax, (i, j) in zip(np.atleast_1d(axes), pairs):
        a, b = PC[i], PC[j]
        ax.scatter(part[a], part[b], s=1, alpha=0.1, color="0.8",
                   rasterized=True, zorder=0)
        ax.scatter(part[a].to_numpy()[shown], part[b].to_numpy()[shown], s=2,
                   alpha=0.3, color="tab:blue", rasterized=True, zorder=1)
        for pop in pops:
            s = kg[kg["pop"] == pop]
            ax.scatter(s[a], s[b], s=30 if pop in ANCHOR_POPS else 16, marker="x",
                       color=cols[pop], label=pop, zorder=2)
        if gate is not None and i < len(gate.pcs) and j < len(gate.pcs):
            ellipse(ax, gate.mu, gate.cov, i, j, gate.threshold,
                    edgecolor="tab:blue", lw=1.4, zorder=3)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
    np.atleast_1d(axes)[0].legend(fontsize=8, markerscale=1.4)
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.show()
    return out_png


def ellipse(ax, centre, cov, i, j, radius, **kw):
    """The gate's cross-section in the (i, j) plane."""
    M = cov[np.ix_([i, j], [i, j])]
    vals, vecs = np.linalg.eigh(M)
    ax.add_patch(Ellipse(
        (centre[i], centre[j]),
        2 * radius * np.sqrt(vals[-1]), 2 * radius * np.sqrt(vals[0]),
        angle=np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1])),
        fill=False, **kw))


def loadings_by_position(pca_prefix, out_png, n_show=4, title=""):
    """Loading^2 against genomic position, one panel per PC.

    A PC carried by one region is an inversion, an LD block or a mapping
    artefact -- not population structure. Read this before trusting any PC.
    """
    import pandas as pd

    L = pd.read_csv(f"{pca_prefix}.eigenvec.allele", sep=r"\s+")
    idc = "#ID" if "#ID" in L.columns else "ID"
    L[["CHROM", "POS"]] = L[idc].str.split(":", n=2, expand=True).iloc[:, :2]
    L["CHROM"] = L["CHROM"].str.replace("chr", "", regex=False)
    autosomes = [str(c) for c in range(1, 23)]
    L = L[L["CHROM"].isin(autosomes)].copy()
    L["CHROM"] = pd.Categorical(L["CHROM"], categories=autosomes, ordered=True)
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

    fig, axes = plt.subplots(n_show, 1, figsize=(14, 2.4 * n_show), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, pc_names(n_show)):
        for k, c in enumerate(centres):
            s = L[L["CHROM"] == c]
            ax.scatter(s["CUM"], s[p] ** 2, s=2, alpha=0.5,
                       color=["tab:blue", "tab:orange"][k % 2], rasterized=True)
        ax.set_ylabel(f"{p} loading²")
    axes[-1].set_xticks(list(centres.values()))
    axes[-1].set_xticklabels(list(centres), fontsize=7)
    axes[-1].set_xlabel("chromosome")
    fig.suptitle(title or "PCA loadings")
    plt.tight_layout()
    plt.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.show()
    return L


def subpop_panel(part_scores, kg_scores, out_png, n_pairs=4, n_pcs=20, title=""):
    """PC pairs 1-8 with the 1000G European subpopulations overlaid."""
    PC = pc_names(n_pcs)
    pairs = [(PC[i], PC[i + 1]) for i in range(0, min(2 * n_pairs, n_pcs - 1), 2)]
    eur_kg = kg_scores[kg_scores["super_pop"] == "EUR"]

    fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 4.8))
    for ax, (a, b) in zip(np.atleast_1d(axes), pairs):
        ax.hexbin(part_scores[a], part_scores[b], gridsize=70, cmap="Greys",
                  mincnt=1, bins="log", linewidths=0, zorder=1)
        for pop, col in POP_COLORS.items():
            sub = eur_kg[eur_kg["pop"] == pop]
            ax.scatter(sub[a], sub[b], s=12, color=col, label=pop,
                       alpha=0.85, linewidths=0, zorder=3)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
    np.atleast_1d(axes)[0].legend(fontsize=7, markerscale=1.5)
    plt.suptitle(title, y=1.01)
    plt.tight_layout()
    plt.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.show()
    return out_png
