"""Diagnostics that gate whether a PCA is usable.

These were copies in four notebooks. One implementation, called where relevant.
"""

import os

import numpy as np
import pandas as pd

from .plink import pc_names

DEFAULT_META_FACTORS = ["data_partner_id", "site_id", "sequencing_center"]


def ld_peaks(loadings, n_pcs=10, quantile=0.9995, flank_kb=250):
    """Find single-region PCs and return exclusion intervals for them.

    Flags loading outliers, widens each to a +/-flank window, and merges. If
    this returns anything, add it to the LD exclusion list and refit -- the PC
    is tracking one locus, not ancestry.

    `loadings` is the frame returned by `plots.loadings_by_position`.
    """
    L = loadings
    hit = np.zeros(len(L), dtype=bool)
    for p in pc_names(min(n_pcs, sum(c.startswith("PC") for c in L.columns))):
        v = L[p].to_numpy() ** 2
        hit |= v > np.quantile(v, quantile)
    print(f"{hit.sum()} loading outliers across PC1-{n_pcs}")
    if not hit.any():
        print("no peaks; PCs look clean")
        return []

    H = L.loc[hit, ["CHROM", "POS"]].sort_values(["CHROM", "POS"])
    regions, k = [], 0
    flank = flank_kb * 1000
    for chrom, sub in H.groupby("CHROM", sort=False, observed=True):
        lo = hi = None
        for pos in sub["POS"]:
            if lo is None or pos > hi + flank:
                if lo is not None:
                    k += 1
                    regions.append(f"chr{chrom} {max(lo - flank, 1)} {hi + flank} peak{k}")
                lo = pos
            hi = pos
        k += 1
        regions.append(f"chr{chrom} {max(lo - flank, 1)} {hi + flank} peak{k}")
    print(f"\n{len(regions)} peak regions -- add to the LD exclusion list and refit:")
    print("\n".join(regions))
    return regions


def batch_effect(part_scores, meta_path, factors=None, n_pcs=10):
    """Variance in each PC explained by batch / site metadata.

    A PC with high R^2 on collection site is measuring the assay, not the
    person; carrying it as a covariate is fine, interpreting it is not.
    """
    factors = factors or DEFAULT_META_FACTORS
    if not meta_path or not os.path.isfile(meta_path):
        print("no metadata table available -- skipping batch check")
        return None
    m = pd.read_csv(meta_path, sep="\t", dtype=str).rename(
        columns={"research_id": "person_id"})
    j = part_scores.merge(m, on="person_id", how="inner")
    PC = pc_names(n_pcs)
    rows = {f: [float(1 - j.groupby(f)[p].transform("var").mean() / j[p].var())
                for p in PC]
            for f in factors if f in j}
    if not rows:
        print(f"none of {factors} present in {meta_path}")
        return None
    out = pd.DataFrame(rows, index=PC).round(3)
    print(f"R² of PC1-{n_pcs} on batch factors:")
    print(out.to_string())
    return out
