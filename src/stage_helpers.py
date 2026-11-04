"""Panel re-checks and plots that never show a group of 1-20 participants."""

from __future__ import annotations

import os
import subprocess
import tempfile

import numpy as np
import pandas as pd

AMBIGUOUS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def verify_panel(rep, prefix: str, plink2: str, maf: float, geno: float, hwe: float,
                 label: str = "the sample set", ambiguous_ok: bool = False) -> pd.DataFrame:
    """Re-check the written panel, rather than trusting the command that made it."""
    bim = pd.read_csv(prefix + ".bim", sep="\t", header=None,
                      names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={"chr": str})
    n_amb = sum((a, b) in AMBIGUOUS for a, b in zip(bim["a1"], bim["a2"]))
    if not ambiguous_ok:
        rep.check("no strand-ambiguous SNPs (A/T, C/G)", n_amb == 0, f"{n_amb} found")
    acgt = list("ACGT")
    rep.check("biallelic SNPs only", bim["a1"].isin(acgt).all() and bim["a2"].isin(acgt).all())
    rep.check("no duplicate variant IDs", not bim["id"].duplicated().any())
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([plink2, "--bfile", prefix, "--freq", "--missing", "variant-only",
                        "--hardy", "--silent", "--out", f"{tmp}/v"], check=True)
        frq = pd.read_csv(f"{tmp}/v.afreq", sep="\t")
        vmiss = pd.read_csv(f"{tmp}/v.vmiss", sep="\t")
        hardy = pd.read_csv(f"{tmp}/v.hardy", sep="\t")
    m = np.minimum(frq["ALT_FREQS"], 1 - frq["ALT_FREQS"])
    rep.check(f"MAF >= {maf} within {label}", m.min() >= maf, f"min {m.min():.4f}")
    rep.check(f"variant missingness <= {geno}", vmiss["F_MISS"].max() <= geno,
              f"max {vmiss['F_MISS'].max():.4f}")
    rep.check(f"HWE p >= {hwe:g}", hardy["P"].min() >= hwe, f"min {hardy['P'].min():.2e}")
    return bim


def hexbin_safe(ax, x, y, **kw):
    """Participant density with no cell holding fewer than 21 people shown."""
    kw = {"gridsize": 120, "bins": "log", "cmap": "Greys", **kw, "mincnt": 21}
    return ax.hexbin(x, y, **kw)


def hist_safe(ax, x, bins=60, **kw):
    """Histogram with bins of 1-20 participants left blank."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    lo, hi = np.quantile(x, [0.005, 0.995])
    counts, edges = np.histogram(x, bins=np.linspace(lo, hi, bins + 1))
    shown = np.where(counts > 20, counts, np.nan)
    ax.stairs(np.nan_to_num(shown), edges, fill=True, alpha=0.6, **kw)
    return edges


def peak_enrichment(loadings: pd.DataFrame, pcs, window_mb: float = 10) -> pd.Series:
    """Max over windows of (share of squared loading / share of variants); ~1 = spread out."""
    w = loadings["chr"].astype(str) + ":" + (loadings["pos"] // (window_mb * 1e6)).astype(int).astype(str)
    share_var = w.value_counts(normalize=True)
    out = {}
    for pc in pcs:
        sq = loadings[pc] ** 2
        share = sq.groupby(w).sum() / sq.sum()
        out[pc] = float((share / share_var.reindex(share.index)).max())
    return pd.Series(out)
