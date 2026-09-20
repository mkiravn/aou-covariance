"""Relatedness-regression estimators on grm_class_tool merge output.

Each estimator regresses the phenotype cross-product y_i*y_j on GRM relatedness
a_ij over a range of relatedness (table in notes.md). The merged output holds
bin sums and pair counts rather than pairs, so every fit is pair-weighted least
squares on bin means -- identical to pair-level OLS with each pair's a_ij taken
as its bin midpoint.

Standard errors are delete-block jackknife: every estimator is refitted on each
of the merge's delete-block replicates (the .jk.tsv), using the same blocks the
accumulators were built with.

Class sets:
  noPO    other + FS    primary -- parent-offspring pairs removed
  pooled  other + FS + PO   what an unclassified run would give
  PO      PO only       for the FS-range estimators, as a contrast

Slope ladder, all slopes through the origin; U = a < 0.02, T = 0.02-0.05,
R = 0.05-0.7:
  h2_OneSlope          one slope over U+T+R
  h2_Rel               separate slope over T+R; its U slope is h2_Unrel
  h2_OneSlopeOffsets   one slope over U+R, plus an offset per degree band
  h2_RelOffsets        separate slope over R, plus an offset per degree band
Offsets are identified only above 0.05, so T is left out of the offset models.
With an offset in every band, the related slope of h2_RelOffsets comes only
from variation in a within each degree. U and R share no bins, so the U slope
of the two-slope models is exactly h2_Unrel. Differences are computed inside
each jackknife replicate, so their SEs carry the covariance of the two terms.
"""
import numpy as np
import pandas as pd

CLS = ["other", "FS", "PO"]
NOPO = ["other", "FS"]
POOLED = ["other", "FS", "PO"]

U_HI, T_HI, R_HI = 0.02, 0.05, 0.7
# KING degree cutoffs (kinship 0.0442/0.0884/0.177/0.354, doubled) rounded to bin edges
DEG_BANDS = {"deg4": (0.05, 0.09), "deg3": (0.09, 0.18),
             "deg2": (0.18, 0.36), "deg1": (0.36, R_HI)}


def load_grid(bins_path):
    """Bin midpoints, in bin_index order (bin_index - 1)."""
    return np.loadtxt(bins_path).mean(axis=1)


def load_arrays(full_path, jk_path, nbins, nblocks):
    """Bin sums S and pair counts N, shape (1 + nblocks, len(CLS), nbins).
    Replicate 0 is the full data; replicate b+1 omits jackknife block b."""
    f = pd.read_csv(full_path, sep="\t")
    j = pd.read_csv(jk_path, sep="\t")
    S = np.zeros((1 + nblocks, len(CLS), nbins))
    N = np.zeros_like(S)
    for ci, c in enumerate(CLS):
        fc = f[f["class"] == c]
        k = fc["bin_index"].to_numpy() - 1
        S[0, ci, k] = fc["full_sum"].to_numpy()
        N[0, ci, k] = fc["full_n"].to_numpy()
        jc = j[j["class"] == c]
        b = jc["block"].to_numpy() + 1
        k = jc["bin_index"].to_numpy() - 1
        S[b, ci, k] = jc["jk_sum"].to_numpy()
        N[b, ci, k] = jc["jk_n"].to_numpy()
    return S, N


def _ind(a, lo, hi):
    return ((a >= lo) & (a < hi)).astype(float)


def _wls(X, m, n):
    if len(m) <= X.shape[1]:
        return None
    w = np.sqrt(n)
    Xw = X * w[:, None]
    if np.linalg.matrix_rank(Xw) < X.shape[1]:
        return None
    beta, *_ = np.linalg.lstsq(Xw, m * w, rcond=None)
    return beta


def estimates(S, N, r, mid):
    """All estimators for replicate r, as {(estimator, classes): value}."""

    def bins(classes, *ranges):
        ci = [CLS.index(c) for c in classes]
        s, n = S[r, ci].sum(0), N[r, ci].sum(0)
        inside = np.zeros(len(mid), dtype=bool)
        for lo, hi in ranges:
            inside |= (mid >= lo) & (mid < hi)
        k = (n > 0) & inside
        return mid[k], s[k] / n[k], n[k]

    out = {}

    def fit(label, classes, ranges, design, readouts):
        a, m, n = bins(classes, *ranges)
        beta = _wls(design(a), m, n) if len(a) else None
        for name, readout in readouts.items():
            out[(name, label)] = readout(beta) if beta is not None else np.nan

    origin = lambda a: a[:, None]
    step = lambda a: np.c_[a, _ind(a, .4, .6), _ind(a, .2, .4)]
    two_slopes = lambda a: np.c_[a * (a < U_HI), a * (a >= U_HI)]
    offsets = lambda a: np.column_stack([_ind(a, lo, hi) for lo, hi in DEG_BANDS.values()])
    unrel = (-np.inf, U_HI)

    # SNP heritability from unrelated pairs; no intercept
    fit("noPO", NOPO, [unrel], origin, {"h2_Unrel": lambda b: b[0]})
    # same fit with an intercept, as a diagnostic for residual structure
    fit("noPO", NOPO, [unrel], lambda a: np.c_[a, np.ones_like(a)],
        {"h2_UnrelInt": lambda b: b[0], "Unrel.intercept": lambda b: b[1]})

    for label, classes in (("noPO", NOPO), ("pooled", POOLED), ("PO", ["PO"])):
        # classic sib regression -- slope conflates h2 and shared environment
        fit(label, classes, [(0.4, 0.6)], origin, {"h2_FS": lambda b: b[0]})
        # shared environment among first-degree pairs
        fit(label, classes, [(0.4, 0.6)], lambda a: np.c_[a, np.ones_like(a)],
            {"b2_FS": lambda b: 2 * b[1], "b2_FS.slope": lambda b: b[0]})

    band_readouts = lambda model, first: {f"{model}.{d}": (lambda b, i=first + j: b[i])
                                          for j, d in enumerate(DEG_BANDS)}
    for label, classes in (("noPO", NOPO), ("pooled", POOLED)):
        # quadratic absorbs non-linearity (Wainschtein et al. 2025)
        fit(label, classes, [(0.05, 0.7)], lambda a: np.c_[a, a ** 2],
            {"h2_PedW25": lambda b: b[0], "PedW25.quad": lambda b: b[1]})
        # intercept offsets at pedigree-class transitions
        fit(label, classes, [(0.05, 0.7)],
            lambda a: np.c_[a, _ind(a, .1, .2), _ind(a, .2, .4), _ind(a, .4, .6)],
            {"h2_Pedf": lambda b: b[0]})
        # separate FS (0.4-0.6) and HS (0.2-0.4) intercepts, common slope
        fit(label, classes, [(0.2, 0.6)], step,
            {"b2_step": lambda b: 4 * (b[1] - b[2]),
             "b2_step_ratio": lambda b: b[1] / b[2] if b[2] != 0 else np.nan})

        fit(label, classes, [(-np.inf, R_HI)], origin, {"h2_OneSlope": lambda b: b[0]})
        fit(label, classes, [(-np.inf, R_HI)], two_slopes,
            {"h2_Rel": lambda b: b[1], "diff_Rel-Unrel": lambda b: b[1] - b[0]})
        fit(label, classes, [unrel, (T_HI, R_HI)], lambda a: np.c_[a, offsets(a)],
            {"h2_OneSlopeOffsets": lambda b: b[0], **band_readouts("OneSlopeOffsets", 1)})
        fit(label, classes, [unrel, (T_HI, R_HI)], lambda a: np.c_[two_slopes(a), offsets(a)],
            {"h2_RelOffsets": lambda b: b[1], "diff_RelOffsets-Unrel": lambda b: b[1] - b[0],
             **band_readouts("RelOffsets", 2)})
        out[("diff_Rel-RelOffsets", label)] = out[("h2_Rel", label)] - out[("h2_RelOffsets", label)]

    # mean cross-product above the additive prediction h2_Unrel * a_ij
    h2u = out[("h2_Unrel", "noPO")]
    for cls in ("FS", "PO"):
        a, m, n = bins([cls], (0.4, 0.6))
        out[("excess", cls)] = (np.sum(n * (m - h2u * a)) / np.sum(n)
                                if len(a) and np.isfinite(h2u) else np.nan)
        out[("n_pairs", cls)] = float(np.sum(n))
    return out


def jackknife(S, N, mid, nblocks):
    """Point estimates on the full data with delete-block jackknife SEs.
    SE is NaN unless every block's replicate produced a finite estimate."""
    full = estimates(S, N, 0, mid)
    reps = [estimates(S, N, r, mid) for r in range(1, nblocks + 1)]
    rows = []
    for key, est in full.items():
        th = np.array([rp[key] for rp in reps], dtype=float)
        th = th[np.isfinite(th)]
        B = len(th)
        se = (np.sqrt((B - 1) / B * np.sum((th - th.mean()) ** 2))
              if B == nblocks and np.isfinite(est) else np.nan)
        rows.append({"estimator": key[0], "classes": key[1], "est": est, "se": se})
    return pd.DataFrame(rows)
