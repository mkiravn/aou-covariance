"""Ancestry gates. Round 1 is 01_round1_gate's code; round 2 uses reweighted
MCD and a chi-square radius (see docs/pipeline_architecture.md for why not OGK)."""

from __future__ import annotations

import ast

import numpy as np
import pandas as pd


# ---- inputs ------------------------------------------------------------------

def parse_pcs(col: pd.Series, n: int = 16) -> np.ndarray:
    arr = np.vstack(col.map(ast.literal_eval).to_numpy())
    assert arr.shape[1] == n, arr.shape
    return arr


def read_1kg_labels(path: str) -> pd.DataFrame:
    """1000G ped-population file -> sample, pop, region, founder."""
    kg = pd.read_csv(path, sep=r"\s+", dtype=str)
    return pd.DataFrame({
        "s": kg["SampleID"], "pop": kg["Population"], "region": kg["Superpopulation"],
        "project": "1000G",
        "founder": (kg["FatherID"] == "0") & (kg["MotherID"] == "0"),
    })


# ---- distances and gates -------------------------------------------------------

def mahalanobis(Z: np.ndarray, mu: np.ndarray, C: np.ndarray) -> np.ndarray:
    d = Z[:, : len(mu)] - mu
    return np.sqrt(np.einsum("ij,jk,ik->i", d, np.linalg.inv(C), d))


def fit_quantile_gate(R: np.ndarray, width: float) -> dict:
    """Round 1: MVN on reference founders; the radius holds `width` of them."""
    assert len(R) > 2 * R.shape[1], f"too few founders ({len(R)}) for {R.shape[1]} PCs"
    mu, C = R.mean(0), np.cov(R, rowvar=False)
    ref_d = mahalanobis(R, mu, C)
    return {"mu": mu, "C": C, "radius": float(np.quantile(ref_d, width)), "ref_d": ref_d}


def fit_robust_gate(A: np.ndarray, q: float, seed: int = 0) -> dict:
    """Round 2: reweighted MCD on the anchor; radius = sqrt(chi2_K quantile q)."""
    from scipy.stats import chi2
    from sklearn.covariance import MinCovDet

    K = A.shape[1]
    assert len(A) > 5 * K, f"anchor of {len(A)} is too small for {K} PCs"
    mcd = MinCovDet(random_state=seed).fit(A)
    mu, C = mcd.location_, mcd.covariance_
    return {"mu": mu, "C": C, "radius": float(np.sqrt(chi2.ppf(q, K))),
            "ref_d": mahalanobis(A, mu, C),
            "anchor_outliers": int((~mcd.support_).sum())}


def coverage_curve(dist: np.ndarray, ref_d: np.ndarray,
                   cover=np.linspace(0.05, 0.999, 60)) -> pd.DataFrame:
    """Participants inside the gate as reference coverage grows."""
    radii = np.quantile(ref_d, cover)
    return pd.DataFrame({"coverage": cover, "radius": radii,
                         "kept": [int((dist <= r).sum()) for r in radii]})


# ---- plotting -----------------------------------------------------------------

def gate_ellipse(ax, mu, C, i: int, j: int, r: float, **kw) -> None:
    """Shadow of the ellipsoid on PCs i, j (the 2x2 marginal covariance)."""
    from matplotlib.patches import Ellipse

    M = C[np.ix_([i, j], [i, j])]
    vals, vecs = np.linalg.eigh(M)
    ax.add_patch(Ellipse((mu[i], mu[j]),
                         2 * r * np.sqrt(vals[-1]), 2 * r * np.sqrt(vals[0]),
                         angle=np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1])),
                         fill=False, **kw))
