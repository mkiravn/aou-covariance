"""The ancestry gate: a multivariate normal fitted to an anchor population.

Mean and covariance come from the projected anchor samples; participants are
kept by Mahalanobis distance under that MVN, with the radius set by a coverage
quantile of the anchor's own distances. The anchor supplies the orientation and
the relative weighting of PCs. Projection shrinks the anchor's scores toward the
origin, but that only rescales the distances, which the coverage quantile
absorbs.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .plink import pc_names


@dataclass
class Gate:
    mu: np.ndarray
    cov: np.ndarray
    threshold: float
    pcs: list
    coverage: float
    n_anchor: int

    def distance(self, frame):
        d = frame[self.pcs].to_numpy() - self.mu
        return np.sqrt(np.einsum("ij,jk,ik->i", d, np.linalg.inv(self.cov), d))

    def keep(self, frame):
        return self.distance(frame) <= self.threshold


def fit(anchor, k_pcs, coverage=0.90):
    """Fit the gate to `anchor` (the projected reference samples)."""
    pcs = pc_names(k_pcs)
    A = anchor[pcs].to_numpy()
    assert len(A) > 2 * k_pcs, f"too few anchor samples for a {k_pcs}x{k_pcs} covariance"
    mu, cov = A.mean(0), np.cov(A, rowvar=False)

    g = Gate(mu, cov, threshold=np.inf, pcs=pcs, coverage=coverage, n_anchor=len(A))
    g.threshold = float(np.quantile(g.distance(anchor), coverage))
    print(f"anchor: {len(A)} samples; SD per PC {np.sqrt(np.diag(cov)).round(4)}")
    print(f"threshold: {coverage:.0%} of anchor distances = {g.threshold:.4f}")
    return g


def self_gate(participants, k_pcs, coverage=0.90):
    """Comparison gate centred on the participants themselves (what eur_D2 did).

    No anchor: the cloud defines its own centre, so the gate cannot be anchored
    to an external population and drifts with whoever happens to be in the set.
    """
    pcs = pc_names(k_pcs)
    A = participants[pcs].to_numpy()
    mu, cov = A.mean(0), np.cov(A, rowvar=False)
    g = Gate(mu, cov, threshold=np.inf, pcs=pcs, coverage=coverage, n_anchor=len(A))
    g.threshold = float(np.quantile(g.distance(participants), coverage))
    return g


def overlap(sets, out_tsv=None):
    """Jaccard overlap between named keep-sets."""
    names = list(sets)
    J = pd.DataFrame(
        [[len(sets[a] & sets[b]) / len(sets[a] | sets[b]) for b in names] for a in names],
        index=names, columns=names,
    )
    if out_tsv:
        J.to_csv(out_tsv, sep="\t")
    return J
