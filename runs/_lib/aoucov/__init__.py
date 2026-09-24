"""Shared pipeline code for the AoU phenotypic-covariance runs.

Usage in a workbench notebook:

    import sys; sys.path.insert(0, os.path.expanduser("~/aou-covariance/runs/_lib"))
    from aoucov import Run, env, plink, refs, gate, plots, checks, provenance

    run = Run("eur_r2")

Nothing here reads or writes participant data to anywhere but the workspace
bucket. Outputs never land in the repo working tree.
"""

from .config import (  # noqa: F401
    CDR_VERSION,
    PCA_DISK,
    PCA_MACHINE,
    PCA_MEM_MB,
    QC_DISK,
    QC_MACHINE,
    QC_MEM_MB,
    Run,
    commit_sha,
)
from . import batch, checks, env, gate, plink, plots, provenance, refs  # noqa: F401

__all__ = [
    "Run", "commit_sha", "CDR_VERSION",
    "QC_MACHINE", "QC_MEM_MB", "QC_DISK",
    "PCA_MACHINE", "PCA_MEM_MB", "PCA_DISK",
    "batch", "checks", "env", "gate", "plink", "plots", "provenance", "refs",
]
