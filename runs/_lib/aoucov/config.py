"""Paths, machine specs, and provenance for one run.

Everything that used to be a 40-line config cell repeated in every notebook.
A notebook says `run = Run("eur_r2")` and gets every path from it.
"""

import os
import subprocess
from dataclasses import dataclass, field

# --- workbench constants -----------------------------------------------------
# These describe the workspace, not any one run. They change when the workspace
# is cloned or the CDR is bumped -- see notes in runs/README.md.

PROJECT_ID = "wb-swift-sprout-7231"
REGION = "us-central1"
SERVICE_ACCOUNT = "pet-27799165194323faf22e2@wb-swift-sprout-7231.iam.gserviceaccount.com"
NETWORK = f"projects/{PROJECT_ID}/global/networks/network"
SUBNETWORK = f"projects/{PROJECT_ID}/regions/{REGION}/subnetworks/subnetwork"
CLOUD_SDK_IMAGE = "gcr.io/google.com/cloudsdktool/cloud-sdk:581.0.0-slim"

CDR_VERSION = "v9"
WS = os.path.expanduser("~/workspace/Data from All of Us Controlled Tier /shared-env-pilot")
WS_GS = "gs://cloned-shared-env-pilot-wb-swift-sprout-7231"

# Default machines for the two heavy Batch steps. Override per call site.
QC_MACHINE, QC_MEM_MB, QC_DISK = "n1-highmem-16", 96_000, 500
PCA_MACHINE, PCA_MEM_MB, PCA_DISK = "n1-highmem-32", 200_000, 500


@dataclass
class Run:
    """One sample set. Bucket paths are the source of truth; `local` is scratch.

    Attribute pairs ending `_gs` are the gs:// form of the mounted path, for
    handing to dsub (Batch workers do not see the fuse mount).
    """

    sample_set: str
    cdr: str = CDR_VERSION

    b: str = field(init=False)
    b_gs: str = field(init=False)
    root: str = field(init=False)
    root_gs: str = field(init=False)
    logs_gs: str = field(init=False)

    def __post_init__(self):
        self.b = f"{WS}/phenotypic_covariance_{self.cdr}"
        self.b_gs = f"{WS_GS}/phenotypic_covariance_{self.cdr}"
        self.root = f"{self.b}/{self.sample_set}"
        self.root_gs = f"{self.b_gs}/{self.sample_set}"
        self.logs_gs = f"{self.root_gs}/dsub_logs"

    # --- output locations ----------------------------------------------------

    def out(self, *parts, make=True):
        """A directory under this run's bucket root, created on request."""
        p = os.path.join(self.root, *parts)
        if make:
            os.makedirs(p, exist_ok=True)
        return p

    def out_gs(self, *parts):
        return "/".join([self.root_gs, *parts])

    def scratch(self, tag, make=True):
        """Ephemeral local disk. Never the source of truth for anything."""
        p = os.path.expanduser(f"~/scratch_{self.sample_set}_{tag}")
        if make:
            os.makedirs(p, exist_ok=True)
        return p

    # --- shared inputs -------------------------------------------------------

    @property
    def panel_gs(self):
        return f"{self.b_gs}/01_ancestry_filtering/unified_panel/unified_panel_{self.cdr}"

    @property
    def plink2_gs(self):
        return f"{self.b_gs}/01_ancestry_filtering/unified_panel/bin/plink2"

    # --- Batch ---------------------------------------------------------------

    @property
    def batch_args(self):
        """The dsub flags that are identical for every job in this workspace."""
        return dict(
            project=PROJECT_ID,
            region=REGION,
            logging=self.logs_gs,
            service_account=SERVICE_ACCOUNT,
            network=NETWORK,
            subnetwork=SUBNETWORK,
            image=CLOUD_SDK_IMAGE,
        )

    def job_name(self, step):
        # dsub job names reject underscores.
        return f"{step}-{self.sample_set}".replace("_", "-")


def commit_sha():
    """The commit this code was pulled from, or 'unknown' outside a checkout.

    Stamped into every provenance file so a bucket artefact can be traced back
    to the exact code that produced it.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        r = subprocess.run(
            ["git", "-C", here, "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        )
        dirty = subprocess.run(
            ["git", "-C", here, "status", "--porcelain"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return r.stdout.strip() + ("-dirty" if dirty else "")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"
