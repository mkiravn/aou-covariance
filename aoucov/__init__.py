"""Run bookkeeping for the AoU phenotypic-covariance pipeline.

Four small pieces, usable from a laptop, a Workbench VM, or inside a dsub task:

  config      a run is a committed JSON/YAML file plus a run id
  store       uniform exists/list/read/write over local paths and gs:// URIs
  manifest    what a stage should produce, and what is actually there
  provenance  a .meta.json sidecar beside every output saying how it was made

Nothing here submits jobs or knows what a GRM is. It answers "what is done,
what is missing, and which code produced it".
"""

from .config import RunConfig, load_run
from .config import _repo_root as repo_root
from .store import get_store, LocalStore, GcsStore
from .manifest import targets, missing, present, stage_status, run_status
from .provenance import collect_meta, write_meta, read_meta, record

__all__ = [
    "RunConfig", "load_run", "repo_root",
    "get_store", "LocalStore", "GcsStore",
    "targets", "missing", "present", "stage_status", "run_status",
    "collect_meta", "write_meta", "read_meta", "record",
]
