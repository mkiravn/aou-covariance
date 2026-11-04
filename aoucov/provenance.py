"""A .meta.json sidecar beside every output.

The point: six months later, "which code, which config and which machine
produced this shard" is a `gcloud storage cat`, not archaeology. It also makes
a whole class of silent bug visible -- a stage whose outputs were produced by
two different commits is exactly how an N-mismatch sneaks between rungs, and
`aoucov status` flags it because the sidecars disagree.
"""

from __future__ import annotations

import datetime as _dt
import getpass
import json
import os
import platform
import socket
import subprocess
import time
from contextlib import contextmanager
from functools import lru_cache
from typing import Any

SCHEMA = 1


def meta_uri(target: str) -> str:
    return f"{target}.meta.json"


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _git(args: list[str], repo: str) -> str | None:
    try:
        proc = subprocess.run(
            ["git", "-C", repo] + args, capture_output=True, text=True,
            check=False, timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


@lru_cache(maxsize=4)
def code_version(repo: str | None = None) -> dict[str, Any]:
    """How the running code is identified.

    Inside a dsub container there is usually no .git, so the job script is
    expected to export AOUCOV_GIT_SHA when it installs the pinned package.
    """
    env_sha = os.environ.get("AOUCOV_GIT_SHA")
    repo = repo or os.environ.get("AOUCOV_REPO") or os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
    )
    sha = env_sha or _git(["rev-parse", "HEAD"], repo)
    dirty = None
    branch = None
    if not env_sha and sha:
        status = _git(["status", "--porcelain"], repo)
        dirty = bool(status) if status is not None else None
        branch = _git(["rev-parse", "--abbrev-ref", "HEAD"], repo)
    return {
        "git_sha": sha,
        "git_sha_source": "env" if env_sha else "git",
        "dirty": dirty,
        "branch": branch,
        "package_version": os.environ.get("AOUCOV_VERSION"),
    }


@lru_cache(maxsize=1)
def machine_info() -> dict[str, Any]:
    info = {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpus": os.cpu_count(),
        "machine_type": os.environ.get("AOUCOV_MACHINE_TYPE"),
        "image": os.environ.get("AOUCOV_IMAGE"),
        "backend": os.environ.get("AOUCOV_BACKEND"),
    }
    try:
        info["user"] = getpass.getuser()
    except Exception:
        info["user"] = None
    # dsub / Batch leave these around; harmless and empty elsewhere
    for key in ("DSUB_JOB_ID", "DSUB_TASK_ID", "BATCH_TASK_INDEX", "BATCH_JOB_UID"):
        if os.environ.get(key):
            info.setdefault("job", {})[key] = os.environ[key]
    if info["machine_type"] is None:
        info["machine_type"] = _gce_machine_type()
    return info


def _gce_machine_type() -> str | None:
    """Best effort, short timeout; returns None off GCE."""
    try:
        import urllib.request

        req = urllib.request.Request(
            "http://metadata.google.internal/computeMetadata/v1/instance/machine-type",
            headers={"Metadata-Flavor": "Google"},
        )
        with urllib.request.urlopen(req, timeout=0.5) as resp:
            return resp.read().decode().rsplit("/", 1)[-1]
    except Exception:
        return None


def collect_meta(
    cfg,
    stage: str,
    target: str,
    item: dict | None = None,
    tools: dict[str, str] | None = None,
    extra: dict | None = None,
) -> dict:
    return {
        "schema": SCHEMA,
        "run_id": cfg.run_id,
        "config_hash": cfg.config_hash(),
        "config_source": cfg.source_path,
        "stage": stage,
        "item": item or {},
        "target": target,
        "code": code_version(),
        "machine": machine_info(),
        "tools": tools or {},
        "written_at": _now(),
        **(extra or {}),
    }


def write_meta(store, target: str, meta: dict) -> str:
    uri = meta_uri(target)
    store.write_text(uri, json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return uri


def read_meta(store, target: str) -> dict | None:
    uri = meta_uri(target)
    try:
        if not store.exists(uri):
            return None
        return json.loads(store.read_text(uri))
    except Exception:
        return None


@contextmanager
def record(cfg, stage: str, target: str, store=None, item: dict | None = None,
           tools: dict[str, str] | None = None):
    """Wrap the body of a job. On success, writes the sidecar with timings.

    Anything the body puts in the yielded dict is merged into the sidecar, so
    a stage can record its own counts (`m["n_pairs"] = ...`).

    On failure nothing is written: an output with no sidecar, or no output at
    all, both read as "not done" to `aoucov status`, which is what you want.
    """
    from .store import get_store

    store = store or get_store(target)
    extra: dict[str, Any] = {}
    started = _now()
    t0 = time.time()
    yield extra
    meta = collect_meta(cfg, stage, target, item=item, tools=tools, extra=extra)
    meta["started_at"] = started
    meta["seconds"] = round(time.time() - t0, 1)
    write_meta(store, target, meta)
