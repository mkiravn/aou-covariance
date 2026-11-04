"""Run a stage. The backend is an argument, not a rewrite.

The same `fn` runs on a laptop, in a notebook on a Workbench VM, and inside a
dsub task; the only difference is who iterates. Outputs that already exist are
skipped, so resubmitting a half-failed array is the resume path.
"""

from __future__ import annotations

import os
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Callable

from .manifest import pending_items, target_items
from .provenance import record
from .store import get_store


def run_stage(
    cfg,
    stage: str,
    fn: Callable,
    jobs: int = 1,
    force: bool = False,
    tools: dict[str, str] | None = None,
    store=None,
    dry_run: bool = False,
    log=print,
) -> list[str]:
    """Call `fn(cfg, item, out_uri, meta)` for every output not already there.

    `meta` is a dict the job may add counts to; it ends up in the sidecar.
    Returns the list of outputs produced this call.
    """
    store = store or get_store(cfg.root)
    todo = pending_items(cfg, stage, store=store, force=force)
    total = len(target_items(cfg, stage))
    log(f"[{stage}] {len(todo)} to do of {total}")
    if dry_run:
        for uri, _ in todo:
            log(f"  would produce {uri}")
        return []

    done: list[str] = []

    def one(pair):
        uri, item = pair
        if not uri.startswith("gs://"):
            os.makedirs(os.path.dirname(uri), exist_ok=True)
        with record(cfg, stage, uri, store=store, item=item, tools=tools) as meta:
            fn(cfg, item, uri, meta)
        return uri

    if jobs <= 1:
        for pair in todo:
            done.append(one(pair))
            log(f"  [{len(done)}/{len(todo)}] {pair[0]}")
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            for uri in pool.map(one, todo):
                done.append(uri)
                log(f"  [{len(done)}/{len(todo)}] {uri}")
    return done


def write_task_tsv(cfg, stage: str, path: str, force: bool = False,
                   store=None, extra_env: dict[str, str] | None = None) -> int:
    """Write a dsub `--tasks` TSV covering exactly the outputs still missing.

    One column per loop variable, plus the output URI, plus the run id and
    stage so the task script can re-derive everything else from the config.
    Returns the number of tasks; zero means the stage is already complete.
    """
    todo = pending_items(cfg, stage, store=store, force=force)
    if not todo:
        with open(path, "w") as fh:
            fh.write("")
        return 0
    loop_vars = list(todo[0][1].keys())
    extra_env = extra_env or {}
    header = (
        [f"--env {v.upper()}" for v in loop_vars]
        + ["--env AOUCOV_RUN", "--env AOUCOV_STAGE"]
        + [f"--env {k}" for k in extra_env]
        + ["--output OUT"]
    )
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for uri, item in todo:
            row = (
                [str(item[v]) for v in loop_vars]
                + [cfg.run_id, stage]
                + [str(v) for v in extra_env.values()]
                + [uri]
            )
            fh.write("\t".join(row) + "\n")
    return len(todo)


def task_env(cfg=None) -> dict:
    """Inside a dsub task: the item this task was handed.

    Reads back the --env columns written by write_task_tsv, so a task script
    does not need to parse anything.
    """
    out = {
        "run": os.environ.get("AOUCOV_RUN"),
        "stage": os.environ.get("AOUCOV_STAGE"),
        "out": os.environ.get("OUT"),
    }
    if cfg is not None and out["stage"] in cfg.stages:
        for var in cfg.resolve_over(cfg.stages[out["stage"]].over):
            val = os.environ.get(var.upper())
            if val is not None:
                out.setdefault("item", {})[var] = val
    return out


def banner(cfg, stage: str, file=sys.stderr) -> None:
    from .provenance import code_version

    cv = code_version()
    print(
        f"aoucov run={cfg.run_id} cfg={cfg.config_hash()} stage={stage} "
        f"sha={cv['git_sha']}{'+dirty' if cv.get('dirty') else ''}",
        file=file,
    )
