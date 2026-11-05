"""`python -m aoucov <command>` -- the thing you run at the start of a session
to find out where a run actually got to."""

from __future__ import annotations

import argparse
import json
import sys

from .config import load_run
from .manifest import run_status, stage_status, target_items, targets
from .provenance import read_meta
from .report import check_tally, read_summary, write_results
from .runner import write_task_tsv
from .store import get_store


def _fmt_status(rows, check_shas: bool, tallies: dict | None = None) -> str:
    tallies = tallies or {}
    head = ["stage", "have", "want", "missing", "no-meta", "blocked-by", "checks"]
    if check_shas:
        head.append("code")
    out = [head]
    for r in rows:
        mark = ("~~ " if r.stale_vs else "ok ") if r.done else ("-- " if r.present == 0 else "** ")
        row = [
            mark + r.stage,
            str(r.present),
            str(r.expected),
            str(r.n_missing),
            str(r.without_meta) if r.without_meta else "",
            ",".join(r.needs_unmet) or ("stale vs " + ",".join(r.stale_vs) if r.stale_vs else ""),
            tallies.get(r.stage, ""),
        ]
        if check_shas:
            row.append(", ".join(f"{k}:{v}" for k, v in sorted(r.shas.items())))
        out.append(row)
    widths = [max(len(r[i]) for r in out) for i in range(len(head))]
    lines = ["  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip() for r in out]
    lines.insert(1, "  ".join("-" * w for w in widths))
    return "\n".join(lines)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="aoucov", description=__doc__)
    p.add_argument("--run", "-r", required=True,
                   help="path to a run config, or a run id under runs/")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("status", help="what is done and what is missing")
    s.add_argument("stages", nargs="*")
    s.add_argument("--shas", action="store_true",
                   help="read every sidecar and report which commits produced "
                        "the outputs (slow; one read per object)")

    t = sub.add_parser("targets", help="list what a stage should produce")
    t.add_argument("stage")
    t.add_argument("--items", action="store_true", help="show loop variables too")

    m = sub.add_parser("missing", help="list what a stage has not produced")
    m.add_argument("stage")

    k = sub.add_parser("tasks", help="write a dsub --tasks TSV for what is missing")
    k.add_argument("stage")
    k.add_argument("--out", "-o", required=True)
    k.add_argument("--force", action="store_true", help="include outputs that exist")

    v = sub.add_parser("meta", help="show the provenance sidecar for one output")
    v.add_argument("uri")

    sub.add_parser("config", help="print the resolved config and its hash")
    sub.add_parser("report", help="stitch stage summaries into report/RESULTS.md")
    r = sub.add_parser("nb", help="execute aucov_across_pops/<name>.ipynb for this run (nbclient)")
    r.add_argument("name")
    r.add_argument("--force", action="store_true")
    y = sub.add_parser("sync", help="copy a local root to params.mirror (gcloud storage rsync)")
    y.add_argument("--dry-run", action="store_true")

    a = p.parse_args(argv)
    cfg = load_run(a.run)

    if a.cmd == "status":
        rows = run_status(cfg, a.stages or None, check_shas=a.shas)
        tallies = {r.stage: check_tally(read_summary(cfg, r.stage)) for r in rows}
        print(f"run {cfg.run_id}  config {cfg.config_hash()}  root {cfg.root}")
        print(_fmt_status(rows, a.shas, tallies))
        if any(r.stale_vs for r in rows):
            print("\n~~ stale: a stage it needs was remade after it; rerun it with FORCE")
        mixed = [r.stage for r in rows if r.mixed_code]
        if mixed:
            print(f"\nwarning: outputs from more than one commit in: {', '.join(mixed)}")
        return 0

    if a.cmd == "targets":
        for uri, item in target_items(cfg, a.stage):
            print(f"{uri}\t{json.dumps(item)}" if a.items else uri)
        return 0

    if a.cmd == "missing":
        st = stage_status(cfg, a.stage)
        for uri in st.missing:
            print(uri)
        print(f"# {st.n_missing} of {st.expected} missing", file=sys.stderr)
        return 0

    if a.cmd == "tasks":
        n = write_task_tsv(cfg, a.stage, a.out, force=a.force)
        print(f"{n} tasks -> {a.out}", file=sys.stderr)
        return 0 if n else 1

    if a.cmd == "meta":
        meta = read_meta(get_store(a.uri), a.uri)
        if meta is None:
            print(f"no sidecar for {a.uri}", file=sys.stderr)
            return 1
        print(json.dumps(meta, indent=2, sort_keys=True))
        return 0

    if a.cmd == "nb":
        # nbclient ships with Jupyter; papermill drags in gcsfs and google-cloud upgrades
        import os
        import nbformat
        from nbclient import NotebookClient
        from .config import _repo_root

        src = f"{_repo_root()}/aucov_across_pops/{a.name}.ipynb"
        out = cfg.path("notebooks", f"{a.name}.ipynb")
        if out.startswith("gs://"):
            print("nb needs a local root", file=sys.stderr)
            return 1
        os.makedirs(os.path.dirname(out), exist_ok=True)
        nb = nbformat.read(src, as_version=4)
        cell = next(c for c in nb.cells if "parameters" in c.metadata.get("tags", []))
        cell.source = f"RUN_ID = {cfg.run_id!r}\nFORCE = {a.force}"
        class Live(NotebookClient):   # cell output to the terminal/log as it happens, not only at the end
            def output(self, outs, msg, display_id, cell_index):
                if msg["msg_type"] == "stream":
                    sys.stdout.write(msg["content"]["text"]); sys.stdout.flush()
                return super().output(outs, msg, display_id, cell_index)

        client = Live(nb, timeout=None, kernel_name="python3",
                      resources={"metadata": {"path": f"{_repo_root()}/aucov_across_pops"}})
        try:
            client.execute()
        finally:
            nbformat.write(nb, out)   # failed runs keep their outputs up to the error
        print(f"executed -> {out}")
        return 0

    if a.cmd == "sync":
        import subprocess
        mirror = cfg.params.get("mirror")
        if not mirror or cfg.root.startswith("gs://"):
            print("needs a local root and params.mirror", file=sys.stderr)
            return 1
        import shutil
        if not shutil.which("gcloud"):
            print("gcloud not found", file=sys.stderr)
            return 1
        cmd = ["gcloud", "storage", "rsync", "--recursive", "--exclude", r".*\.ipynb_checkpoints.*", cfg.root, mirror.rstrip("/")]
        if a.dry_run:
            cmd.append("--dry-run")
        print(" ".join(cmd), file=sys.stderr)
        return subprocess.run(cmd).returncode

    if a.cmd == "report":
        print(write_results(cfg))
        return 0

    if a.cmd == "config":
        print(json.dumps(cfg.canonical(), indent=2, sort_keys=True))
        print(f"\n# config_hash {cfg.config_hash()}", file=sys.stderr)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
