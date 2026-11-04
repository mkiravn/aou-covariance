"""Per-stage reports in {root}/report/{stage}/: summary.json, figures, tables.

Counts of 1-20 are masked (AoU dissemination policy); simulated runs are
stamped on every figure. Not an egress review.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from typing import Any, Iterable

from .provenance import _now, code_version
from .store import get_store, is_gs

SMALL = 20
MASK = "≤20"


def mask_count(n: Any) -> Any:
    """AoU: counts of 1-20 are never shown. Zero is fine."""
    try:
        v = int(n)
    except (TypeError, ValueError):
        return n
    return MASK if 0 < v <= SMALL else v


def mask_small(df, cols: Iterable[str]):
    out = df.copy()
    for c in cols:
        out[c] = out[c].map(mask_count).astype(object)
    return out


def _put(local: str, uri: str) -> None:
    if is_gs(uri):
        subprocess.run(["gcloud", "storage", "cp", "-q", local, uri], check=True)
    else:
        os.makedirs(os.path.dirname(uri), exist_ok=True)
        shutil.copyfile(local, uri)


class Report:
    def __init__(self, cfg, stage: str, store=None):
        if stage not in cfg.stages:
            raise KeyError(f"run '{cfg.run_id}' has no stage '{stage}'")
        self.cfg = cfg
        self.stage = stage
        self.dir = cfg.path("report", stage)
        self.store = store or get_store(self.dir)
        self.simulated = bool(cfg.params.get("simulated"))
        self.notes: dict[str, Any] = {}
        self.checks: list[dict] = []
        self.files: list[dict] = []
        self._tmp = tempfile.mkdtemp(prefix=f"report_{stage}_")

    # ---- content -----------------------------------------------------------

    def note(self, **values: Any) -> None:
        """Parameters and decisions. Not for counts."""
        self.notes.update(values)

    def count(self, name: str, n: Any) -> Any:
        """A participant or pair count, masked if 1-20."""
        self.notes[name] = mask_count(n)
        return self.notes[name]

    def check(self, name: str, ok: bool, detail: str = "") -> bool:
        """Recorded, not raised: a failure shows in `status` and RESULTS.md."""
        ok = bool(ok)
        self.checks.append({"check": name, "ok": ok, "detail": detail})
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
        return ok

    def fig(self, fig, name: str, caption: str = "") -> str:
        if self.simulated:
            fig.text(0.995, 0.005, "SIMULATED DATA - not All of Us results",
                     ha="right", va="bottom", fontsize=9, color="#c0392b",
                     alpha=0.9, weight="bold")
        local = os.path.join(self._tmp, f"{name}.png")
        fig.savefig(local, dpi=130, bbox_inches="tight")
        uri = f"{self.dir}/{name}.png"
        _put(local, uri)
        self.files.append({"file": f"{name}.png", "kind": "figure", "caption": caption})
        return uri

    def table(self, df, name: str, counts: Iterable[str] = (), caption: str = ""):
        """Write a table with `counts` columns masked; returns what to display."""
        shown = mask_small(df, counts)
        local = os.path.join(self._tmp, f"{name}.tsv")
        shown.to_csv(local, sep="\t", index=False)
        _put(local, f"{self.dir}/{name}.tsv")
        self.files.append({"file": f"{name}.tsv", "kind": "table", "caption": caption})
        return shown

    # ---- output ------------------------------------------------------------

    def save(self) -> int:
        """Write summary.json. Returns the number of failed checks."""
        failed = [c for c in self.checks if not c["ok"]]
        summary = {
            "run_id": self.cfg.run_id,
            "stage": self.stage,
            "config_hash": self.cfg.config_hash(),
            "code": code_version(),
            "simulated": self.simulated,
            "written_at": _now(),
            "notes": self.notes,
            "checks": self.checks,
            "files": self.files,
        }
        self.store.write_text(f"{self.dir}/summary.json",
                              json.dumps(summary, indent=2, default=str) + "\n")
        print(f"[{self.stage}] report -> {self.dir}  "
              f"checks {len(self.checks) - len(failed)}/{len(self.checks)} passed")
        return len(failed)


# ---- run-level ---------------------------------------------------------------

def read_summary(cfg, stage: str, store=None) -> dict | None:
    uri = cfg.path("report", stage, "summary.json")
    store = store or get_store(uri)
    try:
        if not store.exists(uri):
            return None
        return json.loads(store.read_text(uri))
    except Exception:
        return None


def check_tally(summary: dict | None) -> str:
    if not summary or not summary.get("checks"):
        return ""
    n = len(summary["checks"])
    ok = sum(1 for c in summary["checks"] if c["ok"])
    return f"{ok}/{n}" + ("" if ok == n else " !")


def write_results(cfg, store=None) -> str:
    """Stitch every stage's summary.json into report/RESULTS.md."""
    lines = [f"# {cfg.run_id}: results", ""]
    sims = []
    body = []
    for stage in cfg.stages:
        s = read_summary(cfg, stage, store=store)
        if s is None:
            continue
        sims.append(s.get("simulated"))
        code = s.get("code") or {}
        sha = (code.get("git_sha") or "unknown")[:8] + ("+dirty" if code.get("dirty") else "")
        body += [f"## {stage}", "",
                 f"`{sha}` · config `{s['config_hash']}` · {s['written_at']}", ""]
        if s["notes"]:
            body += ["| | |", "|---|---|"]
            body += [f"| {k} | {v} |" for k, v in s["notes"].items()]
            body.append("")
        if s["checks"]:
            body += ["| check | | detail |", "|---|---|---|"]
            body += [f"| {c['check']} | {'✓' if c['ok'] else '✗'} | {c['detail']} |"
                     for c in s["checks"]]
            body.append("")
        for f in s["files"]:
            if f["kind"] == "figure":
                body += [f"![{f['caption'] or f['file']}]({stage}/{f['file']})", ""]
            else:
                body += [f"- [{f['file']}]({stage}/{f['file']}) {f['caption']}"]
        body.append("")
    if any(sims):
        lines += ["> **Simulated data.** Nothing here is an All of Us result.", ""]
    lines += body
    uri = cfg.path("report", "RESULTS.md")
    (store or get_store(uri)).write_text(uri, "\n".join(lines) + "\n")
    return uri
