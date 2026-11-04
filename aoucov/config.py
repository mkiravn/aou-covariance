"""A run is a committed config file: an id, a root prefix, parameters, and
declarative stage definitions. Nothing else in the codebase should hardcode a
bucket path."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Any


def _load_file(path: str) -> dict:
    with open(path) as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml  # optional; JSON configs need no dependency at all
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise RuntimeError(
                f"{path} is YAML but PyYAML is not installed. "
                "Either `pip install pyyaml` or write the config as JSON."
            ) from exc
        return yaml.safe_load(text)
    return json.loads(text)


@dataclass
class StageSpec:
    """How to name one stage's outputs.

    dir / pattern are format strings. Available fields are the loop variables
    named in `over`, a `n_<var>` count for each, every entry of the run's
    `params`, and `run_id`.
    """

    name: str
    dir: str = ""
    pattern: str = ""
    over: dict[str, Any] = field(default_factory=dict)
    needs: list[str] = field(default_factory=list)
    description: str = ""


@dataclass
class RunConfig:
    run_id: str
    root: str
    params: dict[str, Any] = field(default_factory=dict)
    lists: dict[str, list] = field(default_factory=dict)
    stages: dict[str, StageSpec] = field(default_factory=dict)
    source_path: str | None = None

    # ---- paths -------------------------------------------------------------

    def path(self, *parts: str) -> str:
        """Join onto the run root. Absolute parts and gs:// URIs pass through."""
        if not parts:
            return self.root
        head = str(parts[0])
        if head.startswith(("gs://", "/")):
            base = head
            rest = parts[1:]
        else:
            base = self.root
            rest = parts
        out = base.rstrip("/")
        for p in rest:
            out = f"{out}/{str(p).strip('/')}"
        return out

    # ---- parameters and inputs ---------------------------------------------

    def param(self, name: str) -> Any:
        """A parameter; null means "not decided yet" and raises."""
        if name not in self.params:
            raise KeyError(f"run '{self.run_id}' has no param '{name}'")
        value = self.params[name]
        if value is None:
            raise ValueError(
                f"param '{name}' is still null in {self.source_path or 'the config'}: "
                f"decide it (the stage that needs it reports what to look at) and "
                f"write it into the config before running further")
        return value

    def input(self, name: str, **fields: Any) -> str:
        """An external input from `params.inputs` (the only paths that differ
        between the Workbench and the simulated run)."""
        inputs = self.params.get("inputs") or {}
        if name not in inputs:
            raise KeyError(f"run '{self.run_id}' has no input '{name}' "
                           f"(have: {', '.join(sorted(inputs)) or 'none'})")
        value = str(inputs[name]).format(**fields)
        return _resolve_local(value)

    def tool(self, name: str) -> str:
        """Path to a binary from `params.tools` (plink2, plink, grm_class_tool)."""
        tools = self.params.get("tools") or {}
        if name not in tools:
            raise KeyError(f"run '{self.run_id}' has no tool '{name}'")
        return _resolve_local(str(tools[name]))

    # ---- identity ----------------------------------------------------------

    def canonical(self) -> dict:
        """The part of the config that affects results. `source_path` is where
        the file happened to sit, so it is deliberately excluded."""
        return {
            "run_id": self.run_id,
            "root": self.root,
            "params": self.params,
            "lists": self.lists,
            "stages": {
                k: {"dir": v.dir, "pattern": v.pattern, "over": v.over}
                for k, v in sorted(self.stages.items())
            },
        }

    def config_hash(self) -> str:
        blob = json.dumps(self.canonical(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode()).hexdigest()[:12]

    # ---- loop expansion ----------------------------------------------------

    def resolve_over(self, over: dict[str, Any]) -> dict[str, list]:
        """Turn each `over` entry into a concrete list of values.

        int or "$param"   -> 1..n
        "@listname"       -> cfg.lists[listname]
        explicit list      -> itself
        """
        out: dict[str, list] = {}
        for var, spec in over.items():
            if isinstance(spec, str) and spec.startswith("@"):
                key = spec[1:]
                if key not in self.lists:
                    raise KeyError(f"stage loop '{var}' references unknown list '{key}'")
                values = list(self.lists[key])
            elif isinstance(spec, str) and spec.startswith("$"):
                key = spec[1:]
                if key not in self.params:
                    raise KeyError(f"stage loop '{var}' references unknown param '{key}'")
                values = list(range(1, int(self.params[key]) + 1))
            elif isinstance(spec, int):
                values = list(range(1, spec + 1))
            elif isinstance(spec, (list, tuple)):
                values = list(spec)
            else:
                raise TypeError(f"stage loop '{var}' has unsupported spec {spec!r}")
            out[var] = values
        return out

    # ---- construction ------------------------------------------------------

    @classmethod
    def from_dict(cls, d: dict, source_path: str | None = None) -> "RunConfig":
        for required in ("run_id", "root"):
            if not d.get(required):
                raise ValueError(f"run config is missing required key '{required}'")
        stages = {}
        for name, spec in (d.get("stages") or {}).items():
            spec = dict(spec or {})
            stages[name] = StageSpec(
                name=name,
                dir=spec.get("dir", ""),
                pattern=spec.get("pattern", ""),
                over=spec.get("over", {}) or {},
                needs=list(spec.get("needs", []) or []),
                description=spec.get("description", ""),
            )
        return cls(
            run_id=d["run_id"],
            root=_resolve_local(d["root"]).rstrip("/"),
            params=d.get("params") or {},
            lists=d.get("lists") or {},
            stages=stages,
            source_path=source_path,
        )


def load_run(path: str) -> RunConfig:
    """Load a run config. `path` may be the file, or a run id resolved against
    $AOUCOV_RUNS (default: <repo>/runs)."""
    if os.path.exists(path):
        return RunConfig.from_dict(_load_file(path), source_path=os.path.abspath(path))
    runs_dir = os.environ.get("AOUCOV_RUNS") or os.path.join(_repo_root(), "runs")
    for ext in (".json", ".yaml", ".yml"):
        candidate = os.path.join(runs_dir, path + ext)
        if os.path.exists(candidate):
            return RunConfig.from_dict(_load_file(candidate), source_path=candidate)
    raise FileNotFoundError(f"no run config at '{path}' or in {runs_dir}")


def _repo_root() -> str:
    env = os.environ.get("AOUCOV_REPO")
    if env:
        return env
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _resolve_local(path: str) -> str:
    """gs:// and absolute paths pass through; relative paths are repo-relative."""
    if path.startswith("gs://"):
        return path
    path = os.path.expanduser(path)
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(_repo_root(), path))
