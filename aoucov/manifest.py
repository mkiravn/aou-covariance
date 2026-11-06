"""What a stage should produce, and what is actually there.

This is the whole dependency engine. A stage is ready when the stages it
`needs` have nothing missing; a stage is done when it has nothing missing.
Re-running is free because every job skips an output that already exists, so
"resume" is just "submit the array again".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Any

from .provenance import meta_uri, read_meta
from .store import GcsStore, LocalStore, get_store, is_gs, parent


def _fields(cfg, spec, combo: dict) -> dict:
    fields = dict(cfg.params)
    fields.update(combo)
    fields["run_id"] = cfg.run_id
    return fields


def target_items(cfg, stage: str) -> list[tuple[str, dict]]:
    """[(uri, item)] for one stage, in a stable order.

    `item` is the loop variables for that output -- the shard index, the
    phenotype name -- which is what a runner passes back into the job and what
    the provenance sidecar records.
    """
    if stage not in cfg.stages:
        raise KeyError(f"run '{cfg.run_id}' has no stage '{stage}' "
                       f"(have: {', '.join(sorted(cfg.stages)) or 'none'})")
    spec = cfg.stages[stage]
    loops = cfg.resolve_over(spec.over)
    counts = {f"n_{k}": len(v) for k, v in loops.items()}
    names = list(loops)
    out: list[tuple[str, dict]] = []
    combos = [dict(zip(names, vals)) for vals in product(*(loops[n] for n in names))] or [{}]
    for combo in combos:
        fields = _fields(cfg, spec, {**combo, **counts})
        try:
            rel = "/".join(p for p in (spec.dir.format(**fields),
                                       spec.pattern.format(**fields)) if p)
        except KeyError as exc:
            raise KeyError(
                f"stage '{stage}' pattern references unknown field {exc}. "
                f"Available: {', '.join(sorted(fields))}"
            ) from None
        out.append((cfg.path(rel), combo))
    return out


def targets(cfg, stage: str) -> list[str]:
    return [uri for uri, _ in target_items(cfg, stage)]


@dataclass
class StageStatus:
    stage: str
    expected: int
    present: int
    missing: list[str]
    without_meta: int
    shas: dict[str, int]          # only populated when check_shas=True
    needs_unmet: list[str]
    stale_vs: list[str] = field(default_factory=list)   # needs remade after this stage

    @property
    def done(self) -> bool:
        return self.expected > 0 and not self.missing

    @property
    def n_missing(self) -> int:
        return len(self.missing)

    @property
    def mixed_code(self) -> bool:
        return len(self.shas) > 1


def stage_status(cfg, stage: str, store=None, check_shas: bool = False,
                 min_size: int = 1, _done_cache: dict | None = None) -> StageStatus:
    uris = targets(cfg, stage)
    if uris and is_gs(uris[0]) and (store is None or isinstance(store, LocalStore)):
        store = GcsStore()   # a stage with a gs:// dir under a local root (acaf_hm3)
    store = store or get_store(cfg.root)
    found = store.exists_many(uris, min_size=min_size)
    missing = [u for u in uris if not found[u]]
    present = [u for u in uris if found[u]]

    # sidecars come from the same listing, so this costs nothing extra
    metas = store.exists_many([meta_uri(u) for u in present], min_size=1) if present else {}
    without_meta = sum(1 for u in present if not metas.get(meta_uri(u), False))

    shas: dict[str, int] = {}
    if check_shas:
        for u in present:
            meta = read_meta(store, u)
            code = (meta or {}).get("code", {}) or {}
            sha = (code.get("git_sha") or "unknown")[:8]
            if code.get("dirty"):
                sha += "+dirty"
            shas[sha] = shas.get(sha, 0) + 1

    needs_unmet = []
    for dep in cfg.stages[stage].needs:
        if _done_cache is not None and dep in _done_cache:
            dep_done = _done_cache[dep]
        else:
            dep_done = stage_status(cfg, dep, store=store).done
            if _done_cache is not None:
                _done_cache[dep] = dep_done
        if not dep_done:
            needs_unmet.append(dep)

    return StageStatus(stage, len(uris), len(present), missing,
                       without_meta, shas, needs_unmet)


def missing(cfg, stage: str, store=None, min_size: int = 1) -> list[str]:
    return stage_status(cfg, stage, store=store, min_size=min_size).missing


def present(cfg, stage: str, store=None, min_size: int = 1) -> list[str]:
    st = stage_status(cfg, stage, store=store, min_size=min_size)
    miss = set(st.missing)
    return [u for u in targets(cfg, stage) if u not in miss]


def _written_range(cfg, stage: str, store) -> tuple[str, str] | None:
    times = []
    for u in targets(cfg, stage):
        meta = read_meta(store, u)
        if meta and meta.get("written_at"):
            times.append(meta["written_at"])
    return (min(times), max(times)) if times else None


def run_status(cfg, stages: list[str] | None = None, store=None,
               check_shas: bool = False, check_stale: bool = True) -> list[StageStatus]:
    """`check_stale` reads every sidecar: a stage is stale when a stage it
    needs was written after it."""
    store = store or get_store(cfg.root)
    cache: dict[str, bool] = {}
    out = []
    for name in (stages or list(cfg.stages)):
        st = stage_status(cfg, name, store=store, check_shas=check_shas,
                          _done_cache=cache)
        cache[name] = st.done
        out.append(st)
    if check_stale:
        ranges = {}
        for st in out:
            for name in [st.stage, *cfg.stages[st.stage].needs]:
                if name not in ranges:
                    ranges[name] = _written_range(cfg, name, store)
            mine = ranges[st.stage]
            st.stale_vs = [d for d in cfg.stages[st.stage].needs
                           if mine and ranges[d] and ranges[d][1] > mine[0]]
    return out


def pending_items(cfg, stage: str, store=None, force: bool = False
                  ) -> list[tuple[str, dict]]:
    """The (uri, item) pairs still to do -- what a runner iterates over."""
    items = target_items(cfg, stage)
    if force:
        return items
    if items and is_gs(items[0][0]) and (store is None or isinstance(store, LocalStore)):
        store = GcsStore()
    store = store or get_store(cfg.root)
    found = store.exists_many([u for u, _ in items])
    return [(u, it) for u, it in items if not found[u]]


def require_ready(cfg, stage: str, store=None) -> None:
    """Raise unless every stage `stage` needs is complete."""
    st = stage_status(cfg, stage, store=store)
    if st.needs_unmet:
        raise RuntimeError(
            f"stage '{stage}' needs {', '.join(st.needs_unmet)}, which "
            f"{'is' if len(st.needs_unmet) == 1 else 'are'} not complete. "
            f"Run `python -m aoucov --run {cfg.run_id} status` to see what is missing.")


def single_target(cfg, stage: str) -> str:
    """The one output of a stage that has no loop."""
    uris = targets(cfg, stage)
    if len(uris) != 1:
        raise ValueError(f"stage '{stage}' has {len(uris)} outputs, not one")
    return uris[0]


def needs_run(uri: str, force: bool = False, store=None) -> bool:
    """True if `uri` is missing or `force`; otherwise say who made it."""
    store = store or get_store(uri)
    if force or not store.exists(uri):
        return True
    meta = read_meta(store, uri) or {}
    sha = ((meta.get("code") or {}).get("git_sha") or "unknown")[:8]
    print(f"exists, keeping it: {uri}\n  made {meta.get('written_at', '?')} by {sha}; "
          f"set FORCE = True to remake it (then remake what depends on it)")
    return False
