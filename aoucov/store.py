"""One existence/listing/read/write interface over local paths and gs:// URIs.

Everything else in the package asks a Store rather than touching the
filesystem or shelling out, which is what lets a stage run unchanged on a
laptop against a temp directory and on a Batch worker against a bucket.

Existence is checked by listing a prefix once and diffing, never one
stat per object: a 300-shard stage is one `gcloud storage ls`, not 300.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Iterable


def is_gs(uri: str) -> bool:
    return str(uri).startswith("gs://")


def parent(uri: str) -> str:
    return uri.rsplit("/", 1)[0]


class Store:
    def listing(self, prefix: str) -> dict[str, int]:
        """Map of uri -> size in bytes for everything under `prefix`."""
        raise NotImplementedError

    def exists(self, uri: str, min_size: int = 1) -> bool:
        raise NotImplementedError

    def read_text(self, uri: str) -> str:
        raise NotImplementedError

    def write_text(self, uri: str, text: str) -> None:
        raise NotImplementedError

    # Shared: resolve many URIs with as few listings as possible.
    def exists_many(self, uris: Iterable[str], min_size: int = 1) -> dict[str, bool]:
        uris = list(uris)
        sizes: dict[str, int] = {}
        for prefix in sorted({parent(u) for u in uris}):
            sizes.update(self.listing(prefix))
        return {u: sizes.get(u, -1) >= min_size for u in uris}


class LocalStore(Store):
    def listing(self, prefix: str) -> dict[str, int]:
        out: dict[str, int] = {}
        if not os.path.isdir(prefix):
            return out
        for root, _dirs, files in os.walk(prefix):
            for f in files:
                p = os.path.join(root, f)
                try:
                    out[p] = os.path.getsize(p)
                except OSError:
                    pass
        return out

    def exists(self, uri: str, min_size: int = 1) -> bool:
        try:
            return os.path.getsize(uri) >= min_size
        except OSError:
            return False

    def read_text(self, uri: str) -> str:
        with open(uri) as fh:
            return fh.read()

    def write_text(self, uri: str, text: str) -> None:
        os.makedirs(os.path.dirname(uri), exist_ok=True)
        tmp = uri + ".tmp"
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, uri)


class GcsStore(Store):
    """Shells out to `gcloud storage` (or `gsutil` if that is all there is).

    Listings are cached for the life of the object, so a status pass over a
    run costs one call per stage directory. Make a fresh GcsStore to re-check.
    """

    def __init__(self, cli: list[str] | None = None):
        self._cache: dict[str, dict[str, int]] = {}
        if cli is not None:
            self.cli = cli
        elif shutil.which("gcloud"):
            self.cli = ["gcloud", "storage"]
        elif shutil.which("gsutil"):
            self.cli = ["gsutil"]
        else:  # pragma: no cover - environment dependent
            raise RuntimeError("neither `gcloud` nor `gsutil` is on PATH")

    def _run(self, args: list[str]) -> str:
        proc = subprocess.run(
            self.cli + args, capture_output=True, text=True, check=False
        )
        if proc.returncode != 0:
            # "no URLs matched" is an empty prefix, not a failure
            blob = (proc.stderr or "").lower()
            if "matched no objects" in blob or "no urls matched" in blob or "not found" in blob:
                return ""
            raise RuntimeError(
                f"{' '.join(self.cli + args)} failed ({proc.returncode}): {proc.stderr.strip()}"
            )
        return proc.stdout

    def listing(self, prefix: str) -> dict[str, int]:
        prefix = prefix.rstrip("/")
        if prefix in self._cache:
            return self._cache[prefix]
        out: dict[str, int] = {}
        for line in self._run(["ls", "-l", f"{prefix}/**"]).splitlines():
            line = line.strip()
            if not line or line.startswith("TOTAL:"):
                continue
            parts = line.split()
            uri = parts[-1]
            if not is_gs(uri):
                continue
            try:
                out[uri] = int(parts[0])
            except ValueError:
                out[uri] = 0
        self._cache[prefix] = out
        return out

    def exists(self, uri: str, min_size: int = 1) -> bool:
        return self.listing(parent(uri)).get(uri, -1) >= min_size

    def read_text(self, uri: str) -> str:
        return self._run(["cat", uri])

    def write_text(self, uri: str, text: str) -> None:
        proc = subprocess.run(
            self.cli + ["cp", "-", uri], input=text, text=True,
            capture_output=True, check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"writing {uri} failed: {proc.stderr.strip()}")


def get_store(uri: str) -> Store:
    return GcsStore() if is_gs(uri) else LocalStore()
