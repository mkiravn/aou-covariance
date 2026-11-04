"""Run a stage script (env vars in, files out) locally or, later, via dsub."""

from __future__ import annotations

import os
import shlex
import subprocess
import sys

from .config import _repo_root


def script_path(name: str) -> str:
    return os.path.join(_repo_root(), "scripts", "stages", name)


def run_script(cfg, name: str, env: dict, log_path: str | None = None) -> None:
    backend = cfg.params.get("backend", "local")
    path = script_path(name)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if backend != "local":
        raise NotImplementedError(
            f"backend '{backend}' is not wired yet; submit {name} with dsub as the "
            f"existing notebooks do, passing these as --env:\n"
            + "\n".join(f"  {k}={shlex.quote(str(v))}" for k, v in env.items()))
    cmd = ["Rscript", path] if name.endswith(".R") else ["bash", path]
    full_env = {**os.environ, **{k: str(v) for k, v in env.items()},
                "AOUCOV_RUN": cfg.run_id, "AOUCOV_REPO": _repo_root()}
    print(f"$ {' '.join(cmd)}", flush=True)
    proc = subprocess.Popen(cmd, env=full_env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    lines = []
    for line in proc.stdout:
        lines.append(line)
        sys.stdout.write(line)
        sys.stdout.flush()
    proc.wait()
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "w") as fh:
            fh.writelines(lines)
    if proc.returncode:
        raise RuntimeError(f"{name} exited {proc.returncode}")
