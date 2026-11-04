#!/bin/bash
# Example of the in-container half of a stage.
#
# Launch with a TSV from:  python -m aoucov --run <id> tasks <stage> -o tasks.tsv
# dsub passes one row per task as --env columns plus --output OUT.
set -euo pipefail

: "${AOUCOV_RUN:?}" "${AOUCOV_STAGE:?}" "${OUT:?}"

# Pin the code, and tell the sidecar what it is. Nothing else should install
# aoucov, so the sha in the meta.json is always the sha that ran.
: "${AOUCOV_GIT_SHA:?set this to the commit the job is pinned to}"
pip install --quiet "git+https://github.com/<org>/aou-covariance@${AOUCOV_GIT_SHA}"
export AOUCOV_BACKEND=dsub
export AOUCOV_IMAGE="${IMAGE:-unknown}"

# The config travels with the code, so the task re-derives its own paths
# instead of trusting parameters pasted into the TSV.
python - <<'PY'
import os
from aoucov import load_run
from aoucov.provenance import record
from aoucov.runner import task_env

cfg = load_run(os.environ["AOUCOV_RUN"])
env = task_env(cfg)
out = env["out"]
item = env.get("item", {})

with record(cfg, env["stage"], out, item=item,
            tools={"grm_class_tool": os.popen("grm_class_tool --version").read().strip()}) as meta:
    # ... do the work, write to `out` ...
    meta["n_pairs"] = 0
PY
