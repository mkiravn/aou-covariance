#!/usr/bin/env bash
# Simulate, then run every stage notebook in order for RUN (default sim).
set -euo pipefail
cd "$(dirname "$0")/../.."
RUN=${1:-sim}
[ "$RUN" = sim ] && .venv/bin/python scripts/sim/simulate.py
for nb in ancestry_round1 hm3_panel hm3_overlap ancestry_round2 covariate_pcs grm phenotypes crossproducts curves; do
  .venv/bin/python -m aoucov --run "$RUN" nb "$nb" 2> >(grep -v IPKernelApp >&2)
done
.venv/bin/python -m aoucov --run "$RUN" report
.venv/bin/python -m aoucov --run "$RUN" status
