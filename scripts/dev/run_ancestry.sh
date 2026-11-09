#!/usr/bin/env bash
# Workbench, start to finish: setup, the shared HM3 list and ACAF extract (dsub), ancestry filtering for every run, sync.
#   cd ~/aou-covariance && git pull && nohup bash scripts/dev/run_ancestry.sh > ~/ancestry.log 2>&1 &
# A failed stage stops that run only; the others carry on. Finished stages are skipped on a rerun.
set -uo pipefail
cd "$(dirname "$0")/../.."
RUNS=${RUNS-"chb_chs yri lwk stu_itu_beb pjl asw mxl"}
SAMPLES=${SAMPLES-"uniform_10k:uniform_sample even_10k:even_sample"}   # run:notebook that replaces round 1
OVERLAP_RUN=${OVERLAP_RUN:-chb_chs}
SETUP=${SETUP:-1}
say() { echo "=== $(date +%T) $*"; }

if [ "$SETUP" = 1 ]; then
  say setup
  bash scripts/dev/setup_vm.sh || { say "setup failed, stopping"; exit 1; }
fi

# the matched HM3 list every run reads; remade if it predates the multi-allelic split
f=$(python -m aoucov --run "$OVERLAP_RUN" targets hm3_overlap)
grep -q "multi-allelic rows split" "$f" 2>/dev/null && force="" || force="--force"
if [ -s "$f" ] && [ -z "$force" ]; then
  say "hm3_overlap already done"
else
  say "hm3_overlap $force"
  python -m aoucov --run "$OVERLAP_RUN" nb hm3_overlap $force || { say "hm3_overlap failed, stopping"; exit 1; }
fi

# ACAF cut to HM3 for everyone, one Batch task per chromosome; every run's panel reads it
say "acaf_hm3 $force"
python -m aoucov --run "$OVERLAP_RUN" nb acaf_hm3 $force || { say "acaf_hm3 failed, stopping"; exit 1; }

stages() {   # run, stage notebooks...
  local run=$1; shift
  for nb in "$@"; do
    say "$run $nb"
    python -m aoucov --run "$run" nb "$nb" || { say "!!! $run $nb failed, skipping the rest of $run"; break; }
  done
  python -m aoucov --run "$run" sync || say "!!! $run sync failed"
}
for run in $RUNS; do stages "$run" ancestry_round1 hm3_panel ancestry_round2; done
for rs in $SAMPLES; do stages "${rs%%:*}" "${rs#*:}" hm3_panel ancestry_round2; done

say "done"
for run in $RUNS $(for rs in $SAMPLES; do echo "${rs%%:*}"; done); do
  python -m aoucov --run "$run" status ancestry_round1 hm3_panel round2_pca ancestry_round2
done
