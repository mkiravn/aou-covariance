# aoucov — run bookkeeping

Four small pieces, no dependencies, same behaviour on a laptop, a Workbench
VM, and inside a dsub task:

| module | answers |
|---|---|
| `config.py` | what is this run? (a committed JSON/YAML file + a run id) |
| `store.py` | uniform `exists / list / read / write` over local paths and `gs://` |
| `manifest.py` | what should this stage produce, and what is actually there? |
| `provenance.py` | which code, config and machine produced this output? |
| `runner.py` | run a stage locally, in parallel, or as a dsub task array |

It does not submit jobs, know what a GRM is, or sequence anything. It answers
*what is done, what is missing, and where did this file come from*.

## A run is a config file

`runs/<run_id>.json` (see [`runs/1kg_eur.json`](../runs/1kg_eur.json)), committed
to git. Every output path is derived from it, so no notebook hardcodes a bucket
path again. Stages are declarative:

```json
"crossproducts": {
  "dir": "06_crossproducts/{pheno}",
  "pattern": "acc.{i}.tsv",
  "over": {"pheno": "@phenotypes", "i": "$n_shards"},
  "needs": ["grm_shards", "phenotypes", "relatedness"]
}
```

`@name` pulls a list from `lists`, `$name` expands an integer param to `1..n`.
Patterns may use any loop variable, any `params` entry, `run_id`, and `n_<var>`.

## The command you run at the start of a session

```bash
python -m aoucov --run 1kg_eur status
```

```
run 1kg_eur  config 4b1c9e2f0a77  root gs://…/1kg_eur
stage          have  want  missing  no-meta  blocked-by
-------------  ----  ----  -------  -------  ------------
ok grm_shards    20    20        0
ok phenotypes     2     2        0
-- relatedness    0     1        1
** crossproducts 31    40        9        4  relatedness
-- merged         0     2        2        
```

`ok` done · `**` partial · `--` not started. `no-meta` is outputs with no
provenance sidecar — usually things made by hand before this layer existed.
Add `--shas` to read every sidecar and see which commits produced a stage;
it warns when one stage's outputs came from more than one commit, which is how
an N-mismatch between rungs goes unnoticed.

## Running a stage

One function, three backends:

```python
from aoucov import load_run
from aoucov.runner import run_stage

cfg = load_run("1kg_eur")

def accumulate(cfg, item, out_uri, meta):
    ...                      # produce out_uri
    meta["n_pairs"] = n      # ends up in the sidecar

run_stage(cfg, "crossproducts", accumulate, jobs=8)   # VM, 8 at a time
```

Outputs that already exist are skipped, so **resubmitting is the resume path** —
there is no resume logic to get wrong. `force=True` redoes everything;
`dry_run=True` just lists.

For dsub, generate a `--tasks` TSV covering exactly what is missing:

```bash
python -m aoucov --run 1kg_eur tasks crossproducts -o tasks.tsv
```

One column per loop variable, plus `AOUCOV_RUN`, `AOUCOV_STAGE` and `OUT`. The
task script reads them back with `aoucov.runner.task_env(cfg)` and re-derives
everything else from the config, so the TSV never carries parameters that
could drift from the committed run.

## Provenance

Every output gets `<output>.meta.json` beside it: run id, config hash, stage,
loop item, git sha (+ dirty flag), branch, container image, machine type,
hostname, dsub job/task id, wall time, tool versions, plus anything the job
added. Nothing is written if the job raises — a missing sidecar reads as
"not done".

```bash
python -m aoucov --run 1kg_eur meta gs://…/06_crossproducts/height/acc.7.tsv
```

In a container there is no `.git`, so **the job script must export
`AOUCOV_GIT_SHA`** when it installs the pinned package. Also useful:
`AOUCOV_IMAGE`, `AOUCOV_MACHINE_TYPE`, `AOUCOV_BACKEND`.

## Tests

```bash
python -m unittest discover -s aoucov/tests -t .
```

24 tests, under a second, all against a temp directory. The rule that makes
code changes safe: every stage must run end to end on synthetic data locally
before it costs a Batch job.
