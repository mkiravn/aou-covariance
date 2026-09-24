# runs/

One directory per sample set. Replaces `pipelines/` + the numbered stage folders
for everything from `eur_r2` onward; earlier runs stay where they are.

Three things are different here.

## 1. Code reaches the VM by `git pull`, not copy-paste

The constraint is that **participant data never leaves the workbench**. That does
not require copy-paste — it only requires the flow to be one-way.

```bash
# once, on the workbench VM
git clone https://github.com/<you>/aou-covariance.git ~/aou-covariance

# every time you change something
git -C ~/aou-covariance pull
```

Two rules keep the guarantee:

- **The VM never pushes.** No write credentials on it at all. If you find
  yourself wanting to commit from the workbench, the thing you want to commit is
  probably an output, and outputs go to the bucket.
- **Nothing is written into the repo tree.** Every path a notebook writes comes
  from `Run.out()` or `Run.scratch()`, both of which live outside the checkout.

Because the steps are committed as `.ipynb`, **`nbstripout` is not optional** —
it is the thing standing between a run against real data and participant values
in git history:

```bash
pip install nbstripout && nbstripout --install   # once, per clone
```

Install it in the VM clone as well as locally. Worth checking with
`git diff --cached` before the first commit after a real run.

## 2. Notebooks are drivers, not programs

Steps are real `.ipynb` files — open and run them, no pasting. Prose for the
decisions, short cells that call the library. The machinery — QC, pruning, PCA,
projection, the LD-peak and batch-effect checks, plotting — lives in
`_lib/aoucov/` and is imported. It used to be copy-pasted into four notebooks,
which is why a parameter change meant four edits and three of them got forgotten.

```python
import sys, os
sys.path.insert(0, os.path.expanduser("~/aou-covariance/runs/_lib"))
from aoucov import Run, plink, refs, gate, plots
```

Cells are named for what they do. There are no `Cell 7`s: a cross-reference like
"rerun cells 3–8" rots the moment you insert a cell, and they always get inserted.

## 3. Every run has a RESULTS.md

Cleared outputs plus bucket-only artefacts means the repo otherwise records what
you *intended* to run and nothing about what happened. `provenance.result_line()`
prints a row to paste in.

Parameters and aggregate counts only. Anything that has not cleared AoU egress
review does not go in the repo — if a number can't be written down, record its
shape (`n=222,4xx`) and leave the exact value in the bucket.

`provenance.write()` additionally drops a stamped TSV beside each artefact in the
bucket, including the commit the code was pulled from, so any output can be
traced back to exact code.

## Layout

```
runs/
  _lib/aoucov/      shared code; the only place plink invocations live
  <sample_set>/
    README.md       what this run is, and the step table
    RESULTS.md      what each step actually produced
    <step>.ipynb    notebooks, named by what they do
```

## Not Snakemake

Considered and rejected. The heavy steps already dispatch to Google Batch through
dsub, which a workflow manager would not manage without rebuilding the networking
setup; and the rest of the work is genuinely interactive — look at the gate plot,
pick a radius, refit. A DAG engine is bad at that loop. The problems worth solving
here were provenance and duplication, and neither needs one.
