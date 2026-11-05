# aucov_across_pops: running locally on simulated data

Everything under `.sim/` is simulated and git-ignored. Figures from it are stamped "SIMULATED DATA".

## Once

```bash
bash scripts/dev/setup_local.sh        # .venv, R packages, plink/plink2 into .tools/bin, grm_class_tool
.venv/bin/python scripts/sim/simulate.py
```

In VS Code: open the repo folder, pick `.venv` as the interpreter/kernel. The tasks in
**Terminal → Run Task** (setup, simulate data, status, run notebook, report, tests) do the same as the commands below.

## The loop

```bash
.venv/bin/python -m aoucov --run sim status
```

Shows each stage: `ok` done, `**` partial, `--` not started, `~~` stale (something it needs was remade after it),
what it's blocked by, and its check tally.

`bash scripts/sim/run_all.sh` simulates and runs every stage in order (about 2 minutes).

Run the next unblocked stage either way:
- open `aucov_across_pops/<stage>.ipynb`, keep `RUN_ID = "sim"`, Run All, or
- `.venv/bin/python -m aoucov --run sim nb <stage>` (executed copy goes to `.sim/runs/sim/notebooks/`).

Order: `ancestry_round1` → `hm3_panel` → (`hm3_overlap`, optional diagnostic) → `ancestry_round2` → `covariate_pcs` → `grm` → `phenotypes` → `crossproducts` → `curves`.

## Decisions go in the config

`ancestry_round2` stops at `cfg.param("r2_k")` while it's `null`. Read the `pc_diagnostics` table and the loadings plot,
set `r2_k` and `r2_q` in `runs/sim.json`, rerun. A finished output is kept on rerun unless `FORCE = True`.

## What you get

- `python -m aoucov --run sim report` → `.sim/runs/sim/report/RESULTS.md` (open with Markdown preview): notes, checks, figures per stage.
- `python -m aoucov --run sim meta <output path>` → which commit, config and machine made it.
- The `curves` stage checks that `h2_Unrel` recovers the simulated h2 (`.sim/inputs/truth.json`).

## On the Workbench

Same notebooks with `RUN_ID = "eur"`. Fill the `TODO` inputs in `runs/eur.json`; run `hm3_overlap` first, its
`hm3_ids_matched.txt` is the `hm3_ids` input. Variant IDs are `chr:pos:ref:alt` without "chr" (plink2 `@:#:$r:$a`). The `%%bash` and R cells
run on the VM; at full scale the panel and GRM-shard cells need Batch (dsub), which isn't wired yet.
