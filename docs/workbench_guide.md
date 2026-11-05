# Running aucov_across_pops on the Workbench

One run = one target population = one config, `runs/<target>.json`, named by its 1000G populations
(`ceu_gbr`, `chb_jpt`, `yri`). `superpop` sets round 1 and round 2's PCA; `target_pops` sets round 2's gate.

## Once

1. Workbench → Apps → new JupyterLab app: `n1-highmem-64`, 1 TB data disk. Open a terminal.
2. Clone and set up:
   ```bash
   git clone https://github.com/mkiravn/aou-covariance.git ~/aou-covariance
   cd ~/aou-covariance && bash scripts/dev/setup_vm.sh
   ```
3. Check the mounts and disk (note the trailing space in `Controlled Tier `):
   ```bash
   ls ~/workspace/cdrv9/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/acaf_threshold/pgen | head -3
   ls ~/workspace/"Data from All of Us Controlled Tier "/shared-env-pilot/1000g_reference/1kg_all_qc.bed
   df -h ~
   ```

## A new target population

On your laptop, not the VM (the VM never pushes):

```bash
cp runs/ceu_gbr.json runs/chb_jpt.json
```

Edit `run_id`, `root`, `mirror`, `superpop` (`EAS`), `target_pops` (`["CHB", "JPT"]`), and set `r2_k`, `r2_q` to
`null`. `hm3_ids` can keep pointing at `ceu_gbr`'s file: it doesn't depend on the population. Then:

```bash
git add runs/chb_jpt.json && git commit -m "Run config: CHB+JPT" && git tag run-chb_jpt-v1 && git push --tags origin HEAD
```

On the VM: `cd ~/aou-covariance && git pull`.

## Running a stage

Run stages inside `tmux` so they survive a closed browser tab:

```bash
tmux new -s aucov
python -m aoucov --run chb_jpt nb ancestry_round1
python -m aoucov --run chb_jpt sync
```

Detach with `Ctrl-b d`, reattach with `tmux attach -t aucov`. Without tmux: `nohup python -m aoucov --run chb_jpt nb grm > grm.log 2>&1 &`. The executed notebook, with outputs, is in
`~/aucov_runs/chb_jpt/notebooks/`; open it in JupyterLab to look at the figures and checks.
`python -m aoucov --run chb_jpt status` shows what's done (`ok`), partial (`**`), missing (`--`) or stale (`~~`).

## Order, and what to look at

| stage | look at | decide |
|---|---|---|
| `hm3_overlap` | SNPs lost at each step. Once per CDR, not per population | nothing |
| `ancestry_round1` | coverage curve, gate on PCs 1–4 | `r1_width` if the gate cuts into the cloud |
| `hm3_panel` | per-chromosome times (reads ACAF through the mount) | if chr1 takes hours, stop: this needs dsub |
| `ancestry_round2` | stops at the gate until `r2_k`, `r2_q` are set; scores by population, loadings by position | set `r2_k`, `r2_q` in the config, commit, pull, rerun |
| `covariate_pcs` | scores, LD regions found | nothing |
| `grm` | diagonal mean ≈ 1, N matches round 2 | nothing; the longest stage |
| `phenotypes` | needs `phenotype_tables` (still TODO) | |
| `crossproducts`, `curves` | curves, h2 checks | bin widths (`bins`) |

`sync` after each stage. `nb <stage> --force` remakes a stage; `status` then marks what depends on it as stale.

## Ancestry filtering only, several populations

Make the configs (above), with a provisional `r2_k` (e.g. 4) and `r2_q` (0.999) so the gate runs, then stop
after round 2. An `n1-standard-16` is enough for this part.

```bash
python -m aoucov --run ceu_gbr nb hm3_overlap   # once
for run in ceu_gbr chb_jpt yri; do
  for nb in ancestry_round1 hm3_panel ancestry_round2; do python -m aoucov --run $run nb $nb || break; done
  python -m aoucov --run $run sync
done
```

`ancestry_round1` alone takes minutes (premade PCs only), and its result depends only on `superpop`. The cost
is `hm3_panel`, which reads ACAF for everyone round 1 keeps. Compare `n_kept`, the gate figures and the
kept-vs-`q` table in each run's `round2_pca` and `ancestry_round2` reports.

## Taking results out

```bash
python -m aoucov --run chb_jpt report
```

writes `~/aucov_runs/chb_jpt/report/RESULTS.md`: aggregate counts (1–20 masked), checks, parameters and the
commit. Review it and the figures before they leave the Workbench. Executed notebooks and anything per-person
stay in the bucket.

## When you stop

Stop the app when nothing is running (n1-highmem-64 is roughly $3.8/h). The data disk survives a stop;
the bucket `mirror` is the durable copy.
