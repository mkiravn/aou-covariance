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

## Runs that exist

`ceu_gbr` (EUR → CEU+GBR), `chb_jpt` (EAS → CHB+JPT), `yri` (AFR → YRI). `r2_k = 4`, `r2_q = 0.999` in all three
are provisional: set them from the round-2 diagnostics.

## A new target population

On your laptop, not the VM (the VM never pushes):

```bash
cp runs/ceu_gbr.json runs/chb_jpt.json
```

Edit `run_id`, `root`, `mirror`, `superpop` (`EAS`), `target_pops` (`["CHB", "JPT"]`). `hm3_ids` can keep pointing at `ceu_gbr`'s file: it doesn't depend on the population. Then:

```bash
git add runs/chb_jpt.json && git commit -m "Run config: CHB+JPT" && git tag run-chb_jpt-v1 && git push --tags origin HEAD
```

On the VM: `cd ~/aou-covariance && git pull`.

## Running a stage

Short stages can run straight in the terminal:

```bash
python -m aoucov --run chb_jpt nb ancestry_round1
python -m aoucov --run chb_jpt sync
```

Long ones (`hm3_overlap`, `hm3_panel`, `covariate_pcs`, `grm`) should survive a closed browser tab. The VM has no
`tmux`, so use `nohup` and follow the log:

```bash
nohup python -m aoucov --run chb_jpt nb hm3_panel > ~/hm3_panel.chb_jpt.log 2>&1 &
tail -f ~/hm3_panel.chb_jpt.log
```

The executed notebook, with outputs, is in
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

The three configs above already have provisional `r2_k`/`r2_q`, so the gate runs. Stop after round 2. An `n1-standard-16` is enough for this part (`threads` is capped at the VM's cores).

```bash
nohup bash -c 'python -m aoucov --run ceu_gbr nb hm3_overlap   # once, writes the shared hm3_ids
for run in ceu_gbr chb_jpt yri; do
  for nb in ancestry_round1 hm3_panel ancestry_round2; do python -m aoucov --run $run nb $nb || break; done
  python -m aoucov --run $run sync
done' > ~/ancestry.log 2>&1 &
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
