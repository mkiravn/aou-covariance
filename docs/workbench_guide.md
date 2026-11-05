# Running aucov_across_pops on the Workbench

One run = one target population = one config, `runs/<target>.json`, named by its 1000G populations
(`chb_chs`, `stu_itu_beb`). `superpop` sets round 1 and round 2's PCA; `target_pops` sets round 2's gate.

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

| run | `superpop` | `target_pops` | 1000G people (2,504 set) |
|---|---|---|---|
| `chb_chs` | EAS | CHB+CHS | 103 + 105 |
| `yri` | AFR | YRI | 108 |
| `lwk` | AFR | LWK | 99 |
| `stu_itu_beb` | SAS | STU+ITU+BEB | 102 + 102 + 86 |

`chb_chs` also holds `hm3_overlap`'s shared SNP list. `yri` and `lwk` share round 1 (same `superpop`), so
their round-1 output is identical. `r2_k = 4`, `r2_q = 0.999` are provisional in all of them: set them from the
round-2 diagnostics.

## Uniform sample (`uniform_10k`)

10,000 people drawn at random from everyone in `ancestry_preds` (`uniform_n`, `seed`), so the ancestry mix is
the cohort's. `uniform_sample` replaces round 1, and round 2 only fits the PCA on all of 1000G
(`superpop: "ALL"`) to show where the sample lands; there is no gate. Order:

```bash
nohup bash -c 'for nb in uniform_sample hm3_panel ancestry_round2 covariate_pcs grm; do
  python -m aoucov --run uniform_10k nb $nb || break; done
python -m aoucov --run uniform_10k sync' > ~/uniform_10k.log 2>&1 &
```

It needs `chb_chs`'s `hm3_overlap` first. A GRM in a mixed sample carries ancestry structure: on simulated data
12% of pairs fell outside the bins (−0.05 to 1.5) and the h2 checks don't apply.

## A new target population

On your laptop, not the VM (the VM never pushes):

```bash
cp runs/chb_chs.json runs/pur.json
```

Edit `run_id`, `root`, `mirror`, `superpop` (`AMR`), `target_pops` (`["PUR"]`). `hm3_ids` can keep pointing at `chb_chs`'s file: it doesn't depend on the population. Then:

```bash
git add runs/pur.json && git commit -m "Run config: PUR" && git tag run-pur-v1 && git push --tags origin HEAD
```

On the VM: `cd ~/aou-covariance && git pull`.

## Running a stage

Short stages can run straight in the terminal:

```bash
python -m aoucov --run chb_chs nb ancestry_round1
python -m aoucov --run chb_chs sync
```

Long ones (`hm3_overlap`, `hm3_panel`, `covariate_pcs`, `grm`) should survive a closed browser tab. The VM has no
`tmux`, so use `nohup` and follow the log:

```bash
nohup python -m aoucov --run chb_chs nb hm3_panel > ~/hm3_panel.chb_chs.log 2>&1 &
tail -f ~/hm3_panel.chb_chs.log
```

The executed notebook, with outputs, is in
`~/aucov_runs/chb_chs/notebooks/`; open it in JupyterLab to look at the figures and checks.
`python -m aoucov --run chb_chs status` shows what's done (`ok`), partial (`**`), missing (`--`) or stale (`~~`).

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

The configs above already have provisional `r2_k`/`r2_q`, so the gate runs. Stop after round 2. An `n1-standard-16` is enough for this part (`threads` is capped at the VM's cores).

```bash
nohup bash -c 'python -m aoucov --run chb_chs nb hm3_overlap   # once, writes the shared hm3_ids
for run in chb_chs yri lwk stu_itu_beb; do
  for nb in ancestry_round1 hm3_panel ancestry_round2; do python -m aoucov --run $run nb $nb || break; done
  python -m aoucov --run $run sync
done' > ~/ancestry.log 2>&1 &
```

`ancestry_round1` alone takes minutes (premade PCs only), and its result depends only on `superpop`. The cost
is `hm3_panel`, which reads ACAF for everyone round 1 keeps. Compare `n_kept`, the gate figures and the
kept-vs-`q` table in each run's `round2_pca` and `ancestry_round2` reports.

## Taking results out

```bash
python -m aoucov --run chb_chs report
```

writes `~/aucov_runs/chb_chs/report/RESULTS.md`: aggregate counts (1–20 masked), checks, parameters and the
commit. Review it and the figures before they leave the Workbench. Executed notebooks and anything per-person
stay in the bucket.

## When you stop

Stop the app when nothing is running (n1-highmem-64 is roughly $3.8/h). The data disk survives a stop;
the bucket `mirror` is the durable copy.
