# resources

Public, non-participant files every run shares. Nothing here comes from All of Us.

| file | what |
|---|---|
| `high_ld_regions_hg38.txt` | long-range LD regions, GRCh38, in `plink --exclude range` format (`chr start end label`, no header). The canonical list plus the narrow LD-peak intervals found in the `1kg_eur` round-2 loadings; copied from `notebooks/02_round2_gate.ipynb` so there is one copy. The round-2 and covariate PCAs also detect long-range LD automatically (`bigsnpr::bed_autoSVD`), so a region missing from this list is still caught. |

Bins are not a file: each run generates them from `params.bins` with `src/make_bins.py`.
