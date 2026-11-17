"""Height bin means of y_i y_j against GRM and Weir-Goudet relatedness, without the kin_pheno_cov jackknife.
Writes <run>/report/height_quick/binned_covariance.tsv (bins with 21+ pairs; SE ignores pair dependence).
    python slides/height_cov_quick.py RUN..."""
import os, sys
import numpy as np, pandas as pd
sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")]
from aoucov import load_run
from aoucov.manifest import single_target, target_items

EDGES = 2 * np.unique(np.round(np.r_[np.arange(-0.05, 0.06, 0.005), np.arange(0.06, 0.2, 0.02), [0.2, 0.3, 0.4, 0.6]], 4))
for run in sys.argv[1:]:
    cfg = load_run(run)
    shards = target_items(cfg, "grm_shards")
    ids = pd.read_csv(shards[0][0].split(".grm.bin.")[0] + ".grm.id", sep="\t", header=None, dtype=str)[1]
    n = len(ids)
    off = np.ones(n * (n + 1) // 2, bool); off[np.arange(1, n + 1) * (np.arange(1, n + 1) + 1) // 2 - 1] = False
    est = {"GRM": np.concatenate([np.fromfile(u, dtype="float32") for u, _ in shards])[off]}
    wg = single_target(cfg, "weir_goudet")
    assert list(pd.read_csv(wg.replace(".bin", ".id"), header=None, dtype=str)[0]) == list(ids)
    est["Weir-Goudet"] = (2 * np.fromfile(wg, dtype="float32"))[off]; del off
    uri = [u for u, _ in target_items(cfg, "phenotypes") if os.path.basename(u).startswith("height__raw__base_pcs")]
    if not uri: print(run, "no height__raw__base_pcs phenotype, skipped"); continue
    y = pd.read_csv(uri[0], sep=" ", dtype={"IID": str}).set_index("IID")["Y"].reindex(ids).to_numpy()
    i, j = np.tril_indices(n, -1)
    c = y[i] * y[j]; del i, j
    have = ~np.isnan(c)
    rows = []
    for lab, k in est.items():
        k, cc = k[have].astype("float64"), c[have]
        b = np.digitize(k, EDGES) - 1; ok = (b >= 0) & (b < len(EDGES) - 1)
        nb = len(EDGES) - 1
        cnt = np.bincount(b[ok], minlength=nb); sk = np.bincount(b[ok], k[ok], nb); sc = np.bincount(b[ok], cc[ok], nb)
        scc = np.bincount(b[ok], cc[ok] ** 2, nb)
        with np.errstate(invalid="ignore", divide="ignore"):
            mc = sc / cnt; se = np.sqrt((scc / cnt - mc ** 2) / cnt)
        with np.errstate(invalid="ignore", divide="ignore"): rows.append(pd.DataFrame({"pheno": "height", "estimator": lab, "lo": EDGES[:-1], "hi": EDGES[1:], "mean_rel": sk / cnt, "mean_cov": mc, "se": se, "pairs": cnt}))
    out = pd.concat(rows); out = out[out["pairs"] >= 21]   # AoU: no bins of 1-20 pairs
    d = os.path.join(cfg.path("report"), "height_quick"); os.makedirs(d, exist_ok=True)
    out.to_csv(f"{d}/binned_covariance.tsv", sep="\t", index=False)
    print(run, "wrote", f"{d}/binned_covariance.tsv", len(out), "bins")
