#!/usr/bin/env python3
"""Synthetic merged full/jk files with known generating values: every estimator
must recover them exactly when noise-free, and jackknife SEs must be finite and
positive once block-level noise is added."""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)
sys.path.insert(0, SRC)
import he_estimators as E  # noqa: E402

NBLOCKS = 50
H2, B_HS, B_FS, B_PO = 0.5, 0.02, 0.05, 0.0     # intercepts on the cross-product scale


def truth(cls, a):
    m = H2 * a
    if cls == "other":
        return m + np.where((a >= 0.2) & (a < 0.4), B_HS, 0.0)
    return m + (B_FS if cls == "FS" else B_PO)


def present(cls, a):
    related = (a >= 0.4) & (a < 0.6)
    return ~related if cls == "other" else related


def build(out_dir, mid, noise):
    rng = np.random.default_rng(1)
    full, jk = [], []
    for cls in E.CLS:
        for k in np.where(present(cls, mid))[0]:
            a = mid[k]
            n = 5_000_000 if a < 0.02 else 400
            m = truth(cls, a)
            full.append({"class": cls, "bin_index": k + 1, "full_sum": m * n, "full_n": n})
            nj = int(n * (1 - 2 / NBLOCKS))
            for b in range(NBLOCKS):
                mj = m + (rng.normal(0, noise) if noise else 0.0)
                jk.append({"class": cls, "block": b, "bin_index": k + 1,
                           "jk_sum": mj * nj, "jk_n": nj})
    fp, jp = f"{out_dir}/full.tsv", f"{out_dir}/jk.tsv"
    pd.DataFrame(full).to_csv(fp, sep="\t", index=False)
    pd.DataFrame(jk).to_csv(jp, sep="\t", index=False)
    return fp, jp


def main():
    tmp = os.path.join(HERE, "tmp", "estimators")
    os.makedirs(tmp, exist_ok=True)
    bins = os.path.join(tmp, "bins_wide.txt")
    subprocess.run([sys.executable, os.path.join(SRC, "make_bins.py"), "--wide", "--out", bins],
                   check=True, capture_output=True)
    mid = E.load_grid(bins)

    S, N = E.load_arrays(*build(tmp, mid, 0.0), len(mid), NBLOCKS)
    got = E.jackknife(S, N, mid, NBLOCKS).set_index(["estimator", "classes"])["est"]

    expect = {
        ("h2_Unrel", "noPO"): H2,
        ("b2_FS", "noPO"): 2 * B_FS,
        ("b2_FS", "PO"): 2 * B_PO,
        ("h2_Pedf", "noPO"): H2,
        ("b2_step", "noPO"): 4 * (B_FS - B_HS),
        ("b2_step_ratio", "noPO"): B_FS / B_HS,
        ("excess", "FS"): B_FS,
        ("excess", "PO"): B_PO,
    }
    fails = 0
    for key, want in expect.items():
        ok = abs(got[key] - want) < 1e-9
        fails += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {key[0]:<14}{key[1]:<7} "
              f"got {got[key]:+.6f}  want {want:+.6f}")

    # biased by design; check the direction of the bias
    for key, cond, msg in (
        (("h2_FS", "noPO"), got[("h2_FS", "noPO")] > H2,
         "slope through the origin absorbs the FS intercept"),
        (("b2_FS", "pooled"), 2 * B_PO < got[("b2_FS", "pooled")] < 2 * B_FS,
         "PO pairs dilute the FS shared-environment estimate"),
    ):
        fails += not cond
        print(f"{'OK  ' if cond else 'FAIL'} {key[0]:<14}{key[1]:<7} "
              f"{got[key]:+.6f}  ({msg})")

    S, N = E.load_arrays(*build(tmp, mid, 0.002), len(mid), NBLOCKS)
    se = E.jackknife(S, N, mid, NBLOCKS).set_index(["estimator", "classes"])["se"]
    bad = [k for k in expect if not (np.isfinite(se[k]) and se[k] > 0)]
    fails += len(bad)
    print(f"\njackknife SEs finite and positive under noise: {len(expect) - len(bad)}/{len(expect)}")

    print("\nPASS -- estimators recover generating values" if not fails else f"\nFAIL ({fails})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
