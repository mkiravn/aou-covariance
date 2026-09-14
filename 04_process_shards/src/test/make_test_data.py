#!/usr/bin/env python3
"""Generate a small GRM + phenotype + pair-class set, split it into plink-style
shards, and compute the expected per-class bin sums directly from the dense
matrix.

Writes into <out_dir>:
  test.grm.bin, test.grm.id, test.pheno, test.bins, test.classes.tsv,
  shard.<k>            (k = 1..n_shards)
  expected.tsv         (class, bin_index, sum, sum_sq, n) -- the reference

Usage: make_test_data.py <out_dir> [n_ids] [n_shards]
"""
import math
import os
import sys

import numpy as np

N_DEFAULT = 60
SHARDS_DEFAULT = 5


def triangle_divide_off_diag(target: int) -> int:
    if target == 0:
        return 1
    v = int(math.sqrt(target)) + 2
    while v > 1 and (v - 1) * (v - 2) >= target:
        v -= 1
    while v * (v - 1) < target:
        v += 1
    return v


def plink_parallel_row_start(n_ids: int, parallel_idx: int, n_shards: int) -> int:
    target = (n_ids * (n_ids - 1) * parallel_idx) // n_shards
    v = triangle_divide_off_diag(target)
    return 0 if v == 1 else v


def get_bin(g: float, bins) -> int:
    for k, (lo, hi) in enumerate(bins):
        last = k == len(bins) - 1
        if (lo <= g <= hi) if last else (lo <= g < hi):
            return k
    return -1


def main() -> None:
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else N_DEFAULT
    n_shards = int(sys.argv[3]) if len(sys.argv) > 3 else SHARDS_DEFAULT
    os.makedirs(out, exist_ok=True)

    rng = np.random.default_rng(0)

    # symmetric relatedness: unrelated noise plus injected related pairs
    A = rng.normal(0.0, 0.01, size=(n, n))
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 1.0)

    related = [
        (5, 3, 0.50), (10, 7, 0.52), (20, 11, 0.48), (40, 2, 0.51),   # 1st degree
        (30, 29, 1.00), (55, 4, 0.99),                                 # dup/MZ
        (25, 13, 0.25), (44, 31, 0.26),                                # 2nd degree
    ]
    for i, j, v in related:
        A[i, j] = v
        A[j, i] = v

    ids = [f"ID{idx:04d}" for idx in range(n)]

    # dense lower triangle, float32: row i holds j = 0..i (last entry is diagonal)
    with open(f"{out}/test.grm.bin", "wb") as f:
        for i in range(n):
            f.write(A[i, : i + 1].astype(np.float32).tobytes())

    # .grm.id with plink's FID = 0 convention
    with open(f"{out}/test.grm.id", "w") as f:
        for iid in ids:
            f.write(f"0\t{iid}\n")

    # phenotype: FID = IID = the id (deliberately NOT matching .grm.id's FID,
    # to confirm the tool keys on IID alone). A few are missing.
    y = rng.normal(0.0, 1.0, size=n)
    missing = {2, 17, 29}
    with open(f"{out}/test.pheno", "w") as f:
        f.write("FID IID Y\n")
        for idx, iid in enumerate(ids):
            val = "NA" if idx in missing else f"{y[idx]:.10f}"
            f.write(f"{iid} {iid} {val}\n")
    y_use = np.array([np.nan if i in missing else y[i] for i in range(n)])

    bins = [(-0.05, 0.05), (0.05, 0.20), (0.20, 0.35), (0.35, 0.70), (0.70, 1.50)]
    with open(f"{out}/test.bins", "w") as f:
        for lo, hi in bins:
            f.write(f"{lo} {hi}\n")

    # pair classes: some related pairs labelled, one unmappable, one duplicate row
    classes = [
        (ids[5], ids[3], "PO"),
        (ids[10], ids[7], "FS"),
        (ids[7], ids[10], "FS"),        # duplicate, reversed order -- must be ignored
        (ids[20], ids[11], "PO"),
        (ids[40], ids[2], "FS"),        # ID2 has a missing phenotype -> never binned
        (ids[30], ids[29], "DUP"),      # ID29 missing phenotype -> never binned
        (ids[55], ids[4], "DUP"),
        ("NOT_AN_ID", ids[0], "PO"),    # unmappable
    ]
    with open(f"{out}/test.classes.tsv", "w") as f:
        f.write("IID1\tIID2\tcls\n")
        for a, b, c in classes:
            f.write(f"{a}\t{b}\t{c}\n")

    class_of = {}
    for a, b, c in classes:
        if a not in ids or b not in ids:
            continue
        i, j = ids.index(a), ids.index(b)
        key = (max(i, j), min(i, j))
        class_of.setdefault(key, c)      # first wins, matching the tool

    # shards
    starts = [plink_parallel_row_start(n, k, n_shards) for k in range(n_shards)] + [n]
    with open(f"{out}/test.grm.bin", "rb") as src:
        data = src.read()
    for k in range(n_shards):
        a, b = starts[k], starts[k + 1]
        off = a * (a + 1) // 2 * 4
        end = b * (b + 1) // 2 * 4
        with open(f"{out}/shard.{k + 1}", "wb") as f:
            f.write(data[off:end])

    # reference: iterate the whole lower triangle exactly as the tool does
    acc = {}
    for i in range(n):
        for j in range(i):
            yi, yj = y_use[i], y_use[j]
            if np.isnan(yi) or np.isnan(yj):
                continue
            g = np.float32(A[i, j])
            k = get_bin(float(g), bins)
            if k < 0:
                continue
            cls = class_of.get((i, j), "other")
            prod = float(yi) * float(yj)
            s, sq, cnt = acc.get((cls, k), (0.0, 0.0, 0))
            acc[(cls, k)] = (s + prod, sq + prod * prod, cnt + 1)

    with open(f"{out}/expected.tsv", "w") as f:
        f.write("class\tbin_index\tsum\tsum_sq\tn\n")
        for (cls, k), (s, sq, cnt) in sorted(acc.items()):
            f.write(f"{cls}\t{k}\t{s:.17g}\t{sq:.17g}\t{cnt}\n")

    print(f"n_ids={n} n_shards={n_shards} bins={len(bins)}")
    print(f"shard row starts: {starts}")
    print(f"classed pairs mapped: {len(class_of)}")
    print(f"reference rows: {len(acc)}")


if __name__ == "__main__":
    main()
