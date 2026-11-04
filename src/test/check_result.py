#!/usr/bin/env python3
"""Compare grm_class_tool merge output against the Python reference.

Checks:
  1. every (class, bin) in the reference appears with matching sum / sum_sq / n
  2. no extra (class, bin) rows beyond the reference
  3. the synthetic "pooled" class equals the sum of all real classes
  4. n_indiv (distinct people, capped at 21) matches for every class and pooled

Usage: check_result.py <expected.tsv> <merged.full.tsv>
"""
import sys

TOL = 1e-9


def read_expected(path):
    """(class, bin) -> (sum, sum_sq, n), and (class, bin) -> n_indiv. The
    reference writes pooled's people as `pooled_indiv` rows, since distinct
    people can't be summed over classes."""
    out, indiv = {}, {}
    with open(path) as f:
        next(f)
        for line in f:
            cls, k, s, sq, n, ni = line.rstrip("\n").split("\t")
            if cls == "pooled_indiv":
                indiv[("pooled", int(k))] = int(ni)
                continue
            out[(cls, int(k))] = (float(s), float(sq), int(n))
            indiv[(cls, int(k))] = int(ni)
    return out, indiv


def read_merged(path):
    out = {}
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        ci = {name: i for i, name in enumerate(header)}
        for line in f:
            f_ = line.rstrip("\n").split("\t")
            cls = f_[ci["class"]]
            k = int(f_[ci["bin_index"]]) - 1          # merge writes 1-based
            out[(cls, k)] = (float(f_[ci["full_sum"]]),
                             float(f_[ci["full_sum_sq"]]),
                             int(f_[ci["full_n"]]),
                             f_[ci["n_indiv"]])
    return out


def close(a, b):
    return abs(a - b) <= TOL * max(1.0, abs(a), abs(b))


def main():
    expected, exp_indiv = read_expected(sys.argv[1])
    merged = read_merged(sys.argv[2])

    fails = []

    for key, (es, esq, en) in sorted(expected.items()):
        if key not in merged:
            fails.append(f"MISSING {key}: expected n={en}")
            continue
        ms, msq, mn, _ = merged[key]
        if mn != en:
            fails.append(f"N MISMATCH {key}: expected {en}, got {mn}")
        if not close(ms, es):
            fails.append(f"SUM MISMATCH {key}: expected {es!r}, got {ms!r}")
        if not close(msq, esq):
            fails.append(f"SUMSQ MISMATCH {key}: expected {esq!r}, got {msq!r}")

    extra = {k for k in merged if k[0] != "pooled"} - set(expected)
    for key in sorted(extra):
        fails.append(f"UNEXPECTED {key}: n={merged[key][2]}")

    # pooled must equal the sum over real classes
    pooled_expected = {}
    for (cls, k), (s, sq, n) in expected.items():
        ps, psq, pn = pooled_expected.get(k, (0.0, 0.0, 0))
        pooled_expected[k] = (ps + s, psq + sq, pn + n)
    for k, (es, esq, en) in sorted(pooled_expected.items()):
        got = merged.get(("pooled", k))
        if got is None:
            fails.append(f"MISSING pooled bin {k}")
            continue
        if got[2] != en or not close(got[0], es):
            fails.append(f"POOLED MISMATCH bin {k}: expected (n={en}, sum={es!r}), got {got!r}")

    for key, ni in sorted(exp_indiv.items()):
        got = merged.get(key)
        if got is None:
            continue                      # already reported as MISSING
        if got[3] != str(ni):
            fails.append(f"N_INDIV MISMATCH {key}: expected {ni}, got {got[3]}")
    capped = sum(1 for v in exp_indiv.values() if v >= 21)

    classes = sorted({c for c, _ in expected})
    print(f"classes in reference : {classes}")
    print(f"reference rows       : {len(expected)}")
    print(f"merged rows (non-pooled): {len({k for k in merged if k[0] != 'pooled'})}")
    print(f"pairs by class       : "
          + ", ".join(f"{c}={sum(v[2] for (cc, _), v in expected.items() if cc == c)}"
                      for c in classes))
    print(f"n_indiv checked      : {len(exp_indiv)} (class, bin) cells, "
          f"{capped} at the cap, {len(exp_indiv) - capped} below it")
    print()

    if fails:
        print(f"FAIL -- {len(fails)} problem(s):")
        for msg in fails[:40]:
            print("  " + msg)
        sys.exit(1)
    print("PASS -- all class/bin sums, counts, people and the pooled total match the reference")


if __name__ == "__main__":
    main()
