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
        return m + np.where((a >= 0.2) & (a < 0.36), B_HS, 0.0)
    return m + (B_FS if cls == "FS" else B_PO)


def present(cls, a):
    # "other" is absent from 0.18-0.2 and 0.36-0.4 so the HS offset fills both
    # the b2_step band (0.2-0.4) and the deg2 band (0.18-0.36) exactly
    related = (a >= 0.4) & (a < 0.6)
    if cls == "other":
        return (a < 0.18) | ((a >= 0.2) & (a < 0.36)) | (a >= 0.7)
    return related


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
        ("h2_UnrelInt", "noPO"): H2,
        ("Unrel.intercept", "noPO"): 0.0,
        ("b2_FS", "noPO"): 2 * B_FS,
        ("b2_FS", "PO"): 2 * B_PO,
        ("b2_step", "noPO"): 4 * (B_FS - B_HS),
        ("b2_step_ratio", "noPO"): B_FS / B_HS,
        ("excess", "FS"): B_FS,
        ("excess", "PO"): B_PO,
        ("h2_OneSlopeOffsets", "noPO"): H2,
        ("h2_RelOffsets", "noPO"): H2,
        ("diff_RelOffsets-Unrel", "noPO"): 0.0,
        # `pooled` is no longer computed, so the checks that used it -- that
        # pooling averages the FS and PO offsets and leaves the slope untouched
        # -- are gone with it. The per-class machinery is still covered by the
        # noPO and PO expectations above.
    }
    # Band names and ranges come from the module, so re-cutting the bands
    # cannot leave this test asserting on one that no longer exists. The truth
    # is keyed by RANGE: the generating intercept is B_HS over the half-sib
    # band and B_FS over the first-degree band, whatever those bands are called.
    TRUE_BY_RANGE = {(0.09, 0.18): 0.0, (0.18, 0.36): B_HS, (0.36, E.R_HI): B_FS}
    expect.update({(f"OneSlopeOffsets.int.{d}", "noPO"): TRUE_BY_RANGE[r]
                   for d, r in E.OFFSET_BANDS.items()})
    expect.update({(f"RelOffsets.int.{d}", "noPO"): TRUE_BY_RANGE[r]
                   for d, r in E.OFFSET_BANDS.items()})
    # a slope per band, each recovering the generating slope
    expect.update({(f"RelOffsets.slope.{d}", "noPO"): H2 for d in E.OFFSET_BANDS})
    fails = 0
    for key, want in expect.items():
        ok = abs(got[key] - want) < 1e-9
        fails += not ok
        print(f"{'OK  ' if ok else 'FAIL'} {key[0]:<24}{key[1]:<7} "
              f"got {got[key]:+.6f}  want {want:+.6f}")

    # biased by design; check the direction of the bias
    for key, cond, msg in (
        (("h2_FS", "noPO"), got[("h2_FS", "noPO")] > H2,
         "slope through the origin absorbs the FS intercept"),
        (("h2_Rel", "noPO"), got[("h2_Rel", "noPO")] > H2,
         "one intercept for the whole related range cannot absorb both the HS "
         "and FS levels, so the slope takes up the slack"),
        (("diff_Rel-RelOffsets", "noPO"), got[("diff_Rel-RelOffsets", "noPO")] > 0,
         "offsets remove that excess"),
        (("h2_OneSlope", "noPO"), H2 < got[("h2_OneSlope", "noPO")],
         "a single origin slope over everything is pulled up by the related "
         "intercepts"),
    ):
        fails += not cond
        print(f"{'OK  ' if cond else 'FAIL'} {key[0]:<24}{key[1]:<7} "
              f"{got[key]:+.6f}  ({msg})")

    S, N = E.load_arrays(*build(tmp, mid, 0.002), len(mid), NBLOCKS)
    noisy = E.jackknife(S, N, mid, NBLOCKS).set_index(["estimator", "classes"])
    se = noisy["se"]
    bad = [k for k in expect if not (np.isfinite(se[k]) and se[k] > 0)]
    fails += len(bad)
    print(f"\njackknife SEs finite and positive under noise: {len(expect) - len(bad)}/{len(expect)}")

    # t critical values, against R's qt(p, df). The module computes these from
    # its own incomplete beta rather than scipy, so they are worth pinning.
    print()
    for df, level, want in ((NBLOCKS - 1, 0.95, 2.009575), (NBLOCKS - 1, 0.99, 2.679952),
                            (19, 0.95, 2.093024), (1, 0.95, 12.706205),
                            (99, 0.95, 1.984217)):
        got = E.t_crit(df, level)
        ok = abs(got - want) < 1e-5
        fails += not ok
        print(f"{'OK  ' if ok else 'FAIL'} t_crit(df={df}, {level}) = {got:.6f}  "
              f"R qt = {want:.6f}")
    # and the interval is that multiplier applied to the SE, both-sided
    tc = E.t_crit(NBLOCKS - 1, E.CI_LEVEL)
    width = ((noisy["ci_hi"] - noisy["ci_lo"]) / noisy["se"]).dropna()
    brackets = ((noisy["ci_lo"] <= noisy["est"]) & (noisy["est"] <= noisy["ci_hi"]))
    for name, cond, msg in (
        ("CI width is 2*t*se", np.allclose(width, 2 * tc), f"expected {2 * tc:.6f} SEs wide"),
        ("CI brackets est", bool(brackets[noisy["se"].notna()].all()),
         "every interval contains its own point estimate"),
        ("CI df is nblocks-1", bool((noisy["ci_df"] == NBLOCKS - 1).all()),
         "the SE comes from nblocks replicates"),
    ):
        fails += not cond
        print(f"{'OK  ' if cond else 'FAIL'} {name:<24}({msg})")

    # the data carry band offsets, so the offset models must fit and predict better
    fit, contrasts = E.compare_models(S, N, mid, NBLOCKS)
    fit = fit.set_index("model")
    contrasts = contrasts.set_index("contrast")
    print("\n" + fit.round(4).to_string())
    print(contrasts.round(3).to_string())
    for name, cond, msg in (
        ("chi2 favours offsets",
         fit.loc["h2_RelOffsets", "chi2_per_bin"] < fit.loc["Rel_banded", "chi2_per_bin"],
         "adding band offsets lowers misfit, on the same bins"),
        ("cv no worse",
         fit.loc["h2_RelOffsets", "cv_error"] <= fit.loc["h2_OneSlope", "cv_error"] * 1.01,
         "held-out error is small here: each replicate leaves out 2/nblocks of the pairs"),
        ("band intercepts earn their place",
         contrasts.loc["band intercepts", "delta_chi2_per_bin"] > 0,
         "the data were generated with a per-band level shift"),
        ("band slopes add nothing",
         abs(contrasts.loc["band slopes too", "delta_chi2_per_bin"]) < 1e-6,
         "one slope generated every band, so freeing the slopes cannot help"),
    ):
        fails += not cond
        print(f"{'OK  ' if cond else 'FAIL'} {name:<24}({msg})")

    print("\nPASS -- estimators recover generating values" if not fails else f"\nFAIL ({fails})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
