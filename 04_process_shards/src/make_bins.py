#!/usr/bin/env python3
"""Generate a relatedness bins file for grm_class_tool / grm_shard_tool.

Three regions, each with its own width:

  fine    -0.05 .. 0.02   the unrelated region, where h2_Unrel's slope is fit
  mid      0.02 .. 0.20   distant relatedness
  coarse   0.20 .. 1.50   2nd degree (0.25), 1st degree (0.5), dup/MZ (1.0)

Defaults reproduce GRM-pairs/full_grm_bin/bins.txt exactly (236 bins).
--wide doubles the mid and coarse widths, halving the bin count above 0.02 so
each bin holds roughly twice as many pairs. That matters once pairs are split
by class: the ~8k first-degree pairs spread over ~30 bins at width 0.01, and
splitting them into PO and FS leaves few per bin.

Widths are handled in integer thousandths so the output has no float noise.

Usage:
    make_bins.py --out bins.txt                     # 236 bins (default)
    make_bins.py --out bins_wide.txt --wide         # 153 bins
    make_bins.py --out b.txt --mid-width 0.01 --coarse-width 0.02
"""
import argparse
import sys

FINE_LO, FINE_HI = -50, 20        # thousandths
MID_HI = 200
COARSE_HI = 1500

DEFAULTS = {"fine": 1, "mid": 5, "coarse": 10}


def build_edges(fine_w: int, mid_w: int, coarse_w: int) -> list:
    for name, w, lo, hi in (("fine", fine_w, FINE_LO, FINE_HI),
                            ("mid", mid_w, FINE_HI, MID_HI),
                            ("coarse", coarse_w, MID_HI, COARSE_HI)):
        if w <= 0:
            sys.exit(f"{name} width must be positive")
        if (hi - lo) % w:
            sys.exit(f"{name} width {w/1000:g} does not divide "
                     f"[{lo/1000:g}, {hi/1000:g}] evenly")

    edges = list(range(FINE_LO, FINE_HI, fine_w))
    edges += list(range(FINE_HI, MID_HI, mid_w))
    edges += list(range(MID_HI, COARSE_HI + 1, coarse_w))
    return edges


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--wide", action="store_true",
                   help="double the mid and coarse widths")
    p.add_argument("--fine-width", type=float)
    p.add_argument("--mid-width", type=float)
    p.add_argument("--coarse-width", type=float)
    a = p.parse_args()

    mult = 2 if a.wide else 1
    fine_w = round(a.fine_width * 1000) if a.fine_width else DEFAULTS["fine"]
    mid_w = round(a.mid_width * 1000) if a.mid_width else DEFAULTS["mid"] * mult
    coarse_w = round(a.coarse_width * 1000) if a.coarse_width else DEFAULTS["coarse"] * mult

    edges = build_edges(fine_w, mid_w, coarse_w)

    with open(a.out, "w") as f:
        for lo, hi in zip(edges[:-1], edges[1:]):
            f.write(f"{lo/1000:.3f} {hi/1000:.3f}\n")

    n_fine = (FINE_HI - FINE_LO) // fine_w
    n_mid = (MID_HI - FINE_HI) // mid_w
    n_coarse = (COARSE_HI - MID_HI) // coarse_w
    print(f"wrote {a.out}: {len(edges) - 1} bins")
    print(f"  {fine_w/1000:<6.3f} x {n_fine:>4}   {FINE_LO/1000:g} to {FINE_HI/1000:g}")
    print(f"  {mid_w/1000:<6.3f} x {n_mid:>4}   {FINE_HI/1000:g} to {MID_HI/1000:g}")
    print(f"  {coarse_w/1000:<6.3f} x {n_coarse:>4}   {MID_HI/1000:g} to {COARSE_HI/1000:g}")


if __name__ == "__main__":
    main()
