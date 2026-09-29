"""Shared plot style, palettes and trait grouping for the notebooks.

Every notebook that draws anything used to carry its own copy of the rcParams
block, the palettes and a mirror_plots() helper -- five copies, which drifted
(the ramp was changed in one place and not another more than once). They live
here instead.

The palettes are not arbitrary. They were validated with the data-viz
palette checker, and the numbers below are the reason each is what it is:

  CLASS_COL   the three pair classes, as scatter marks, so the ALL-PAIRS
              colourblind floors apply: worst pair CVD dE 9.2, normal-vision
              24.0. A fourth hue alongside them does not clear those floors --
              every combination tried failed, worst cases green/orange dE 3.2
              under protanopia and red/orange dE 7.1 with normal vision.
  MODEL_COL   MetBrewer "Egypt", with its yellow darkened from #FAB255 to
              #E89C33 because the original sits above the lightness band on a
              light surface (L 0.815 against a 0.77 ceiling). These are LINE
              marks, so the adjacent pairlist applies: worst adjacent CVD
              dE 8.8, normal-vision 19.0. Aqua and yellow fall below 3:1
              against the surface, so the relief rule applies and every model
              keeps a dash pattern and a legend entry.
  RAMP        covariate sets are an ordered staircase, so a sequential ramp,
              not a categorical set. The dark end stops at #215596 rather than
              #104281: darker than that and a point was hard to tell from the
              black reference line. Passes the ordinal checks -- monotone
              lightness, adjacent dL >= 0.06, light end clear of the surface.
"""
import os

import numpy as np

SURFACE, INK, INK2, MUTED, GRID, AXIS = ("#ffffff", "#0b0b0b", "#52514e",
                                         "#898781", "#e1e0d9", "#c3c2b7")

CLASS_COL = {"other": "#2a78d6", "FS": "#eb6834", "PO": "#1baf7a",
             "noPO": "#eb6834", "pooled": "#2a78d6"}
CLASS_MARK = {"other": "o", "FS": "s", "PO": "^", "noPO": "o", "pooled": "s"}
CLASS_LABEL = {"other": "unrelated", "FS": "full sibs",
               "PO": "parent-offspring"}

MODEL_COL = {"h2_OneSlope": "#DD5129", "h2_Rel": "#0F7BA2",
             "h2_OneSlopeOffsets": "#43B284", "h2_RelOffsets": "#E89C33"}

RAMP = ["#86b6ef", "#4b93e8", "#2f76cf", "#215596"]

# diverging, for signed quantities: two hues either side of a neutral midpoint
DIVERGING_STEPS = ["#104281", "#f0efec", "#e34948"]

FONT_STACK = ["Helvetica Neue", "Helvetica", "Liberation Sans", "Nimbus Sans",
              "Arial", "DejaVu Sans"]

# y is the standardised residual written by 05_phenotype_residualize.ipynb
# Spelled out, not just the symbol: on a slide the axis is read before the
# caption, and "E[y~i y~j]" alone does not say what quantity it is.
YLABEL = r"Phenotypic cross-product ($\tilde{y}_i\,\tilde{y}_j$)"
XLABEL = r"Genetic relatedness ($a_{ij}$)"
# `other` is drawn in ink rather than its class hue: it is the dense bulk of
# the data, and the fits have to read ON TOP of it. Hue is then carried only by
# the two classes that are actually being contrasted, and by the fit lines.
POINT_COL = dict(CLASS_COL, other=INK)

CATEGORY = {
    **dict.fromkeys(["height", "weight", "bmi", "waist_circumference",
                     "hip_circumference", "waist_hip_ratio"], "anthropometric"),
    **dict.fromkeys(["systolic_bp", "diastolic_bp", "heart_rate"], "cardiovascular"),
    **dict.fromkeys(["glucose", "hba1c", "hdl_cholesterol", "ldl_cholesterol",
                     "triglycerides", "total_cholesterol"], "metabolic"),
    **dict.fromkeys(["hemoglobin", "wbc_count", "platelet_count", "rbc_count",
                     "hematocrit", "mch", "mcv", "mchc", "neutrophil_count",
                     "lymphocyte_count", "monocyte_count", "basophil_count",
                     "eosinophil_count"], "blood cells"),
    **dict.fromkeys(["creatinine", "alt"], "liver & kidney"),
    **dict.fromkeys(["alcohol_audit_c_score", "cigarettes_per_day", "pss_score",
                     "social_support_score", "loneliness_score", "eds_score"],
                    "behavioural"),
}
CAT_ORDER = ["anthropometric", "cardiovascular", "metabolic", "blood cells",
             "liver & kidney", "behavioural", "uncategorised"]
CAT_COL = dict(zip(CAT_ORDER, ["#2a78d6", "#e34948", "#008300", "#4a3aa7",
                               "#8c564b", "#eb6834", "#898781"]))


def resolve_font():
    """The first FONT_STACK entry the machine actually has, else DejaVu Sans.

    Passing the whole stack to rcParams makes matplotlib warn about every
    missing family on every text draw -- tens of lines per figure on a Linux
    image with none of them.
    """
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    return next((f for f in FONT_STACK if f in have), "DejaVu Sans")


def apply_style(transparent=False, verbose=True):
    """Set rcParams for every figure in the pipeline. Returns the font used.

    theme_bw-ish: a full panel border, major gridlines only and few of them,
    white (or transparent) background.
    """
    import matplotlib.pyplot as plt
    font = resolve_font()
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "savefig.transparent": transparent,
        "figure.dpi": 120,
        "font.family": font, "font.size": 9,
        "mathtext.fontset": "custom", "mathtext.rm": font,
        "mathtext.it": f"{font}:italic", "mathtext.bf": f"{font}:bold",
        "axes.edgecolor": INK2, "axes.linewidth": 0.7, "axes.labelcolor": INK2,
        "axes.titlecolor": INK, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.grid": True, "axes.grid.which": "major", "axes.axisbelow": True,
        "grid.color": GRID, "grid.linewidth": 0.5, "grid.alpha": 0.7,
        "xtick.minor.visible": False, "ytick.minor.visible": False,
        "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "legend.frameon": False,
    })
    if verbose:
        have = {f for f in FONT_STACK if f == font}
        print(f"font in use: {font}" + ("" if have else "  (fallback)"))
    return font


def diverging_cmap(name="cov_diverging"):
    """Two hues with a neutral midpoint, for signed quantities like z."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, DIVERGING_STEPS)


def mirror_plots(src_dir, tag, run, root="~/plots"):
    """Copy every PNG under src_dir into <root>/<run>/<tag>/.

    The bucket mount is the durable copy; this is a flat local directory so a
    whole notebook's figures can be flicked through without walking the bucket.
    """
    import glob
    import shutil
    dst = os.path.join(os.path.expanduser(root), run, tag)
    os.makedirs(dst, exist_ok=True)
    found = sorted(glob.glob(os.path.join(src_dir, "**", "*.png"), recursive=True))
    keep = {os.path.basename(f) for f in found}
    # Prune first. The bucket accumulates figures written under filenames that
    # earlier versions of a notebook used, and a stale one in here is
    # indistinguishable from a current one -- the whole point of this directory
    # is to be flicked through without checking provenance. Only this mirror is
    # pruned; the bucket copy is the durable record and is left alone.
    stale = [f for f in glob.glob(os.path.join(dst, "*.png"))
             if os.path.basename(f) not in keep]
    for f in stale:
        os.remove(f)
    for f in found:
        shutil.copy2(f, os.path.join(dst, os.path.basename(f)))
    print(f"{len(found)} PNG -> {dst}"
          + (f"  ({len(stale)} stale removed)" if stale else ""))
    return dst


def collect_run_plots(run_root, run, root="~/plots"):
    """Every PNG anywhere under run_root, into <root>/<run>/all/.

    Prefixed with the stage directory: several stages reuse filenames, and the
    per-notebook mirrors only cover notebooks that have been re-run.
    """
    import glob
    import shutil
    dst = os.path.join(os.path.expanduser(root), run, "all")
    os.makedirs(dst, exist_ok=True)
    found = sorted(glob.glob(os.path.join(run_root, "**", "*.png"), recursive=True))
    per_stage, keep = {}, set()
    for s in found:
        stage = os.path.relpath(s, run_root).split(os.sep)[0]
        name = f"{stage}__{os.path.basename(s)}"
        keep.add(name)
        shutil.copy2(s, os.path.join(dst, name))
        per_stage[stage] = per_stage.get(stage, 0) + 1
    # same pruning as mirror_plots: a figure from an earlier version of a
    # notebook, under a filename it no longer writes, is worse than missing
    stale = [f for f in glob.glob(os.path.join(dst, "*.png"))
             if os.path.basename(f) not in keep]
    for f in stale:
        os.remove(f)
    print(f"{len(found)} PNG -> {dst}"
          + (f"  ({len(stale)} stale removed)" if stale else ""))
    for stage, n in sorted(per_stage.items()):
        print(f"  {n:>4}  {stage}")
    if not found:
        print(f"nothing under {run_root} -- has any notebook been run for this run?")
    return dst, per_stage


def si(n):
    """Pair counts for a legend key: 1996971 -> "2.0M". Spelled-out counts were
    the widest thing in the legend and their last five digits carry nothing."""
    n = float(n)
    for cut, suf in ((1e9, "B"), (1e6, "M"), (1e3, "k")):
        if abs(n) >= cut:
            return f"{n / cut:.1f}{suf}".replace(".0", "")
    return f"{n:,.0f}"


def bulk_limits(centres, los, his, q=75, pad=0.06, floor=None):
    """Axis limits driven by the bulk of the intervals, not the widest one.

    A coefficient fitted on a handful of bins can be an order of magnitude less
    precise than the rest of a forest panel, and letting its interval set the
    axis squashes every other row into a sliver around zero. These limits cover
    every point estimate plus the `q`-th percentile of the half-widths, so most
    intervals fit and only the outliers need capping.

    `floor` is a range always included -- (0, 1) for a heritability panel, so
    the reader always has that scale even when every estimate is small.
    Returns (lo, hi), or None if there is nothing finite to bound.
    """
    c = np.asarray(centres, float)
    half = np.fmax(np.asarray(his, float) - c, c - np.asarray(los, float))
    c, half = c[np.isfinite(c)], half[np.isfinite(half)]
    if not len(c):
        return None
    reach = np.percentile(half, q) if len(half) else 0.0
    lo, hi = c.min() - reach, c.max() + reach
    if floor is not None:
        lo, hi = min(lo, floor[0]), max(hi, floor[1])
    span = hi - lo or max(abs(hi), 1.0)
    return lo - pad * span, hi + pad * span


def capped_errorbar(ax, x, y, lo, hi, horizontal=True, cap_ms=6.0, **kw):
    """`errorbar` that draws an arrowhead where an interval leaves the axis.

    Call after the limits are set. An interval clipped at the axis edge without
    a marker reads as a narrower interval than it is, which is the one thing an
    interval must not do.

    The arrowheads are drawn by hand rather than with matplotlib's
    `xlolims`/`xuplims`: those put the arrow *at the data point* and drop the
    bar entirely, so an interval clipped on both sides renders as a bare marker
    with no interval at all -- the opposite of what it needs to say.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    a0, a1 = ax.get_xlim() if horizontal else ax.get_ylim()
    c = x if horizontal else y          # the centre on the error axis
    have = np.isfinite(lo) & np.isfinite(hi)
    under, over = have & (lo < a0), have & (hi > a1)
    # np.fmax/fmin IGNORE NaN, so clipping a missing interval against the axis
    # would silently produce a full-width bar. Zero those rows explicitly.
    err = np.array([c - np.fmax(lo, a0), np.fmin(hi, a1) - c])
    err[:, ~have] = 0.0
    err[~np.isfinite(err)] = 0.0
    np.clip(err, 0.0, None, out=err)    # a centre outside the axis can flip a side
    art = ax.errorbar(x, y, **{"xerr" if horizontal else "yerr": err}, **kw)

    col = kw.get("color") or kw.get("ecolor") or INK
    z = kw.get("zorder", 3)
    for flag, mark, edge in ((under, "<" if horizontal else "v", a0),
                             (over, ">" if horizontal else "^", a1)):
        if not np.any(flag):
            continue
        px = np.full(int(flag.sum()), edge) if horizontal else x[flag]
        py = y[flag] if horizontal else np.full(int(flag.sum()), edge)
        ax.plot(px, py, mark, ms=cap_ms, color=col, mec="none", ls="none",
                clip_on=False, zorder=z + 0.1)
    return art


def shades(hex_colour, k, lo=0.45, hi=1.0):
    """`k` lightness steps of one hue, darkest first.

    For a model whose line is drawn in pieces that mean different things -- the
    per-band slopes -- so the pieces stay one model to the eye while still
    being told apart. A second hue would claim to be a second model.
    """
    import matplotlib.colors as mcolors
    r, g, b = mcolors.to_rgb(hex_colour)
    out = []
    for t in (np.linspace(lo, hi, k) if k > 1 else [hi]):
        # t < 1 darkens toward black, t = 1 is the colour itself
        out.append(mcolors.to_hex((r * t, g * t, b * t)))
    return out


def repel_labels(ax, x, y, labels, fontsize=6.0, color=None, iters=400,
                 link_lw=0.35):
    """Scatter labels nudged apart, with a leader line back to the point.

    matplotlib has no `geom_text_repel`, and on the phenotype scatters the
    labels sat on top of each other and on the markers. This is the usual
    iterative scheme, with one thing that matters: labels repel as BOXES sized
    from the text, not as discs. A disc of the text's height leaves wide labels
    overlapping horizontally; a disc of its width pushes short ones absurdly
    far apart. Pairs are separated along whichever axis they overlap least on,
    which is what keeps a row of labels from being flung vertically.

    Works in axes fraction internally, so the data units cannot bias which
    direction things move.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    labels = [str(l) for l, k in zip(labels, ok) if k]
    x, y = x[ok], y[ok]
    if not len(x):
        return []
    color = color or MUTED
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    sx, sy = (x1 - x0) or 1.0, (y1 - y0) or 1.0

    # label half-extents in axes fraction, from the axes size in points
    bb = ax.get_window_extent()
    ax_w = max(bb.width, 1.0) * 72.0 / ax.figure.dpi
    ax_h = max(bb.height, 1.0) * 72.0 / ax.figure.dpi
    hw = np.array([0.5 * len(l) * 0.58 * fontsize / ax_w for l in labels])
    hh = np.full(len(labels), 0.5 * 1.35 * fontsize / ax_h)
    mk = 0.5 * 5.0 / ax_h                       # marker half-size, points -> frac

    px, py = (x - x0) / sx, (y - y0) / sy       # points
    lx, ly = px.copy(), py + hh + mk            # labels start just above
    rng = np.random.default_rng(0)
    lx += rng.uniform(-1e-4, 1e-4, len(lx))     # break exact ties

    for _ in range(iters):
        moved = False
        for i in range(len(lx)):
            ddx, ddy = 0.0, 0.0
            # label vs label, as boxes
            ox = (hw + hw[i]) - np.abs(lx - lx[i])
            oy = (hh + hh[i]) - np.abs(ly - ly[i])
            ox[i] = oy[i] = -1.0
            hit = (ox > 0) & (oy > 0)
            for j in np.where(hit)[0]:
                if ox[j] < oy[j]:               # cheaper to separate on x
                    ddx += np.sign(lx[i] - lx[j] or 1e-6) * ox[j]
                else:
                    ddy += np.sign(ly[i] - ly[j] or 1e-6) * oy[j]
            # label vs every marker
            mox = (hw[i] + mk) - np.abs(lx[i] - px)
            moy = (hh[i] + mk) - np.abs(ly[i] - py)
            for j in np.where((mox > 0) & (moy > 0))[0]:
                ddy += np.sign(ly[i] - py[j] or 1.0) * moy[j]
            if ddx or ddy:
                moved = True
                lx[i] += 0.5 * ddx
                ly[i] += 0.5 * ddy
        if not moved:
            break
    lx = np.clip(lx, hw, 1 - hw)
    ly = np.clip(ly, hh, 1 - hh)

    arts = []
    for i, lab in enumerate(labels):
        if np.hypot(lx[i] - px[i], ly[i] - py[i]) > 2 * hh[i]:
            ax.plot([x0 + px[i] * sx, x0 + lx[i] * sx],
                    [y0 + py[i] * sy, y0 + ly[i] * sy],
                    lw=link_lw, color=color, alpha=0.55, zorder=2)
        arts.append(ax.text(x0 + lx[i] * sx, y0 + ly[i] * sy, lab,
                            fontsize=fontsize, color=color, ha="center",
                            va="center", zorder=6))
    return arts
