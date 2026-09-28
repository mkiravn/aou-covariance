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

SURFACE, INK, INK2, MUTED, GRID, AXIS = ("#ffffff", "#0b0b0b", "#52514e",
                                         "#898781", "#e1e0d9", "#c3c2b7")

CLASS_COL = {"other": "#2a78d6", "FS": "#eb6834", "PO": "#1baf7a",
             "noPO": "#eb6834", "pooled": "#2a78d6"}
CLASS_MARK = {"other": "o", "FS": "s", "PO": "^", "noPO": "o", "pooled": "s"}
CLASS_LABEL = {"other": "unrelated / distant", "FS": "full sibs",
               "PO": "parent-offspring"}

MODEL_COL = {"h2_OneSlope": "#DD5129", "h2_Rel": "#0F7BA2",
             "h2_OneSlopeOffsets": "#43B284", "h2_RelOffsets": "#E89C33"}

RAMP = ["#86b6ef", "#4b93e8", "#2f76cf", "#215596"]

# diverging, for signed quantities: two hues either side of a neutral midpoint
DIVERGING_STEPS = ["#104281", "#f0efec", "#e34948"]

FONT_STACK = ["Helvetica Neue", "Helvetica", "Liberation Sans", "Nimbus Sans",
              "Arial", "DejaVu Sans"]

# y is the standardised residual written by 05_phenotype_residualize.ipynb
YLABEL = r"$\mathrm{E}[\tilde{y}_i\,\tilde{y}_j]$"

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
    for f in found:
        shutil.copy2(f, os.path.join(dst, os.path.basename(f)))
    print(f"{len(found)} PNG -> {dst}")
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
    per_stage = {}
    for s in found:
        stage = os.path.relpath(s, run_root).split(os.sep)[0]
        shutil.copy2(s, os.path.join(dst, f"{stage}__{os.path.basename(s)}"))
        per_stage[stage] = per_stage.get(stage, 0) + 1
    print(f"{len(found)} PNG -> {dst}")
    for stage, n in sorted(per_stage.items()):
        print(f"  {n:>4}  {stage}")
    if not found:
        print(f"nothing under {run_root} -- has any notebook been run for this run?")
    return dst, per_stage
