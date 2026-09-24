"""1000 Genomes reference data and the long-range LD exclusion lists."""

import os

import pandas as pd

from .config import WS, WS_GS
from .env import sh

KG_DIR = f"{WS}/1000g_reference"
KG_BFILE = f"{KG_DIR}/1kg_all_qc"
KG_ACOUNT = f"{KG_DIR}/1kg_all_qc.acount"
KG_PANEL = f"{KG_DIR}/integrated_call_samples_v3.20130502.ALL.panel"

KG_BFILE_GS = f"{WS_GS}/1000g_reference/1kg_all_qc"
KG_ACOUNT_GS = f"{WS_GS}/1000g_reference/1kg_all_qc.acount"

EUR_POPS = ["CEU", "GBR", "FIN", "TSI", "IBS"]
ANCHOR_POPS = ["CEU", "GBR"]  # the Kemper et al. anchor

POP_COLORS = {"CEU": "#e6194b", "GBR": "#f58231", "FIN": "#3cb44b",
              "TSI": "#4363d8", "IBS": "#911eb4"}

# Canonical long-range LD / inversion regions (GRCh38), 1 Mb buffer added around
# each; small singleton intervals from the original list are subsumed by their
# neighbours.
HIGH_LD_CANONICAL = """\
chr1 47761740 51761740 1
chr2 85919365 100517106 2
chr2 182427027 189427029 3
chr3 47483505 49987563 4
chr3 83368158 86868160 5
chr5 44464140 51168409 6
chr5 129636407 132636409 7
chr6 25391792 33424245 8
chr6 57788603 58453888 9
chr6 61109122 61357029 10
chr6 139637169 142137170 11
chr7 54964812 66897578 12
chr8 8105067 12105082 13
chr8 43025699 48924888 14
chr8 110918594 113918595 15
chr10 36671065 43184546 16
chr11 88127183 91127184 17
chr12 32955798 41319931 18
chr20 33948532 36438183 19
"""

# The canonical list plus the narrow peaks found by `checks.ld_peaks` on the
# round-2 gate PCA. Prefer this one: the extra intervals are cheap and the peaks
# were real (single-region PCs in the first pass).
HIGH_LD_WITH_PEAKS = HIGH_LD_CANONICAL + """\
chr1 125169943 125170022 20
chr1 144106678 144106709 21
chr1 181955019 181955047 22
chr2 87416141 87416186 23
chr2 87417804 87417863 24
chr2 87418924 87418981 25
chr2 89917298 89917322 26
chr2 135275091 135275210 27
chr2 207609786 207609808 28
chr6 26726947 26726981 29
chr6 61424410 61424451 30
chr7 62182500 62277073 31
chr8 47303500 47317337 32
chr9 40365644 40365693 33
chr9 64198500 64200392 34
chr9 88958735 88959017 35
chr10 41693521 41885273 36
chr12 34639034 34639084 37
chr14 87391719 87391996 38
chr14 94658026 94658080 39
chr17 43159541 43159574 40
chr20 4031884 4032441 41
chr22 30060084 30060162 42
chr22 42980497 42980522 43
"""


def write_ld_regions(path, regions=HIGH_LD_WITH_PEAKS):
    """Write an exclusion list for `plink2 --exclude bed1`."""
    with open(path, "w") as fh:
        fh.write(regions)
    return path


def panel():
    """1000G sample -> population / superpopulation."""
    return pd.read_csv(KG_PANEL, sep=r"\s+")[["sample", "pop", "super_pop"]]


def eur_keep(path):
    """Write the 1000G European sample IDs, for a Europeans-only PCA fit."""
    p = panel()
    p.loc[p["super_pop"].eq("EUR"), "sample"].to_csv(path, index=False, header=False)
    return path


def shared_variants(pvar_or_alleles, out_ids, kg_acount=KG_ACOUNT, cols="3,4,5"):
    """Variants matching 1000G on ID *and* both alleles.

    Matching on ID alone silently keeps strand-ambiguous and re-mapped sites
    that then project into nonsense. `cols` selects the ID/REF/ALT columns of
    the participant-side file (a .pvar is 3,4,5; an .eigenvec.allele is 2,4,5).
    """
    c = [f"${x}" for x in cols.split(",")]
    sh(f"""
grep -v '^##' "{pvar_or_alleles}" | awk 'NR>1 {{print {", ".join(c)}}}' \
  | LC_ALL=C sort > "{out_ids}.lhs"
awk 'NR>1 {{print $2, $3, $4}}' "{kg_acount}" \
  | LC_ALL=C sort > "{out_ids}.kg"
LC_ALL=C comm -12 "{out_ids}.lhs" "{out_ids}.kg" | awk '{{print $1}}' > "{out_ids}"
echo "shared with 1000G: $(wc -l < "{out_ids}") variants"
""")
    return out_ids


def assert_inputs():
    """Fail early and loudly if the reference data is not mounted."""
    for p in (f"{KG_BFILE}.bed", KG_ACOUNT, KG_PANEL):
        assert os.path.isfile(p), f"missing reference file: {p}"
