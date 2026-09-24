"""plink2 steps: QC, pruning, PCA, and projection scoring.

Each function returns the output prefix so calls chain readably in a notebook.
All of them are skip-if-exists, because these are minutes-to-hours long and
notebooks get rerun.
"""

import os

import numpy as np
import pandas as pd

from .env import sh

QC_FILTERS = "--maf 0.01 --hwe 1e-6 0 keep-fewhet --geno 0.05"


def _exists(prefix, ext=".pgen"):
    return os.path.isfile(prefix + ext)


def qc(pfile, out, keep=None, extra=QC_FILTERS, force=False):
    """Site + sample QC, producing a new pfile.

    `keep` restricts to a sample list first, so MAF/HWE are evaluated in the
    cohort actually being analysed rather than inherited from a wider one.
    """
    if _exists(out) and not force:
        print(f"qc exists: {out}")
        return out
    k = f'--keep "{keep}"' if keep else ""
    sh(f"""
plink2 --pfile "{pfile}" {k} --nonfounders {extra} \
  --max-alleles 2 --rm-dup exclude-all \
  --threads $(nproc) --make-pgen --out "{out}"
echo "{os.path.basename(out)}: $(wc -l < "{out}.psam") samples, $(grep -vc '^##' "{out}.pvar") variants"
""")
    return out


def prune(pfile, out, extract=None, ld_regions=None, params="1000kb 1 0.1",
          maf=None, force=False):
    """LD pruning, excluding long-range LD regions.

    `params` is plink's `--indep-pairwise` triple. r2=0.1 is the default here:
    at 0.05 the HM3 common set loses ~77% of its variants for little gain.
    """
    if os.path.isfile(f"{out}.prune.in") and not force:
        print(f"prune exists: {out}.prune.in "
              f"({sum(1 for _ in open(f'{out}.prune.in')):,} variants)")
        return out
    e = f'--extract "{extract}"' if extract else ""
    x = f'--exclude bed1 "{ld_regions}"' if ld_regions else ""
    m = f"--maf {maf}" if maf else ""
    sh(f"""
plink2 --pfile "{pfile}" --nonfounders {e} {x} {m} \
  --indep-pairwise {params} \
  --threads $(nproc) --out "{out}"
echo "after pruning: $(wc -l < "{out}.prune.in") variants"
""")
    return out


def pca(pfile, out, extract=None, n_pcs=20, approx=True, keep=None, force=False):
    """Fit a PCA, keeping allele weights.

    `allele-wts` is not optional for this pipeline: every downstream score --
    participants included -- goes through these loadings, because plink's
    `--pca` eigenvectors and a `--score` projection do not land on the same
    coordinates. Re-scoring everyone puts both sets on identical axes.
    """
    if os.path.isfile(f"{out}.eigenvec.allele") and not force:
        print(f"pca exists: {out}")
        return out
    e = f'--extract "{extract}"' if extract else ""
    k = f'--keep "{keep}"' if keep else ""
    a = "approx " if approx else ""
    sh(f"""
plink2 --pfile "{pfile}" --nonfounders {e} {k} --freq counts \
  --pca {a}{n_pcs} allele-wts \
  --threads $(nproc) --out "{out}"
echo "pca done: {out}.eigenvec.allele"
""")
    return out


def _score_cols(weights, n_pcs):
    """Column numbers plink's --score wants, read from the weights header."""
    h = open(weights).readline().split()
    idc = "#ID" if "#ID" in h else "ID"
    return (h.index(idc) + 1, h.index("A1") + 1,
            h.index("PC1") + 1, h.index(f"PC{n_pcs}") + 1)


def score(weights, out, pfile=None, bfile=None, freq=None, extract=None,
          n_pcs=20, force=False):
    """Project samples through a PCA's allele weights.

    Pass the fit's own `.acount` as `freq` so every cohort is centred on the
    same allele frequencies -- otherwise each projected set gets its own
    centring and the coordinates are not comparable.
    """
    if os.path.isfile(f"{out}.sscore") and not force:
        print(f"scores exist: {out}.sscore")
        return out
    assert (pfile is None) != (bfile is None), "pass exactly one of pfile/bfile"
    src = f'--pfile "{pfile}"' if pfile else f'--bfile "{bfile}"'
    e = f'--extract "{extract}"' if extract else ""
    f = f'--read-freq "{freq}"' if freq else ""
    i, a1, p1, pk = _score_cols(weights, n_pcs)
    sh(f"""
plink2 {src} --nonfounders {e} {f} \
  --score "{weights}" {i} {a1} header-read no-mean-imputation variance-standardize \
  --score-col-nums {p1}-{pk} \
  --threads $(nproc) --out "{out}"
echo "scored: $(($(wc -l < "{out}.sscore") - 1)) samples"
""")
    return out


def pc_names(n_pcs):
    return [f"PC{k}" for k in range(1, n_pcs + 1)]


def read_scores(sscore, id_col, n_pcs=20):
    """Read a .sscore into a tidy frame with columns [id_col, PC1..PCn]."""
    d = pd.read_csv(sscore, sep=r"\s+")
    idc = "#IID" if "#IID" in d.columns else "IID"
    d = d.rename(columns={idc: id_col,
                          **{f"PC{k}_AVG": f"PC{k}" for k in range(1, n_pcs + 1)}})
    d[id_col] = d[id_col].astype(str)
    return d[[id_col] + pc_names(n_pcs)]


def read_eigenval(prefix):
    """Eigenvalues as a percentage-of-variance series."""
    ev = np.loadtxt(f"{prefix}.eigenval")
    return ev, ev / ev.sum() * 100
