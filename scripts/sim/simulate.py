#!/usr/bin/env python3
"""Simulated AoU-shaped inputs for running the pipeline locally. Not participant data.

Writes, under --out (default .sim/inputs), the files a real run reads from the CDR:
  acaf/acaf_threshold.chr{c}.{pgen,pvar,psam} (IDs chr:pos:REF, multi-allelic sites on one row, as AoU ships them),
  aux/ancestry/{ancestry_preds,training_pca}.tsv,
  aux/relatedness/samples_relatedness.tsv, reference/1kg_sim.* (+ labels, hm3_snps.txt, w_hm3_sim.snplist, 1kg_sim_rs.pvar),
  phenotypes/{height,bmi}.tsv, truth.json

Planted for the checks to catch: an unsampled European-like population, an
African-European continuum, nuclear families, an unlisted long-range LD block
(chr6:70-90 Mb), A/T and C/G SNPs, high-missingness and HWE-failing SNPs, unit errors.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess

import numpy as np
import pandas as pd

SUPER = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN"],
    "EUR": ["CEU", "GBR", "FIN", "TSI", "IBS"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
}
F_SUPER, F_SUB = 0.10, 0.006
F_SUB_SPECIAL = {"FIN": 0.015}
N_REF_PER_POP = 100
N_REF_KIDS = 10                  # trio children in the reference (not founders)

N_SNPS = 24_000
INV_CHR, INV_LO, INV_HI, N_INV = 6, 70_000_000, 90_000_000, 300
FRAC_HM3 = 0.7

# cohort founders by origin
COHORT = {"CEU": 1100, "GBR": 1100, "IBS": 300, "TSI": 300, "FIN": 150,
          "ASH": 300,                       # unsampled European-like population
          "ADMIX": 1100,                    # African-European continuum
          "YRI": 150, "ESN": 100,
          "CHB": 200, "JPT": 150, "CHS": 100,
          "GIH": 150, "PJL": 100}
N_COUPLES = 260                   # nuclear families within CEU/GBR, 2 children each
N_SIBSHIPS = 100                  # sib pairs whose parents are not in the cohort


def bn(p, F, rng):
    """Balding-Nichols draw around ancestral frequencies p."""
    a = p * (1 - F) / F
    b = (1 - p) * (1 - F) / F
    return np.clip(rng.beta(a, b), 0.002, 0.998)


def genotypes(freq, n, rng):
    return rng.binomial(2, freq, size=(n, len(freq))).astype(np.int8)


def write_bed(prefix, G, fam, bim):
    """PLINK 1 bed, variant-major. G counts the A1 allele (= ALT, column 5)."""
    n, m = G.shape
    code = np.full(G.shape, 1, dtype=np.uint8)            # 01 missing
    code[G == 2] = 0                                      # 00 hom A1
    code[G == 1] = 2                                      # 10 het
    code[G == 0] = 3                                      # 11 hom A2
    pad = (-n) % 4
    if pad:
        code = np.vstack([code, np.zeros((pad, m), dtype=np.uint8)])
    c = code.T.reshape(m, -1, 4)
    packed = (c[..., 0] | (c[..., 1] << 2) | (c[..., 2] << 4) | (c[..., 3] << 6)).astype(np.uint8)
    with open(prefix + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.tobytes())
    fam.to_csv(prefix + ".fam", sep="\t", header=False, index=False)
    bim.to_csv(prefix + ".bim", sep="\t", header=False, index=False)


def write_acaf_pgen(prefix, G, iids, chrom, pos, ref, alt, alt2, carriers, plink2):
    """ACAF as AoU ships it: pgen, IDs chr:pos:REF, and a second ALT joined onto the row (REF -> ALT,ALT2)."""
    gt = np.array(["0/0", "0/1", "1/1"])
    S = np.where((G >= 0) & (G <= 2), gt[np.clip(G, 0, 2)], "./.").astype("<U3")
    S[carriers] = "0/2"
    alts = np.where(alt2 != "", alt + "," + alt2, alt)
    with open(prefix + ".vcf", "w") as fh:
        fh.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(iids) + "\n")
        for j in range(len(pos)):
            fh.write(f"chr{chrom}\t{pos[j]}\tchr{chrom}:{pos[j]}:{ref[j]}\t{ref[j]}\t{alts[j]}\t.\tPASS\t.\tGT\t"
                     + "\t".join(S[:, j]) + "\n")
    subprocess.run([plink2, "--vcf", prefix + ".vcf", "--output-chr", "chrM", "--make-pgen",
                    "--out", prefix, "--silent"], check=True)
    os.remove(prefix + ".vcf")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=".sim/inputs")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--plink2", default=os.environ.get("PLINK2") or (".tools/bin/plink2" if os.path.exists(".tools/bin/plink2")
                                                                      else shutil.which("plink2")))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = a.out
    for d in ("acaf", "aux/ancestry", "aux/relatedness", "reference", "phenotypes"):
        os.makedirs(f"{out}/{d}", exist_ok=True)

    # ---- variants ---------------------------------------------------------------
    chr_len = np.array([248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133,
                        114, 107, 102, 90, 83, 80, 59, 64, 47, 51]) * 1_000_000
    n_per = np.round(chr_len / chr_len.sum() * (N_SNPS - N_INV)).astype(int)
    n_per[0] += (N_SNPS - N_INV) - n_per.sum()
    chrom, pos = [], []
    for c, (L, k) in enumerate(zip(chr_len, n_per), start=1):
        p = np.sort(rng.choice(np.arange(1_000_000, L - 1_000_000, 37), k, replace=False))
        if c == INV_CHR:                      # keep the inversion region for the block
            p = p[(p < INV_LO) | (p > INV_HI)]
            extra = k - len(p)
            p = np.sort(np.r_[p, rng.choice(np.arange(INV_HI + 10, L - 1_000_000, 37),
                                            extra, replace=False)])
        chrom += [c] * k
        pos += list(p)
    inv_pos = np.sort(rng.choice(np.arange(INV_LO, INV_HI, 37), N_INV, replace=False))
    chrom += [INV_CHR] * N_INV
    pos += list(inv_pos)
    var = pd.DataFrame({"chr": chrom, "pos": pos}).sort_values(["chr", "pos"]).reset_index(drop=True)
    m = len(var)
    is_inv = ((var["chr"] == INV_CHR) & var["pos"].between(INV_LO, INV_HI)).to_numpy()
    pairs = [("A", "C"), ("A", "G"), ("A", "T"), ("C", "A"), ("C", "G"), ("C", "T"),
             ("G", "A"), ("G", "C"), ("G", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
    al = np.array(pairs)[rng.integers(0, 12, m)]
    var["ref"], var["alt"] = al[:, 0], al[:, 1]
    var["id"] = (var["chr"].astype(str) + ":" + var["pos"].astype(str)
                 + ":" + var["ref"] + ":" + var["alt"])
    bim = pd.DataFrame({"chr": var["chr"], "id": var["id"], "cm": 0, "pos": var["pos"],
                        "a1": var["alt"], "a2": var["ref"]})

    # ---- frequencies ---------------------------------------------------------------
    p_anc = rng.uniform(0.05, 0.95, m)
    f_super = {s: bn(p_anc, F_SUPER, rng) for s in SUPER}
    f_pop = {}
    for s, pops in SUPER.items():
        for pop in pops:
            f_pop[pop] = bn(f_super[s], F_SUB_SPECIAL.get(pop, F_SUB), rng)
    nwe = bn(f_super["EUR"], 0.004, rng)          # CEU and GBR share a NW-European parent
    f_pop["CEU"], f_pop["GBR"] = bn(nwe, 0.0005, rng), bn(nwe, 0.0005, rng)
    f_pop["ASH"] = bn(f_super["EUR"], 0.02, rng)
    # the inversion: each haplotype is standard or inverted; inverted haplotypes
    # carry shifted allele frequencies at every block SNP
    shift = rng.choice([-1, 1], N_INV) * rng.uniform(0.25, 0.4, N_INV)
    inv_freq = {"AFR": 0.10, "EUR": 0.30, "EAS": 0.05, "SAS": 0.20}
    pop_super = {pop: s for s, pops in SUPER.items() for pop in pops}
    pop_super["ASH"] = "EUR"

    def draw(pop, n):
        g = genotypes(f_pop[pop], n, rng)
        q = inv_freq[pop_super[pop]]
        base = f_pop[pop][is_inv]
        alt = np.clip(base + shift, 0.01, 0.99)
        for h in range(2):
            inv = rng.random(n) < q
            fr = np.where(inv[:, None], alt, base)
            hap = (rng.random((n, N_INV)) < fr).astype(np.int8)
            if h == 0:
                g[:, is_inv] = hap
            else:
                g[:, is_inv] += hap
        return g

    # ---- reference panel -------------------------------------------------------------
    ref_G, ref_rows = [], []
    sid = 0
    for s, pops in SUPER.items():
        for pop in pops:
            g = draw(pop, N_REF_PER_POP)
            ref_G.append(g)
            for _ in range(N_REF_PER_POP):
                sid += 1
                ref_rows.append({"SampleID": f"SIM{sid:05d}", "FatherID": "0", "MotherID": "0",
                                 "Sex": int(rng.integers(1, 3)), "Population": pop,
                                 "Superpopulation": s})
    ref_G = np.vstack(ref_G)
    ref = pd.DataFrame(ref_rows)
    kids = []
    ceu = np.flatnonzero(ref["Population"] == "CEU")
    for k in range(N_REF_KIDS):
        fa, mo = ceu[2 * k], ceu[2 * k + 1]
        kids.append(rng.binomial(1, ref_G[fa] / 2) + rng.binomial(1, ref_G[mo] / 2))
        sid += 1
        ref_rows.append({"SampleID": f"SIM{sid:05d}", "FatherID": ref.at[fa, "SampleID"],
                         "MotherID": ref.at[mo, "SampleID"], "Sex": int(rng.integers(1, 3)),
                         "Population": "CEU", "Superpopulation": "EUR"})
    ref_G = np.vstack([ref_G, np.array(kids, dtype=np.int8)])
    ref = pd.DataFrame(ref_rows)
    ref[["SampleID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"]].to_csv(
        f"{out}/reference/1kg_sim_ped_population.txt", sep=" ", index=False)
    write_bed(f"{out}/reference/1kg_sim", ref_G,
              pd.DataFrame({"fid": ref["SampleID"], "iid": ref["SampleID"], "pat": 0, "mat": 0,
                            "sex": ref["Sex"], "phe": -9}), bim)

    # ---- cohort founders ----------------------------------------------------------------
    G, origin, afr_frac = [], [], []
    for pop, n in COHORT.items():
        if pop == "ADMIX":
            alpha = rng.beta(6, 1.8, n)                    # African ancestry fraction
            fa = (f_pop["YRI"] + f_pop["ESN"]) / 2
            fe = (f_pop["CEU"] + f_pop["GBR"]) / 2
            g = np.zeros((n, m), dtype=np.int8)
            for _ in range(2):
                from_afr = rng.random((n, m)) < alpha[:, None]
                g += (rng.random((n, m)) < np.where(from_afr, fa, fe)).astype(np.int8)
            G.append(g)
            afr_frac += list(alpha)
        else:
            G.append(draw(pop, n))
            afr_frac += [1.0 if pop_super[pop] == "AFR" else 0.0] * n
        origin += [pop] * n
    G = np.vstack(G)
    origin = np.array(origin, dtype=object)
    afr_frac = np.array(afr_frac)
    household = np.arange(len(G))

    # ---- families --------------------------------------------------------------------------
    rel = []
    nw = np.flatnonzero(np.isin(origin, ["CEU", "GBR"]))
    parents = rng.choice(nw, 2 * N_COUPLES, replace=False)
    kids_G, kid_origin, kid_afr, kid_house = [], [], [], []
    n0 = len(G)

    def child(g1, g2):
        return (rng.binomial(1, g1 / 2) + rng.binomial(1, g2 / 2)).astype(np.int8)

    for c in range(N_COUPLES):
        fa, mo = parents[2 * c], parents[2 * c + 1]
        idx = []
        for _ in range(2):
            kids_G.append(child(G[fa], G[mo]))
            kid_origin.append(origin[fa])
            kid_afr.append(0.0)
            kid_house.append(household[fa])
            idx.append(n0 + len(kids_G) - 1)
        household[mo] = household[fa]
        rel += [(fa, idx[0], "PO"), (fa, idx[1], "PO"), (mo, idx[0], "PO"),
                (mo, idx[1], "PO"), (idx[0], idx[1], "FS")]
    for c in range(N_SIBSHIPS):
        pop = rng.choice(["CEU", "GBR", "IBS", "TSI"])
        g1, g2 = draw(pop, 2)
        idx = []
        for _ in range(2):
            kids_G.append(child(g1, g2))
            kid_origin.append(pop)
            kid_afr.append(0.0)
            kid_house.append(10_000_000 + c)
            idx.append(n0 + len(kids_G) - 1)
        rel.append((idx[0], idx[1], "FS"))
    G = np.vstack([G, np.array(kids_G, dtype=np.int8)])
    origin = np.r_[origin, kid_origin]
    afr_frac = np.r_[afr_frac, kid_afr]
    household = np.r_[household, kid_house]
    n = len(G)
    perm = rng.permutation(n)                     # don't leave families adjacent
    inv_perm = np.argsort(perm)
    G, origin, afr_frac, household = G[perm], origin[perm], afr_frac[perm], household[perm]
    rel = [(inv_perm[i], inv_perm[j], t) for i, j, t in rel]
    pid = (1_000_000 + rng.choice(9_000_000, n, replace=False)).astype(str)

    # ---- genotype artefacts ------------------------------------------------------------------
    G = G.copy()
    G[rng.random(G.shape) < 0.002] = -1
    bad_miss = rng.choice(m, int(0.01 * m), replace=False)
    G[:, bad_miss] = np.where(rng.random((n, len(bad_miss))) < 0.08, -1, G[:, bad_miss])
    bad_hwe = rng.choice(np.setdiff1d(np.arange(m), bad_miss), int(0.005 * m), replace=False)
    sub = G[:, bad_hwe]
    sub[(sub == 2) & (rng.random(sub.shape) < 0.5)] = 1
    sub[(sub == 0) & (rng.random(sub.shape) < 0.15)] = 1
    G[:, bad_hwe] = sub

    # ---- HapMap3 ---------------------------------------------------------------------------------
    hm3 = rng.random(m) < FRAC_HM3
    hm3 |= is_inv                                  # the inversion SNPs are "HM3" too
    var.loc[hm3, "id"].to_csv(f"{out}/reference/hm3_snps.txt", index=False, header=False)
    # rsID versions for the hm3_overlap diagnostic; separate RNG so nothing else changes
    r2 = np.random.default_rng(a.seed + 1000)
    rs = "rs" + (1_000_000 + np.arange(m)).astype(str)
    pv = pd.DataFrame({"#CHROM": var["chr"], "POS": var["pos"], "ID": rs, "REF": var["ref"], "ALT": var["alt"]})
    swap = hm3 & (r2.random(m) < 0.005)             # allele order disagrees with ACAF
    pv.loc[swap, ["REF", "ALT"]] = pv.loc[swap, ["ALT", "REF"]].to_numpy()
    multi = hm3 & ~swap & (r2.random(m) < 0.003)     # multi-allelic in the reference
    pv.loc[multi, "ALT"] = pv.loc[multi, "ALT"] + ",N"
    pv.to_csv(f"{out}/reference/1kg_sim_rs.pvar", sep="\t", index=False)
    unmapped = "rs" + (90_000_000 + np.arange(int(0.02 * hm3.sum()))).astype(str)
    pd.DataFrame({"SNP": np.r_[rs[hm3], unmapped]}).to_csv(
        f"{out}/reference/w_hm3_sim.snplist", sep="\t", index=False)

    sex = rng.integers(1, 3, n)
    fam = pd.DataFrame({"fid": 0, "iid": pid, "pat": 0, "mat": 0, "sex": sex, "phe": -9})
    # a rare second ALT at 30% of HM3 sites, on the same row: only splitting recovers those SNPs
    r3 = np.random.default_rng(a.seed + 2000)
    joined = hm3 & (r3.random(m) < 0.30)
    other = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
    alt2 = np.array([next(x for x in other[r] if x != t) if j else "" for r, t, j in zip(var["ref"], var["alt"], joined)])
    carriers = (G == 0) & joined[None, :] & (r3.random(G.shape) < 0.005)
    for c in range(1, 23):
        k = (var["chr"] == c).to_numpy()
        write_acaf_pgen(f"{out}/acaf/acaf_threshold.chr{c}", G[:, k], pid, c, var["pos"].to_numpy()[k],
                        var["ref"].to_numpy()[k], var["alt"].to_numpy()[k], alt2[k], carriers[:, k], a.plink2)

    # ---- premade ancestry PCs: fit on the reference, project the cohort ---------------------------------
    founders = (ref["FatherID"] == "0").to_numpy()
    ok = ~is_inv
    p_ref = ref_G[founders][:, ok].mean(0) / 2
    keep = (p_ref > 0.05) & (p_ref < 0.95)
    p_ref = p_ref[keep]
    sd = np.sqrt(2 * p_ref * (1 - p_ref))
    X = (ref_G[:, ok][:, keep] - 2 * p_ref) / sd
    U, S, Vt = np.linalg.svd(X[founders], full_matrices=False)
    V = Vt[:16].T
    ref_scores = X @ V
    Gc = G[:, ok][:, keep].astype(float)
    Gc[Gc < 0] = np.nan
    Xc = np.nan_to_num((Gc - 2 * p_ref) / sd)
    cohort_scores = Xc @ V
    centroid = {s: ref_scores[founders & (ref["Superpopulation"] == s).to_numpy(), :4].mean(0)
                for s in SUPER}
    pred = [min(centroid, key=lambda s: np.sum((x[:4] - centroid[s]) ** 2)).lower()
            for x in cohort_scores]
    fmt = lambda row: "[" + ", ".join(f"{v:.6f}" for v in row) + "]"
    pd.DataFrame({"research_id": pid, "ancestry_pred": pred,
                  "pca_features": [fmt(r) for r in cohort_scores]}).to_csv(
        f"{out}/aux/ancestry/ancestry_preds.tsv", sep="\t", index=False)
    pd.DataFrame({"s": ref["SampleID"], "scores": [fmt(r) for r in ref_scores]}).to_csv(
        f"{out}/aux/ancestry/training_pca.tsv", sep="\t", index=False)

    # ---- relatedness (AoU reports kinship > 0.1 only) --------------------------------------------------------
    kin = {"PO": 0.25, "FS": 0.25}
    pd.DataFrame({"i.s": [pid[i] for i, _, _ in rel], "j.s": [pid[j] for _, j, _ in rel],
                  "kin": [kin[t] + rng.normal(0, 0.01) for *_, t in rel]}).to_csv(
        f"{out}/aux/relatedness/samples_relatedness.tsv", sep="\t", index=False)

    # ---- phenotypes -------------------------------------------------------------------------------------------
    # causal SNPs on HapMap3 only, so a HapMap3 GRM can in principle recover h2
    causal = rng.choice(np.flatnonzero(hm3 & ~is_inv), 2000, replace=False)
    Gq = G[:, causal].astype(float)
    Gq[Gq < 0] = np.nan
    Gq = np.where(np.isnan(Gq), np.nanmean(Gq, 0), Gq)
    age = rng.uniform(18, 85, n).round(1)
    house_ids, house_idx = np.unique(household, return_inverse=True)

    def trait(h2, c2):
        beta = rng.normal(0, 1, len(causal))
        g = (Gq - Gq.mean(0)) @ beta
        nw_ = np.isin(origin, ["CEU", "GBR"])          # h2 is defined within CEU+GBR
        g = (g - g[nw_].mean()) / g[nw_].std()
        c = rng.normal(0, 1, len(house_ids))[house_idx]
        e = rng.normal(0, 1, n)
        return np.sqrt(h2) * g + np.sqrt(c2) * c + np.sqrt(1 - h2 - c2) * e

    truth = {"height": {"h2": 0.65, "c2": 0.05}, "bmi": {"h2": 0.30, "c2": 0.10}}
    height = 163 + 13 * (sex == 1) - 0.04 * (age - 50) - 2.0 * afr_frac + 6.5 * trait(truth["height"]["h2"], truth["height"]["c2"])
    bmi = 28 + 0.03 * (age - 50) + 1.5 * afr_frac + 5 * trait(truth["bmi"]["h2"], truth["bmi"]["c2"])
    height[rng.random(n) < 0.003] *= 10                 # unit errors (mm entered as cm)
    bmi[rng.random(n) < 0.002] = 250
    sexlab = np.where(sex == 1, "Male", "Female")
    for name, y in (("height", height), ("bmi", bmi)):
        has = rng.random(n) > 0.08
        pd.DataFrame({"person_id": pid[has], "phenotype": y[has].round(2),
                      "age": age[has], "sex_at_birth": sexlab[has]}).to_csv(
            f"{out}/phenotypes/{name}.tsv", sep="\t", index=False)

    import json
    with open(f"{out}/truth.json", "w") as fh:
        json.dump({**truth, "note": "h2 and c2 are variance fractions within CEU+GBR, after sex and age; "
                   "c2 is a household effect shared by parents and children"}, fh, indent=2)
    print(f"simulated {n:,} participants, {len(ref):,} reference samples, {m:,} variants -> {out}")
    print(f"  {len(rel)} related pairs, inversion {N_INV} SNPs on chr{INV_CHR}:{INV_LO//10**6}-{INV_HI//10**6} Mb")


if __name__ == "__main__":
    main()
