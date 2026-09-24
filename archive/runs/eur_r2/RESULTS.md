# eur_r2 — results

What each step actually produced. `provenance.result_line()` prints the row;
paste it here.

Parameters and aggregate counts only — nothing that has not cleared AoU egress
review. If a number can't be written down, record its shape (`n=222,4xx`) and
leave the exact value in the bucket provenance file.

| date | step | commit | numbers |
|---|---|---|---|
| 2026-09-23 | gate | 2d5c9c0 | round1 n=275,717 → round2 n=223,209 (52,508 dropped by CEU+GBR gate); K_PCS=5; round2 is a strict subset |
| 2026-09-24 | covariate_pcs | 2d5c9c0 | n=223,209; 64,379 HM3-common variants (r²=0.1); K_PCS=5; PC1 16.90%, PC2 13.95%, PC6–20 noise |

## Notes

Anything that needed a judgement call, a rerun, or an explanation — one short
entry per event, newest last.

- 2026-09-24: covariate PCA uses r²=0.1 (not 0.05) — at 0.05 HM3 common loses 77% of variants
  (97K→23K) because HM3 was designed to tag LD blocks. At 0.1 we get 64,379 variants and stable
  PC estimates. Loadings checked — no LD artefact peaks in PCs 1–4.
- 2026-09-23: gate drops 52,508 participants (19%) from round 1. Round 2 is a strict subset — no participant appears in round 2 that wasn't in round 1.
