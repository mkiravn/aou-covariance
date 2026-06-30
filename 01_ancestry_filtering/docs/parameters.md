# Analysis Parameters

## Round 2 ancestry filtering

Validated on 1000G EUR chr22 (2026-06-30).

| Setting | Value |
|---|---|
| N PCs | 6 |
| Threshold | chi2(0.9, df=6) = Mahalanobis 3.263 |
| Reference cluster | CEU + GBR |

### Validation results (1000G EUR chr22)
| Population | N | Pass | % |
|---|---|---|---|
| CEU | 99 | 91 | 91.9% |
| GBR | 91 | 77 | 84.6% |
| FIN | 99 | 5 | 5.1% |
| TSI | 107 | 32 | 29.9% |
| IBS | 107 | 57 | 53.3% |

CEU+GBR retention: 88.4%
FIN exclusion: 94.9%

Note: chr22-only validation — genome-wide expected to be tighter.

