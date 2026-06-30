import pandas as pd
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# ── Load data ─────────────────────────────────────────────────────────

# PCA scores — has #FID and IID columns
scores = pd.read_csv('eur_pca.eigenvec', sep='\t')
# Rename IID to sample for the merge
scores = scores.rename(columns={'IID': 'sample'})
scores = scores.drop(columns=['#FID'])

# Population labels — has sample, pop, super_pop, gender
labels = pd.read_csv(
    'integrated_call_samples_v3.20130502.ALL.panel',
    sep='\t'
)

scores = scores.merge(labels[['sample', 'pop', 'super_pop']], on='sample')

print("Samples per population:")
print(scores['pop'].value_counts())

# ── Fit ellipsoid to CEU+GBR ──────────────────────────────────────────

N_PCS    = 2    # can increase this to 3 or more, but visualization will be limited to 2D
pc_cols  = [f'PC{i}' for i in range(1, N_PCS + 1)]

ref_mask = scores['pop'].isin(['CEU', 'GBR'])
ref_pcs  = scores.loc[ref_mask, pc_cols].values

ref_mean    = ref_pcs.mean(axis=0)
ref_cov     = np.cov(ref_pcs, rowvar=False)
ref_cov_inv = np.linalg.inv(ref_cov)

def mahal(x, mean, cov_inv):
    d = x - mean
    return np.sqrt(d @ cov_inv @ d)

all_pcs       = scores[pc_cols].values
scores['mahal'] = [mahal(row, ref_mean, ref_cov_inv) for row in all_pcs]

# Chi-squared threshold for 90% confidence ellipsoid
threshold        = np.sqrt(chi2.ppf(0.9, df=N_PCS))
scores['pass']   = scores['mahal'] <= threshold

print(f"\nMahalanobis threshold ({N_PCS} PCs, 90%): {threshold:.3f}")

# ── Validate ──────────────────────────────────────────────────────────

print("\nPass rate by population:")
summary = scores.groupby('pop')['pass'].agg(['sum', 'count'])
summary['pct'] = (summary['sum'] / summary['count'] * 100).round(1)
print(summary.sort_values('pct', ascending=False))

# Key checks
ceu_gbr_rate = scores.loc[ref_mask, 'pass'].mean()
fin_rate     = scores.loc[scores['pop'] == 'FIN', 'pass'].mean()
tsi_rate     = scores.loc[scores['pop'] == 'TSI', 'pass'].mean()
ibs_rate     = scores.loc[scores['pop'] == 'IBS', 'pass'].mean()

print(f"\nCEU+GBR retention : {ceu_gbr_rate:.1%}  (expect ~90%)")
print(f"FIN exclusion      : {1-fin_rate:.1%}  (want high)")
print(f"TSI exclusion      : {1-tsi_rate:.1%}  (want high)")
print(f"IBS exclusion      : {1-ibs_rate:.1%}  (want high)")

# ── Plot ──────────────────────────────────────────────────────────────

colors = {
    'CEU': '#2166ac',
    'GBR': '#4dac26',
    'FIN': '#d01c8b',
    'TSI': '#f1a340',
    'IBS': '#998ec3',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for pop, grp in scores.groupby('pop'):
    for ax in axes:
        ax.scatter([], [])  # placeholder for consistent colors

# Left plot: PC1 vs PC2 colored by population
for pop, grp in scores.groupby('pop'):
    axes[0].scatter(
        grp['PC1'], grp['PC2'],
        c=colors.get(pop, 'grey'),
        label=pop, alpha=0.7, s=20
    )

# Draw 90% ellipse on PC1/PC2 plane
cov_2d   = np.cov(ref_pcs[:, :2], rowvar=False)
mean_2d  = ref_pcs[:, :2].mean(axis=0)
eigvals, eigvecs = np.linalg.eigh(cov_2d)
angle    = np.degrees(np.arctan2(*eigvecs[:, 1][::-1]))
scale_2d = np.sqrt(chi2.ppf(0.9, df=2))
width    = 2 * scale_2d * np.sqrt(eigvals[1])
height   = 2 * scale_2d * np.sqrt(eigvals[0])
ellipse  = Ellipse(
    xy=mean_2d, width=width, height=height, angle=angle,
    fill=False, edgecolor='black', linewidth=2,
    linestyle='--', label='90% CEU+GBR ellipsoid'
)
axes[0].add_patch(ellipse)
axes[0].set_xlabel('PC1'); axes[0].set_ylabel('PC2')
axes[0].set_title('PC1 vs PC2 — 1000G EUR')
axes[0].legend(markerscale=2, fontsize=8)

# Right plot: Mahalanobis distance distribution by population
for pop, grp in scores.groupby('pop'):
    axes[1].hist(
        grp['mahal'], bins=30, alpha=0.5,
        label=pop, color=colors.get(pop, 'grey')
    )
axes[1].axvline(
    threshold, color='black', linestyle='--', linewidth=2,
    label=f'Threshold ({threshold:.2f})'
)
axes[1].set_xlabel('Mahalanobis distance from CEU+GBR centroid')
axes[1].set_ylabel('Count')
axes[1].set_title(f'Ellipsoid distance distribution ({N_PCS} PCs)')
axes[1].legend(fontsize=8)

plt.suptitle(
    f'Round 2 filtering validation — 1000G EUR chr22\n'
    f'CEU+GBR retention: {ceu_gbr_rate:.1%}  |  '
    f'FIN exclusion: {1-fin_rate:.1%}  |  '
    f'TSI exclusion: {1-tsi_rate:.1%}',
    fontsize=11
)
plt.tight_layout()
plt.savefig('round2_validation.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nPlot saved to round2_validation.png")