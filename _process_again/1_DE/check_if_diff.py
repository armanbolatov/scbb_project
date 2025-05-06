import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1) Load both result sets
base = pd.read_csv(
    '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/results/condition_basic.csv',
    index_col=0
).rename(columns={
    'log2FoldChange':'LFC_noGender',
    'padj':'padj_noGender'
})

with_gender = pd.read_csv(
    '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/results/gender.csv',
    index_col=0
).rename(columns={
    'log2FoldChange':'LFC_withGender',
    'padj':'padj_withGender'
})
base = base.rename(columns={
    'baseMean':'baseMean_noG',
    'lfcSE':  'lfcSE_noG',
    'stat':   'stat_noG',
    'pvalue': 'pvalue_noG'
    # LFC and padj already renamed earlier
})

with_gender = with_gender.rename(columns={
    'baseMean':'baseMean_wG',
    'lfcSE':  'lfcSE_wG',
    'stat':   'stat_wG',
    'pvalue': 'pvalue_wG'
})


# 2) Merge on gene name
cmp = base.join(with_gender, how="inner")

# 3) Compute differences and correlations
cmp['ΔLFC']   = cmp['LFC_withGender'] - cmp['LFC_noGender']
corr = cmp[['LFC_noGender','LFC_withGender']].corr().iloc[0,1]
print(f"Pearson r between LFCs: {corr:.4f}")
print(f"Mean absolute ΔLFC: {cmp['ΔLFC'].abs().mean():.4f}")

# 4) How many genes changed significance?
sig1 = set(base.query('padj_noGender < 0.05').index)
sig2 = set(with_gender.query('padj_withGender < 0.05').index)
print(f"sig w/o gender: {len(sig1)}, sig w/ gender: {len(sig2)}")
print(f"Genes only uncovered by gender model: {len(sig2 - sig1)}")
print(f"Genes only lost by gender model:      {len(sig1 - sig2)}")

# 5) Quick scatter of LFCs
plt.figure(figsize=(5,5))
plt.scatter(cmp['LFC_noGender'], cmp['LFC_withGender'], alpha=0.4, s=5)
lims = np.percentile(np.concatenate([cmp['LFC_noGender'],cmp['LFC_withGender']]), [1,99])
plt.plot(lims, lims, 'k--', linewidth=1)
plt.xlabel('LFC (no gender)')
plt.ylabel('LFC (+ gender)')
plt.title('Comparison of log₂FC')
plt.tight_layout()
plt.show()
