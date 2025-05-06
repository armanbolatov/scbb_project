import os, glob
import scanpy as sc
import pandas as pd
import memento

# 1) Load sample metadata
meta = (
    pd.read_csv('/l/users/darya.taratynova/scbb_project/_process_again/dataset_meta.csv')
      .set_index('sample_id')
)

# 2) Read & concatenate all your AnnData, annotating with gender & a numeric covariate
paths = glob.glob('/l/users/darya.taratynova/scbb_project/_process_again/**/ann_data_mapped/*.h5ad', recursive=True)
adatas = []
for p in paths:
    sample = os.path.basename(p).replace('.h5ad','')
    ad = sc.read_h5ad(p)
    ad.obs['gender']       = meta.loc[sample, 'gender']
    # example covariate: total UMI count per cell
    ad.obs['umi_counts']   = ad.X.sum(axis=1).A1 if hasattr(ad.X, 'A1') else ad.X.sum(axis=1)
    adatas.append(ad)
print(f'Loaded {len(adatas)} samples')
adata_all = sc.concat(adatas, join='inner', index_unique=None)

# 1) Fix duplicate obs names
adata_all.obs_names_make_unique()
print(f'Concatenated {len(adata_all.obs_names)} cells')
# 2) Compute a capture_rate in [0,1]
tc = adata_all.X.sum(axis=1)
rates = tc.A1 if hasattr(tc,'A1') else tc

# get the max once
max_rate = rates.max()

# scale so max is 0.999 (or any number <1)
adata_all.obs['capture_rate'] = rates / (max_rate * 1.001)
print(f'Capture rate: {adata_all.obs["capture_rate"].min():.2f} - {adata_all.obs["capture_rate"].max():.2f}')
# 3) Setup memento

sc.pp.normalize_total(adata_all)
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, n_top_genes=2000, flavor='seurat_v3')
adata_small = adata_all[:, adata_all.var['highly_variable']].copy()

from memento.wrappers import binary_test_1d

from memento.wrappers import binary_test_1d

res_gender = binary_test_1d(
    adata_small,                         # your AnnData
    adata_small.obs['capture_rate'],     # per-cell q
    'gender',                            # column to compare
    1,                                   # num_cpus (positional)
    num_boot=200,                        # override default bootstraps
    verbose=1
)

res_gender.to_csv('DE_gender_memento.csv', index=False)
print(res_gender.sort_values('de_pval').head(10))


res_gender.to_csv('DE_gender_memento.csv', index=False)
print(res_gender.sort_values('de_pval').head(10))
