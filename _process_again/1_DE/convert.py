import scanpy as sc
import pandas as pd
import glob, os

# 1) Load your DE results (indexed by gene name)
de = pd.read_csv(
    '/l/users/darya.taratynova/scbb_project/_process_again/1_DE/results/condition_basic.csv',
    index_col=0
)

# 2) List of folders containing your mapped AnnDatas
base_dirs = [
    '/l/users/darya.taratynova/scbb_project/_process_again/GSE147528/ann_data_mapped',
    '/l/users/darya.taratynova/scbb_project/_process_again/GSE160936/ann_data_mapped',
    '/l/users/darya.taratynova/scbb_project/_process_again/GSE174367/ann_data_mapped',
]

for base in base_dirs:
    for path in glob.glob(os.path.join(base, '*.h5ad')):
        # load
        adata = sc.read_h5ad(path)

        # compute intersection of genes
        common = adata.var_names.intersection(de.index)

        # annotate var with your DE metrics
        adata.var['log2FoldChange'] = pd.NA
        adata.var['padj']           = pd.NA
        adata.var.loc[common, 'log2FoldChange'] = de.loc[common, 'log2FoldChange']
        adata.var.loc[common, 'padj']           = de.loc[common, 'padj']

        # optionally add a per-cell “DE signature score”:
        # from scanpy import tl
        # sigs = de.query('padj<0.05 and abs(log2FoldChange)>1').index.tolist()
        # tl.score_genes(adata, gene_list=sigs, score_name='DE_score', use_raw=False)

        # write out
        out = path.replace('.h5ad','_de_annotated.h5ad')
        adata.write_h5ad(out)
        print("Wrote:", out)
