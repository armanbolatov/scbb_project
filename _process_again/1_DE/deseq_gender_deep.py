import os, glob
import scanpy as sc
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import multiprocessing as mp

mp.set_start_method('spawn', force=True)

def main():
    meta = (
        pd.read_csv('/l/users/darya.taratynova/scbb_project/_process_again/dataset_meta.csv')
        .set_index('sample_id')
    )

    input_root = '/l/users/darya.taratynova/scbb_project/_process_again'
    h5ad_files = glob.glob(os.path.join(input_root, '**', 'ann_data_mapped', '*.h5ad'), recursive=True)

    counts_list = []
    sample_ids   = []

    for f in h5ad_files:
        sample = os.path.basename(f).replace('.h5ad','')
        adata  = sc.read_h5ad(f)

        if hasattr(adata.X, 'sum'): 
            gene_counts = adata.X.sum(axis=0)
            if hasattr(gene_counts, 'A1'):
                gene_counts = gene_counts.A1
        else:
            gene_counts = adata.X.sum(axis=0)
        counts_list.append(pd.Series(gene_counts, index=adata.var_names))
        sample_ids.append(sample)
        print(f"Loaded {f} with {adata.shape[0]} cells and {adata.shape[1]} genes")
    print(f"Loaded {len(h5ad_files)} files with {len(counts_list)} samples")

    counts_df = pd.DataFrame(counts_list, index=sample_ids)
    counts_df = counts_df.dropna(axis=1)
    counts_df.index.name = 'sample_id'

    col_data = meta.loc[counts_df.index]
    disease = col_data['condition']=='Control'
    counts_df_d   = counts_df.loc[disease]
    col_data_d    = col_data.loc[disease]
    mask = col_data['condition']=='Control'
    counts_df_d = counts_df.loc[mask]
    col_data_d  = col_data.loc[mask]

    # Build a simple gender‐only model
    dds = DeseqDataSet(
        counts=counts_df_d,
        metadata=col_data_d,
        design="~ gender",
        refit_cooks=True
    )
    dds.fit_size_factors()
    dds.deseq2()
    dds.fit_LFC()

    inf   = DefaultInference(n_cpus=1)
    # Male vs Female
    stats = DeseqStats(
        dds,
        contrast=["gender", "Male", "Female"],
        inference=inf
    )
    stats.summary()
    res = stats.results_df
    res.to_csv('DE_gender_within_control.csv')
if __name__ == "__main__":
    main()