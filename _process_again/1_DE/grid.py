import os, glob
import scanpy as sc
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import multiprocessing as mp
   
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
from pydeseq2.dds            import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds             import DeseqStats
pb_df   = counts_df
pb_meta = meta.loc[counts_df.index]

# define your four comparisons: (subset_str, design_formula, contrast, out_filename)
jobs = [
    # Female AD vs Female Control
    ("Female",   "~ condition", ["condition", "Disease", "Control"], "F_AD_vs_F_Control.csv"),
    # Male   AD vs Male   Control
    ("Male",     "~ condition", ["condition", "Disease", "Control"], "M_AD_vs_M_Control.csv"),
    # Female vs Male within Disease
    ("Disease",  "~ gender",    ["gender",    "Female",  "Male"   ], "F_vs_M_Disease.csv"),
    # Female vs Male within Control
    ("Control",  "~ gender",    ["gender",    "Female",  "Male"   ], "F_vs_M_Control.csv"),
]

for subset_var, design, contrast, outfn in jobs:
    # subset samples
    if design.startswith("~ condition"):
        # gender‐stratified: only Female or only Male
        mask = (pb_meta["gender"] == subset_var)
    else:
        # condition‐stratified: only Disease or only Control
        mask = (pb_meta["condition"] == subset_var)

    counts_sub = pb_df.loc[mask]
    meta_sub   = pb_meta.loc[mask]

    # build and fit the DESeq2 dataset
    dds = DeseqDataSet(
        counts      = counts_sub,
        metadata    = meta_sub,
        design      = design,
        refit_cooks = True
    )
    dds.fit_size_factors()
    dds.deseq2()
    dds.fit_LFC()

    # run the contrast
    inf   = DefaultInference(n_cpus=1)
    stats = DeseqStats(dds, contrast=contrast, inference=inf)
    stats.summary()

    # write out a CSV of all genes + DE stats
    stats.results_df.to_csv(outfn)
    print(f"Wrote {outfn}")