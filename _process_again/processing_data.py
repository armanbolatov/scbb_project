import os
import gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scvi
from scvi.external import MRVI
from scipy import sparse
import anndata as ad
from tqdm import tqdm
from anndata import AnnData

# 1) Read metadata
meta = pd.read_csv("dataset_meta_needed.csv", dtype=str)

# 2) Load each .h5ad, add obs columns, collect into list
adatas = []
for _, row in meta.iterrows():
    dataset = row["dataset"]
    sample = row["sample_id"]
    in_path = os.path.join(dataset, "ann_data", f"{sample}.h5ad")
    if not os.path.exists(in_path):
        print(f"⚠️ File not found: {in_path}, skipping")
        continue

    print(f"Loading {in_path}")
    adata = sc.read_h5ad(in_path)

    # make unique barcodes
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # 3) Inject metadata columns into .obs
    for col in ["gender", "batch", "condition"]:
        adata.obs[col] = row[col]
    
    # Add dataset/sample info
    adata.obs["dataset"] = dataset
    adata.obs["sample_id"] = sample
    adata.obs["group"] = adata.obs["gender"] + "_" + adata.obs["condition"]
    
    # Convert problematic columns to string if they exist
    if "nCount_RNA" in adata.obs:
        adata.obs["nCount_RNA"] = adata.obs["nCount_RNA"].astype(str)

    adatas.append(adata)

# 4) Merge all together
if len(adatas) == 0:
    raise RuntimeError("No AnnData objects loaded – check your paths/CSV")

print(f"Merging {len(adatas)} datasets")
merged = sc.concat(
    adatas,
    join="outer",     # keep all genes
    merge="same",     # require same .obs cols
    fill_value=0,     # fill missing genes with 0
    index_unique=None # keep original cell barcodes
)

# Ensure consistent data types in obs columns
for col in merged.obs.columns:
    try:
        # First try to convert to numeric (for columns like counts)
        merged.obs[col] = pd.to_numeric(merged.obs[col])
    except:
        # If that fails, convert to string
        merged.obs[col] = merged.obs[col].astype(str)

print(f"Result: {merged.n_obs} cells × {merged.n_vars} genes")

# Save with compression to reduce file size
merged.write("merged_all.h5ad", compression="gzip")