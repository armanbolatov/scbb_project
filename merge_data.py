import anndata as ad
import numpy as np
import os
import pandas as pd
from scipy import sparse
from tqdm import tqdm

#Find common genes across datasets
print("Finding common genes...")
dataset_paths = os.listdir("/l/users/shahad.hardan/ann_data/")
dataset_paths = [f"/l/users/shahad.hardan/ann_data/{path}" for path in dataset_paths]

ids = [path.split("/")[-1].split(".")[0] for path in dataset_paths]

# common_genes = pd.read_csv("common_genes_1.csv", header=None).squeeze("columns").astype(str).str.strip().tolist()


all_genes = []
for path in dataset_paths:
    adata = ad.read_h5ad(path)

    # Some samples have the indices stored in a different column
    if "_index" in adata.var.columns:
        gene_names = adata.var["_index"].astype(str)
    else:
        gene_names = adata.var.index.astype(str)

    all_genes.append(set(gene_names))


common_genes = set.intersection(*all_genes)
common_genes = list(map(str, common_genes))
print(f"Number of common genes: {len(common_genes)}")

# Save common genes
# pd.Series(common_genes).to_csv("common_genes_1.csv", index=False)

# Step 2: Save obs matrices
print("Saving obs metadata...")
obs_dict = {}
for i, path in tqdm(enumerate(dataset_paths), total=len(dataset_paths)):
    adata = ad.read_h5ad(path)

    # Some samples have the indices stored in a different column
    if "_index" in adata.var.columns:
        adata.var.index = adata.var["_index"].astype(str)

    adata.var.index = adata.var.index.astype(str).str.strip()
    common_genes_filtered = [gene.strip() for gene in common_genes if gene in adata.var.index]

    adata = adata[:, common_genes_filtered].copy()
    
    obs_dict[f"dataset_{ids[i]}"] = adata.obs.copy()

np.save("obs_metadata.npy", obs_dict)


# Step 3: Save var matrix 
print("Saving var metadata...")
var_dict = {}
for i, path in tqdm(enumerate(dataset_paths), total=len(dataset_paths)):
    adata = ad.read_h5ad(path)

    if "_index" in adata.var.columns:
        adata.var.index = adata.var["_index"].astype(str)

    adata.var.index = adata.var.index.astype(str).str.strip()
    common_genes_filtered = [gene.strip() for gene in common_genes if gene in adata.var.index]

    adata = adata[:, common_genes_filtered].copy()

    var_dict[f"dataset_{ids[i]}"] = adata.var.copy()

# Save var metadata
np.save("var_metadata.npy", var_dict)

#take again all ids in empty_datasets.csv and re run 
# Step 1: Load empty datasets
empty_df = pd.read_csv("empty_datasets.csv")
empty_datasets = empty_df["Dataset"].unique()
print(f"Total empty datasets: {len(empty_datasets)}")

# Step 4: Saving X data according to the common genes and file format
print("Saving X data...")
for i, path in enumerate(dataset_paths):
    print(i)
    print(path)
    print(f"\nProcessing dataset: {ids[i]}")
    if ids[i] in empty_datasets:
    
        adata = ad.read_h5ad(path)

        if "_index" in adata.var.columns:
            adata.var.index = adata.var["_index"].astype(str)

        adata.var.index = adata.var.index.astype(str).str.strip()
        print(f"Total genes before filtering: {adata.shape[1]}")
        
        if adata.X.shape[0] == 0:
            print(f"⚠️ Warning: Dataset {ids[i]} already has an empty X matrix before filtering!")

        # Ensure only valid genes are used
        common_genes_filtered = [gene for gene in common_genes if gene in adata.var.index]
        print(f"Matching genes found: {len(common_genes_filtered)} / {len(common_genes)}")

        if len(common_genes_filtered) == 0:
            print(f"No common genes found in dataset {ids[i]}!")

        # Subset using clean common_genes
        adata = adata[:, common_genes_filtered].copy()

        print(f"Genes after filtering: {adata.shape[1]}")

        # Check if adata.X is sparse
        if sparse.issparse(adata.X):
            print(f"Saving dataset {ids[i]} as sparse .npz format.")
            sparse.save_npz(f"/l/users/shahad.hardan/filtered_data/X_data_{ids[i]}.npz", adata.X)
        else:
            print(f"Saving dataset {ids[i]} as dense .npy format.")
            np.save(f"/l/users/shahad.hardan/filtered_data/X_data_{ids[i]}.npy", adata.X)

print("Done!")
