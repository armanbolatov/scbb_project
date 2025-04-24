import os
import glob
import scanpy as sc

# Directories to inspect
ann_dirs = [
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE147528/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE157827/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE160936/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE174367/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE181279/ann_data",
]

# Loop through each directory, pick the first .h5ad, and summarize
for d in ann_dirs:
    files = glob.glob(os.path.join(d, "*.h5ad"))
    if not files:
        print(f"No .h5ad files found in {d}")
        continue
    f = files[0]  # take the first file
    adata = sc.read_h5ad(f)
    print(f"Directory: {d}")
    print(f"File: {os.path.basename(f)}")
    print(f"Shape (cells × genes): {adata.shape}")
    print("obs columns:", list(adata.obs_keys()))
    print("var columns:", list(adata.var_keys()))
    print("First 10 var_names:", adata.var_names[:10].tolist())
    print("-" * 60)

