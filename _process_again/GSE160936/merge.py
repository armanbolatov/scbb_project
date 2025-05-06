import os
import glob
import scanpy as sc
import pandas as pd

# Base directories
base_dir = "/l/users/darya.taratynova/scbb_project/GSE160936/filtered_matrices"
results_dir = os.path.join(base_dir, "results")
out_dir = "/l/users/darya.taratynova/scbb_project/GSE160936/filtered_matrices/ann_data_per_sample"
os.makedirs(out_dir, exist_ok=True)

# 1) Find all sample IDs by looking at results CSV files
ct_files = glob.glob(os.path.join(results_dir, "*_cell_types.csv"))
samples = [os.path.basename(f).replace("_cell_types.csv", "") for f in ct_files]

for s in samples:
    print(f"Processing sample {s}...")
    # Paths to 10x flat files
    mat_dir = os.path.join(base_dir, s, "rds")
    # Read in count matrix
    adata = sc.read_10x_mtx(mat_dir, var_names="gene_symbols", cache=False)
    
    # Load full metadata
    meta_full = pd.read_csv(os.path.join(results_dir, f"{s}_metadata_full.csv"), index_col=0)
    # Load cell type calls
    ct = pd.read_csv(os.path.join(results_dir, f"{s}_cell_types.csv"), index_col="barcode")
    
    # Align obs indices
    adata.obs = adata.obs.join(meta_full, how="left")
    adata.obs["cell_type_marker"] = adata.obs_names.map(ct["cell_type"])
    
    # Ensure obs columns are strings
    adata.obs = adata.obs.fillna("").astype(str)
    
    # Save per-sample AnnData
    out_f = os.path.join(out_dir, f"{s}.h5ad")
    adata.write_h5ad(out_f)
    print(f"  → wrote {out_f} ({adata.n_obs} cells × {adata.n_vars} genes)")
