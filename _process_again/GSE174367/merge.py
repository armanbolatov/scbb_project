import os
import pandas as pd
import scanpy as sc

# Paths to your files
h5_file   = "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
meta_file = "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/GSE174367_snRNA-seq_cell_meta.csv"
ct_file   = "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/GSE174367_snRNA_cell_types.csv"
out_dir   = "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/ann_data_per_sample"
os.makedirs(out_dir, exist_ok=True)

# 1) Load the raw counts into an AnnData
adata = sc.read_10x_h5(h5_file)

# 2) Load and merge the author metadata
meta = pd.read_csv(meta_file)
meta = meta.set_index("Barcode")
adata.obs = adata.obs.join(meta, how="left")

# 3) Load and merge your cell-type calls
ct = pd.read_csv(ct_file)
ct = ct.set_index("barcode")
adata.obs["cell_type_marker"] = adata.obs_names.map(ct["cell_type"].to_dict())

# 4) Ensure all obs columns are strings (no NaNs)
adata.obs = adata.obs.fillna("").astype(str)

# 5) Split by SampleID and save each subset
for sample in adata.obs["SampleID"].unique():
    ad_sub = adata[adata.obs["SampleID"] == sample].copy()
    out_path = os.path.join(out_dir, f"{sample}.h5ad")
    ad_sub.write_h5ad(out_path)
    print(f"Wrote {out_path} ({ad_sub.n_obs} cells × {ad_sub.n_vars} genes)")
