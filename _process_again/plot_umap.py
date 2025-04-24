#!/usr/bin/env python3
import os, glob
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

ann_dirs = [
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE147528/ann_data",
    # "/l/users/darya.taratynova/scbb_project/_process_again/GSE157827/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE160936/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE174367/ann_data",
    # "/l/users/darya.taratynova/scbb_project/_process_again/GSE181279/ann_data",
]
out_png = "combined_umap.png"
h5ad_files = []
for d in ann_dirs:
    h5ad_files += glob.glob(os.path.join(d, "*.h5ad"))
if not h5ad_files:
    raise RuntimeError("No .h5ad files found in your ann_data folders!")

print(f"Found {len(h5ad_files)} AnnData files:\n")
print(f"{'File':<60} {'Cells':>8} {'Genes':>8} {'Approx MB':>12}")
print("-"*90)

total_cells = 0
total_bytes = 0
for fn in h5ad_files:
    ad = sc.read_h5ad(fn, backed="r")
    n_obs, n_vars = ad.shape

    try:
        X = ad.X
        if hasattr(X, "data"):
            nnz = X.data.shape[0]
            bytes_size = nnz * np.dtype(X.data.dtype).itemsize
        else:
            bytes_size = n_obs * n_vars * 8
    except Exception:
        bytes_size = n_obs * n_vars * 8

    mb = bytes_size/1e6
    print(f"{os.path.basename(fn):<60} {n_obs:8,} {n_vars:8,} {mb:12.1f}")
    total_cells += n_obs
    total_bytes += bytes_size

    # close if backed
    try: ad.file.close()
    except: pass

print("-"*90)
print(f"{'TOTAL':<60} {total_cells:8,} {'—':>8} {total_bytes/1e6:12.1f} MB\n")

adatas = []
for fn in h5ad_files:
    ad = sc.read_h5ad(fn)
    ad.var_names_make_unique()
    series = os.path.basename(os.path.dirname(os.path.dirname(fn)))
    ad.obs["dataset"] = series
    adatas.append(ad)

print(f"Loaded {len(adatas)} AnnData objects into memory → concatenating…")

combined = sc.concat(
    adatas,
    join="outer",    
)

sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
sc.pp.highly_variable_genes(
    combined, n_top_genes=2000, flavor="seurat", subset=True
)
sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, svd_solver="arpack")
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=20)
sc.tl.umap(combined)

fig = sc.pl.umap(
    combined,
    color=["dataset","cell_type_marker"],
    wspace=0.4,
    show=False
)
plt.savefig(out_png, dpi=300, bbox_inches="tight")
