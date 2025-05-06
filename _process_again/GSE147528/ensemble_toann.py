import os
import scanpy as sc
import pandas as pd
from scipy import io

# directory where your flat files live
RESULTS_DIR = "results"

for fname in os.listdir(RESULTS_DIR):
    # pick up each counts.mtx
    if not fname.endswith("_counts.mtx"):
        continue

    base = fname.replace("_counts.mtx", "")
    mtx_path   = os.path.join(RESULTS_DIR, f"{base}_counts.mtx")
    genes_path = os.path.join(RESULTS_DIR, f"{base}_genes.tsv")
    bc_path    = os.path.join(RESULTS_DIR, f"{base}_barcodes.tsv")
    meta_path  = os.path.join(RESULTS_DIR, f"{base}_metadata.csv")

    print(f"Reassembling AnnData for sample '{base}'…")
    # 1) Read the count matrix (genes × cells) and transpose to cells × genes
    X = io.mmread(mtx_path).T.tocsc()

    # 2) Read gene & barcode labels
    genes    = pd.read_csv(genes_path, header=None)[0].values
    barcodes = pd.read_csv(bc_path,    header=None)[0].values

    # 3) Read metadata (must have row index matching barcodes)
    meta = pd.read_csv(meta_path, index_col=0)

    # 4) Read UMAP coords (row index = barcodes)

    # 5) Build AnnData
    adata = sc.AnnData(X)
    adata.var_names = genes
    adata.obs_names = barcodes

    # 6) Attach metadata and embeddings
    #    Align obs to ensure same order
    adata.obs = meta.loc[adata.obs_names]

    # 7) Save final .h5ad
    out_file = os.path.join("results_anndata", f"{base}.h5ad")
    adata.write_h5ad(out_file)
    print(f"  → Written {out_file}")