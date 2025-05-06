import scanpy as sc
from collections import defaultdict
import os, glob

input_folder  = '/l/users/darya.taratynova/scbb_project/_process_again/ann_data_new'
output_folder = '/l/users/darya.taratynova/scbb_project/_process_again/ann_data_collapsed'
os.makedirs(output_folder, exist_ok=True)

for h5ad_path in glob.glob(os.path.join(input_folder, '*.h5ad')):
    # 1. load
    adata = sc.read_h5ad(h5ad_path)

    # 2. build ens→genes map
    ens_to_genes = defaultdict(list)
    for gene, ens in adata.var['ensembl_id'].dropna().items():
        ens_to_genes[ens].append(gene)

    # 3. build gene→primary map
    gene_to_primary = {}
    for ens, genes in ens_to_genes.items():
        primary = genes[0]            # the “first” name in each group
        for g in genes:
            gene_to_primary[g] = primary

    # 4. rename var_names
    new_var_names = [
        gene_to_primary.get(g, g)    # default to itself if no ens or unique
        for g in adata.var_names
    ]
    adata.var_names = new_var_names

    # 5. (opt) drop duplicate var entries, keeping the first
    adata = adata[:, ~adata.var_names.duplicated()]

    # 6. save out
    base = os.path.basename(h5ad_path)
    out_path = os.path.join(output_folder, base)
    adata.write_h5ad(out_path)
    print(f"Wrote collapsed AnnData to {out_path}")
