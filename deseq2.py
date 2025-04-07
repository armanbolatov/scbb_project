import os
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

data_dir = "/l/users/darya.taratynova/scbb_project/filtered_data"
var_metadata_path = "/l/users/darya.taratynova/scbb_project/var_metadata.npy"
obs_metadata_path = "/l/users/darya.taratynova/scbb_project/obs_metadata.npy"
meta_csv_path = "/l/users/darya.taratynova/scbb_project/all_grouped.csv"

var_metadata_dict = np.load(var_metadata_path, allow_pickle=True).item()
obs_metadata_dict = np.load(obs_metadata_path, allow_pickle=True).item()

meta_info = pd.read_csv(meta_csv_path)

all_dfs = []
i = 0

for fname in os.listdir(data_dir):
    if not fname.startswith("X_data_AD") or not (fname.endswith(".npz") or fname.endswith(".npy")):
        continue

    dataset_id = fname.replace("X_data_", "").replace(".npz", "").replace(".npy", "")
    dataset_key = f"dataset_{dataset_id}"
    fpath = os.path.join(data_dir, fname)

    print(f"Loading ({i}): {fname}")
    i += 1

    if fname.endswith(".npz"):
        loader = np.load(fpath)
        if {'indices', 'indptr', 'format', 'shape', 'data'}.issubset(loader.keys()) and b'csr' in loader['format']:
            X_data = sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])
        else:
            print(f"Skipping {fname} (invalid format)")
            continue
    elif fname.endswith(".npy"):
        X_data = sparse.csr_matrix(np.load(fpath))

    if dataset_key not in var_metadata_dict or dataset_key not in obs_metadata_dict:
        print(f"Metadata missing for {dataset_id}, skipping.")
        continue

    df_var = var_metadata_dict[dataset_key]
    df_obs = obs_metadata_dict[dataset_key]

    if df_obs.index.isnull().any():
        df_obs.index = [f"{idx}_{dataset_id}" for idx in range(X_data.shape[0])]
    if df_var.index.isnull().any():
        df_var.index = [f"gene_{idx}" for idx in range(X_data.shape[1])]

    df_X = pd.DataFrame.sparse.from_spmatrix(X_data, index=df_obs.index, columns=df_var.index)
    df_combined = pd.concat([df_obs, df_X], axis=1)
    df_combined['dataset'] = dataset_id  

    all_dfs.append(df_combined)

df_full_sparse = pd.concat(all_dfs, axis=0)
print(f"\nMerged shape: {df_full_sparse.shape}")

df_full_sparse = df_full_sparse.merge(meta_info[['Dataset', 'Gender', 'Condition', 'Age_Group']],
                                      left_on='dataset', right_on='Dataset', how='left')
df_full_sparse.drop(columns=['Dataset'], inplace=True)

cell_type = 'prediction.score.Microglia'
selected_cells = df_full_sparse[df_full_sparse[cell_type] == 1.0]
grouped = selected_cells.groupby('dataset')

pseudo_expr = []
meta_rows = []
dataset_ids = []

gene_columns = df_var.index

for dataset_id, group in grouped:
    X_sparse = sparse.csr_matrix(group[gene_columns].sparse.to_coo())
    expr_sum = pd.Series(np.array(X_sparse.sum(axis=0)).flatten(), index=gene_columns)
    if (expr_sum < 0).any():
        print(f"Negative values in dataset {dataset_id}")
        print(expr_sum[expr_sum < 0].head())

    meta = group.iloc[0][['Gender', 'Age_Group', 'Condition']]  
    pseudo_expr.append(expr_sum)
    meta_rows.append(meta)
    dataset_ids.append(str(dataset_id))  
    print(f"Processed {dataset_id}")

expr_df = pd.DataFrame(pseudo_expr, index=dataset_ids).clip(lower=0).astype(int)
meta_df = pd.DataFrame(meta_rows, index=dataset_ids)

meta_df['Gender'] = meta_df['Gender'].astype('category')
meta_df['Condition'] = meta_df['Condition'].astype('category')

expr_df.index = expr_df.index.astype(str)
meta_df.index = meta_df.index.astype(str)

common_ids = expr_df.index.intersection(meta_df.index)
expr_df = expr_df.loc[common_ids].sort_index()
meta_df = meta_df.loc[common_ids].sort_index()

print("Building DESeq2 dataset...")
dds = DeseqDataSet(
    counts=expr_df,
    metadata=meta_df,
    design="~ Gender + Age_Group + Condition",  
    refit_cooks=True
)

dds.deseq2()
stat_res = DeseqStats(dds, contrast=["Condition", "Disease", "Control"])
stat_res.summary()
results_df = stat_res.results_df.sort_values("padj").dropna()
results_df = results_df.reset_index() 
print(results_df.head())
results_df.to_csv("deseq2_results_microglia.csv", index=False)
