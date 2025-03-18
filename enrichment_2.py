import os
import numpy as np
import pandas as pd
from scipy import sparse
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def gsea_plot_for_comparison(df_full_sparse, df_var, gene_sets, 
                              gene_set_name, 
                              cell_type_score_column,
                              group1, group2, 
                              output_filename):
    """
    Perform GSEA between two groups for a specific cell type and gene set.
    Parameters:
    - df_full_sparse: merged full DataFrame
    - df_var: variable metadata for gene names
    - gene_sets: dictionary of gene sets (from gseapy.get_library)
    - gene_set_name: e.g., 'neurogenesis (GO:0022008)'
    - cell_type_score_column: column to filter e.g., 'prediction.score.Microglia'
    - group1/group2: dicts with 'Gender' and 'Condition'
    - output_filename: for saving the plot
    """

    selected_cells = df_full_sparse[df_full_sparse[cell_type_score_column] == 1.0]

    group1_df = selected_cells[(selected_cells['Gender'] == group1['Gender']) & 
                               (selected_cells['Condition'] == group1['Condition'])]
    group2_df = selected_cells[(selected_cells['Gender'] == group2['Gender']) & 
                               (selected_cells['Condition'] == group2['Condition'])]

    print(f"{group1['Gender']} {group1['Condition']} {cell_type_score_column}: {group1_df.shape[0]}, "
          f"{group2['Gender']} {group2['Condition']} {cell_type_score_column}: {group2_df.shape[0]}")

    gene_columns = df_var.index
    mean1 = group1_df[gene_columns].sparse.to_dense().mean()
    mean2 = group2_df[gene_columns].sparse.to_dense().mean()

    logfc = np.log2(mean1 + 1) - np.log2(mean2 + 1)
    ranking = logfc.sort_values(ascending=False)
    rnk_df = pd.DataFrame(ranking).reset_index()
    rnk_df.columns = ['Gene', 'logFC']

    rnk_df['Gene'] = rnk_df['Gene'].str.upper()
    geneset_genes = gene_sets.get(gene_set_name)

    print(f"Gene count in '{gene_set_name}': {len(geneset_genes)}")

    gsea_res = gp.prerank(rnk=rnk_df,
                          gene_sets={gene_set_name: geneset_genes},
                          permutation_num=100,
                          seed=42,
                          outdir=None,
                          verbose=True)

    res = gsea_res.results[gene_set_name]
    RES = res['RES']
    hits = res['hits']

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(range(len(RES)), RES, color='green', linewidth=1.5)
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
    for hit in hits:
        ax.axvline(hit, ymin=0.95, ymax=1.0, color='black', linewidth=0.5)

    ax.set_xlabel('rank', fontsize=10)
    ax.set_ylabel('enrichment score', fontsize=10)
    title_text = f"{gene_set_name.split()[0]}\n{group1['Gender']} {group1['Condition']} vs {group2['Gender']} {group2['Condition']} in {cell_type_score_column.split('.')[-1]}"
    ax.set_title(title_text, fontsize=12, loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, len(RES))
    ax.set_ylim(min(RES) - 0.1, max(RES) + 0.1)
    plt.tight_layout()

    fig.savefig(f"{output_filename}.png", dpi=300)
    fig.savefig(f"{output_filename}.pdf", bbox_inches='tight')
    plt.close(fig)


data_dir = "/l/users/darya.taratynova/scbb_project/filtered_data_2/"
var_metadata_path = "/l/users/darya.taratynova/scbb_project/var_metadata.npy"
obs_metadata_path = "/l/users/darya.taratynova/scbb_project/obs_metadata.npy"

var_metadata_dict = np.load(var_metadata_path, allow_pickle=True).item()
obs_metadata_dict = np.load(obs_metadata_path, allow_pickle=True).item()

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

    # --- Load the matrix ---
    if fname.endswith(".npz"):
        loader = np.load(fpath)
        if {'indices', 'indptr', 'format', 'shape', 'data'}.issubset(set(loader.keys())) and b'csr' in loader['format']:
            X_data = sparse.csr_matrix((loader['data'], loader['indices'], loader['indptr']), shape=loader['shape'])
        else:
            print(f"Skipping {fname} (invalid format)")
            continue
    elif fname.endswith(".npy"):
        X_data = np.load(fpath)
        X_data = sparse.csr_matrix(X_data)  

    # --- Load metadata ---
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
meta_info = pd.read_csv("/l/users/darya.taratynova/scbb_project/all.csv")

df_full_sparse = df_full_sparse.merge(meta_info[['Dataset', 'Gender', 'Condition']],
                                      left_on='dataset', right_on='Dataset', how='left')
df_full_sparse.drop(columns=['Dataset'], inplace=True)

gene_sets = gp.get_library(name='GO_Biological_Process_2021')
# gsea_plot_for_comparison(df_full_sparse, df_var, gene_sets,
#                          gene_set_name='neurogenesis (GO:0022008)',
#                          cell_type_score_column='prediction.score.Microglia',
#                          group1={'Gender': 'Male', 'Condition': 'Disease'},
#                          group2={'Gender': 'Male', 'Condition': 'Control'},
#                          output_filename='gsea_male_ad_vs_control_microglia')
# gsea_plot_for_comparison(df_full_sparse, df_var, gene_sets,
#                          gene_set_name='neurogenesis (GO:0022008)',
#                          cell_type_score_column='prediction.score.Microglia',
#                          group1={'Gender': 'Female', 'Condition': 'Disease'},
#                          group2={'Gender': 'Female', 'Condition': 'Control'},
#                          output_filename='gsea_female_ad_vs_control_microglia')
# gsea_plot_for_comparison(df_full_sparse, df_var, gene_sets,
#                          gene_set_name='neurogenesis (GO:0022008)',
#                          cell_type_score_column='prediction.score.Microglia',
#                          group1={'Gender': 'Male', 'Condition': 'Disease'},
#                          group2={'Gender': 'Female', 'Condition': 'Disease'},
#                          output_filename='gsea_male_vs_female_ad_microglia')
scavenger_receptor_genes = [
    'SCARA1', 'SCARA2', 'MARCO', 'CD36', 'SCARB1', 'CD68',
    'OLR1', 'SCARF1', 'MEGF10', 'RAGE', 'CD163', 'CXCL16'
]

scavenger_gene_set = {'Scavenger_Receptors': scavenger_receptor_genes}

gsea_plot_for_comparison(df_full_sparse, df_var, scavenger_gene_set,
                         gene_set_name='Scavenger_Receptors',
                         cell_type_score_column='prediction.score.Microglia',
                         group1={'Gender': 'Male', 'Condition': 'Disease'},
                         group2={'Gender': 'Female', 'Condition': 'Disease'},
                         output_filename='gsea_male_vs_female_ad_scavenger_microglia')
