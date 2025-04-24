import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

upregulated_genes = ["ACSL1", "CREM", "PADI2", "SRGN", "USP9X"]
downregulated_genes = ["APBB1IP", "BCYRN1", "GLDN", "LHFPL2", "LRRK1", "RASSF8", "RTTN", "XIST" ]

# genes_from_papers = ["KLHL21", "WDR82", "DTX3L", "UBTD2", "CISH", "ATXN3L", "DGKG", "MAP3K7IP2", "NFKBIE", "VIP", "PCCB"]
genes_from_papers = ["KLHL21", "WDR82", "PCCB"]
# 

# sample IDs from folder file extract sample IDs
dataset_paths = os.listdir("/l/users/shahad.hardan/filtered_data/")
dataset_paths = [f"/l/users/shahad.hardan/filtered_data/{path}" for path in dataset_paths]

ids = [path.split("/")[-1].split(".")[0].split("_")[-1] for path in dataset_paths]

metadata = pd.read_csv("/home/shahad.hardan/CB803/data/metadata.csv")

var_dict = np.load("/home/shahad.hardan/CB803/data/var_metadata.npy", allow_pickle=True).item()
obs_dict = np.load("/home/shahad.hardan/CB803/data/obs_metadata.npy", allow_pickle=True).item()

# Step 2: Extract Expression Values for Each Gene
for gene in genes_from_papers:
    expression_data = []

    for sample_id in ids:
        print(f"Processing sample {sample_id}...")
        # Check if X file exists in npz or npy format
        X_file_npz = f"/l/users/shahad.hardan/filtered_data/X_data_{sample_id}.npz"
        X_file_npy = f"/l/users/shahad.hardan/filtered_data/X_data_{sample_id}.npy"
        if os.path.exists(X_file_npz):
            X = sparse.load_npz(X_file_npz)
        elif os.path.exists(X_file_npy):
            X = np.load(X_file_npy, allow_pickle=True)
        else:
            print(f"Warning: X file not found for sample {sample_id}. Skipping...")
            continue

        var = var_dict["dataset_" + sample_id]
        obs = obs_dict["dataset_" + sample_id] 
        if "prediction.score.Microglia" in obs.columns:
            obs_filtered = obs[obs["prediction.score.Microglia"] == 1] # Filter for Microglia cells
            obs_indices = np.where(obs.index.isin(obs_filtered.index))[0]
            print(f"Microglia prediction score found in obs for sample {sample_id}.")
            gene_expr = X[obs_indices, :]  # Get expression values for Microglia cells
        else: 
            print(f"Warning: Microglia prediction score not found in obs for sample {sample_id}. Skipping...")
            continue

        # Check if gene exists in the current sample
        if gene in var.index:
            gene_idx = np.where(var.index == gene)[0][0]  # Find gene index
            print(f"Gene {gene} found in sample {sample_id}")
            gene_expr = gene_expr[:, gene_idx].toarray().flatten() if sparse.issparse(X) else X[:, gene_idx]  # Get expression values

            # Store expression data
            expression_data.append({
                "Gene": gene,
                "Expression": gene_expr,
                "SampleID": sample_id,
                "Gender": metadata.loc[metadata["Dataset"] == sample_id, "Gender"].iloc[0]  # Get Male/Female label
            })


    df = pd.DataFrame(expression_data)
    df["Gender"] = df["Gender"].astype(str) 
    df["Gender"] = pd.Categorical(df["Gender"], categories=["Male", "Female"])  
    df = df.explode("Expression")

    df["Expression"] = pd.to_numeric(df["Expression"], errors="coerce")
    # Apply log transformation
    df["Expression"] = df["Expression"].apply(lambda x: np.log1p(x) if x >= 0 else np.nan)

    # Step 3: Generate Violin Plot for This Gene
    plt.figure(figsize=(6, 5))
    sns.violinplot(data=df.dropna(subset=["Expression"]), x="Gender", y="Expression", palette={"Male": "blue", "Female": "red"})

    plt.title(f"Violin Plot for {gene} in Microglia (Male vs. Female)")
    plt.xlabel("Gender")
    plt.ylabel("Expression Level")
    plt.savefig(f"violin_plot_{gene}.png") 
    plt.show()
