import numpy as np
import os
import pandas as pd

X_dir = "/l/users/shahad.hardan/filtered_data/"
obs_path = "obs_metadata.npy"
var_path = "var_metadata.npy"

obs_dict = np.load(obs_path, allow_pickle=True).item()
var_dict = np.load(var_path, allow_pickle=True).item()

X_files = sorted(os.listdir(X_dir))

empty_datasets = []

# Check obs and var
for dataset_id in obs_dict.keys():
    obs_df = obs_dict[dataset_id]
    var_df = var_dict.get(dataset_id, None)  

    if obs_df.empty:
        empty_datasets.append({"Dataset": dataset_id, "Type": "obs", "Issue": "Empty DataFrame"})
    
    if var_df is not None and var_df.empty:
        empty_datasets.append({"Dataset": dataset_id, "Type": "var", "Issue": "Empty DataFrame"})

from scipy import sparse
import numpy as np

for X_file in X_files:
    dataset_id = X_file.replace("X_data_", "").replace(".npy", "").replace(".npz", "")  # Extract dataset ID
    file_path = os.path.join(X_dir, X_file)

    if X_file.endswith(".npz"):
        print(f"Checking sparse matrix: {X_file}")
        X_matrix = sparse.load_npz(file_path)
    elif X_file.endswith(".npy"):
        print(f"Checking dense matrix: {X_file}")
        X_matrix = np.load(file_path, allow_pickle=True)
    else:
        print(f"⚠️ Unknown file format: {X_file}, skipping.")
        continue

    if sparse.issparse(X_matrix):
        if X_matrix.nnz == 0:  
            empty_datasets.append({"Dataset": dataset_id, "Type": "X", "Issue": "Empty Matrix (Sparse)"})
    elif X_matrix.ndim == 0 or X_matrix.size == 0:  
        empty_datasets.append({"Dataset": dataset_id, "Type": "X", "Issue": "Empty Matrix (Dense)"})


    # Check if X_matrix has no rows or columns
    elif len(X_matrix.shape) < 2 or X_matrix.shape[0] == 0 or X_matrix.shape[1] == 0:
        empty_datasets.append({"Dataset": dataset_id, "Type": "X", "Issue": "No Rows/Cols"})

    # Check if X_matrix contains only zeros 
    elif sparse.issparse(X_matrix):
        if X_matrix.nnz == 0:  #
            empty_datasets.append({"Dataset": dataset_id, "Type": "X", "Issue": "All Zeros (Sparse)"})
    else:
        if np.all(X_matrix == 0):  
            empty_datasets.append({"Dataset": dataset_id, "Type": "X", "Issue": "All Zeros (Dense)"})


empty_df = pd.DataFrame(empty_datasets)
empty_df.to_csv("empty_datasets.csv", index=False)

print("Check complete. Empty datasets saved to 'empty_datasets.csv'.")
