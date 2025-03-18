# import numpy as np
# import pandas as pd

# # File paths
# npz_path = "/l/users/darya.taratynova/scbb_project/filtered_data/X_data_AD05701.npz"
# var_metadata_path = "/l/users/darya.taratynova/scbb_project/var_metadata.npy"
# obs_metadata_path = "/l/users/darya.taratynova/scbb_project/obs_metadata.npy"

# loaded_npz = np.load(npz_path)

# # List keys (arrays stored inside)
# print("Keys in NPZ file:")
# print(list(loaded_npz.keys()))

# # Inspect shapes/types of each key
# for key in loaded_npz.files:
#     array = loaded_npz[key]
#     print(f"\nKey: {key}")
#     print(f"  Shape: {array.shape}")
#     print(f"  Data type: {array.dtype}")
#     print(f"  Sample values: {array.ravel()[:5]}")
# # Load metadata
# var_metadata = np.load(var_metadata_path, allow_pickle=True).item()
# obs_metadata = np.load(obs_metadata_path, allow_pickle=True).item()

# # Inspect keys in var_metadata and obs_metadata
# print("var_metadata keys (datasets):")
# print(list(var_metadata.keys())[:5], "...")  
# print(f"Total datasets in var_metadata: {len(var_metadata)}\n")

# print("obs_metadata keys (datasets):")
# print(list(obs_metadata.keys())[:5], "...")  
# print(f"Total datasets in obs_metadata: {len(obs_metadata)}\n")

# # Dataset key for X_data
# dataset_key = 'dataset_AD05701'

# # --- var_metadata ---
# if dataset_key in var_metadata:
#     df_var = var_metadata[dataset_key]
#     print(f"Metadata for {dataset_key} (features):")
#     print(df_var.head())
#     print(f"Feature count: {df_var.shape[0]}, Columns: {df_var.shape[1]}\n")
# else:
#     print(f"{dataset_key} not found in var_metadata")

# # --- obs_metadata ---
# if dataset_key in obs_metadata:
#     df_obs = obs_metadata[dataset_key]
#     print(f"Metadata for {dataset_key} (cells):")
#     print(df_obs.head())
#     print(f"Cell count: {df_obs.shape[0]}, Columns: {df_obs.shape[1]}\n")
# else:
#     print(f"{dataset_key} not found in obs_metadata")

# # --- Align and create full DataFrame ---
# # Confirm matching dimensions
# if df_obs.shape[0] == X_data.shape[0] and df_var.shape[0] == X_data.shape[1]:
#     # Set index and column names if not already set
#     if df_obs.index.isnull().any():
#         df_obs.index = [f"cell_{i}" for i in range(df_obs.shape[0])]
#     if df_var.index.isnull().any():
#         df_var.index = [f"gene_{i}" for i in range(df_var.shape[0])]

#     # Create DataFrame from X_data with appropriate index and columns
#     df_X = pd.DataFrame(X_data, index=df_obs.index, columns=df_var.index)
#     print("Aligned DataFrame df_X (expression matrix):")
#     print(df_X.head())

#     # Optionally, concatenate metadata for exploration
#     df_full = pd.concat([df_obs, df_X], axis=1)
#     print("\nFull DataFrame with metadata and expression data:")
#     print(df_full.head())

# else:
#     print("Mismatch in dimensions between X_data and metadata.")
#     print(f"X_data shape: {X_data.shape}")
#     print(f"df_obs shape: {df_obs.shape}")
#     print(f"df_var shape: {df_var.shape}")
import os
import numpy as np

filtered_data_dir = "/l/users/darya.taratynova/scbb_project/filtered_data/"
results = []

for filename in os.listdir(filtered_data_dir):
    filepath = os.path.join(filtered_data_dir, filename)

    if filename.endswith(".npy"):
        try:
            data = np.load(filepath, mmap_mode='r')  # mmap for big files
            dtype = data.dtype
            sample_vals = data.ravel()[:5]
            min_val = sample_vals.min()
            is_raw = (np.issubdtype(dtype, np.integer) or (min_val >= 0 and sample_vals.max() < 10000 and np.all(sample_vals == sample_vals.astype(int))))
            if is_raw:
                results.append((filename, "npy", f"Raw Count Matrix? dtype={dtype}, min={min_val}"))
            else:
                results.append((filename, "npy", f"Processed Matrix dtype={dtype}, min={min_val}"))
        except Exception as e:
            results.append((filename, "npy", f"Error reading: {str(e)}"))

    elif filename.endswith(".npz"):
        try:
            loader = np.load(filepath)
            keys = list(loader.keys())
            if {'indices', 'indptr', 'format', 'shape', 'data'}.issubset(set(keys)) and b'csr' in loader['format']:
                results.append((filename, "npz (CSR)", "Raw Count Matrix (sparse)"))
            else:
                results.append((filename, "npz", "Unknown content"))
        except Exception as e:
            results.append((filename, "npz", f"Error reading: {str(e)}"))

# Display results
print(f"{'File':<40} {'Type':<12} {'Content'}")
print("-" * 80)
for fname, ftype, content in sorted(results):
    print(f"{fname:<40} {ftype:<12} {content}")
