import numpy as np
from scipy import sparse

# Load the NPZ file and metadata
npz_path = "/l/users/darya.taratynova/scbb_project/filtered_data_2/X_data_AD06301.npz" 
var_metadata_path = "/l/users/darya.taratynova/scbb_project/var_metadata.npy"
obs_metadata_path = "/l/users/darya.taratynova/scbb_project/obs_metadata.npy"

loaded_npz = np.load(npz_path)

# Print keys in NPZ file
print("Keys in NPZ file:")
print(list(loaded_npz.keys()))

# Load metadata
var_metadata = np.load(var_metadata_path, allow_pickle=True).item()
obs_metadata = np.load(obs_metadata_path, allow_pickle=True).item()

# Get dataset key
dataset_key = 'dataset_AD05701'

if dataset_key in var_metadata:
    df_var = var_metadata[dataset_key]
    print(f"Metadata for {dataset_key} (features):")
    print(df_var.head())

# Load X_data
X_data = loaded_npz['data']  
X_data = sparse.csr_matrix((X_data, loaded_npz['indices'], loaded_npz['indptr']), shape=loaded_npz['shape'])

# Define scavenger receptor genes
scavenger_receptor_genes = [
    'SCARA1', 'SCARA2', 'MARCO', 'CD36', 'SCARB1', 'CD68',
    'OLR1', 'SCARF1', 'MEGF10', 'RAGE', 'CD163', 'CXCL16'
]

# Get the first 10 genes from df_var index
first_10_genes = df_var.index[:10]
print(f"First 10 genes in the dataset: {first_10_genes}")
# print(df_var)
# Check if any scavenger receptor genes are in the first 10 genes
for gene in scavenger_receptor_genes:
    if gene in df_var:
        print(f"{gene} is present!")
    else:
        print(f"{gene} is not present.")
