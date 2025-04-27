import os
import gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scvi
from scvi.external import MRVI
from scipy import sparse
import anndata as ad
from tqdm import tqdm

# Create output directories
os.makedirs("figs", exist_ok=True)
os.makedirs("results", exist_ok=True)

# Load the merged AnnData object
print("Loading merged AnnData file...")
adata = sc.read_h5ad("merged_all.h5ad")
print(f"\nDataset dimensions: {adata.n_obs} cells × {adata.n_vars} genes")

# Subsample for memory efficiency
SAMPLE_SIZE = 10000
sc.pp.subsample(adata, n_obs=SAMPLE_SIZE, random_state=42)
print(f"Subsampled adata to {adata.n_obs} cells and {adata.n_vars} genes")
gc.collect()

# Select highly variable genes
sc.pp.highly_variable_genes(
    adata, n_top_genes=1000, inplace=True, subset=True, flavor="seurat_v3"
)

# Setup and train model
MRVI.setup_anndata(
    adata,
    sample_key="sample_id",
    batch_key="batch",
    labels_key="group"
)

model = MRVI(adata)
model.train(max_epochs=1)

# Save model after training
print("Saving trained model...")
model.save("mrvi_trained", overwrite=True)
gc.collect()

# Plot and save training metrics
output_dir = "figs"
plt.figure(figsize=(8, 6))
plt.plot(model.history["elbo_validation"].iloc[5:])
plt.xlabel("Epoch")
plt.ylabel("Validation ELBO")
plt.savefig(os.path.join(output_dir, "mrvi_elbo.png"), dpi=300, bbox_inches="tight")
plt.close()

# Extract latent representation
print("Computing latent representations...")
u = model.get_latent_representation()
z = model.get_latent_representation(give_z=True)

# Get the h representation
print("Computing h representation...")
# Extract necessary indices
sample_indices = model.adata.obs.sample_id.cat.codes.values
batch_indices = model.adata.obs.batch.cat.codes.values
# Get the raw count data
x = model.adata.X.toarray() if sparse.issparse(model.adata.X) else model.adata.X
# Generate random epsilon for the compute_h function
np.random.seed(42)
extra_eps = np.random.normal(0, 1, size=(x.shape[0], model.module.n_latent))

# Access the compute_h_from_x_eps method directly from the module
h = model.module.compute_h_from_x_eps(
    x=x, 
    sample_index=sample_indices, 
    batch_index=batch_indices, 
    extra_eps=extra_eps
)

# Store representations in AnnData
adata.obsm["u"] = u
adata.obsm["z"] = z
adata.obsm["h"] = h

# Save representations for later use
representation_df = pd.DataFrame(
    index=adata.obs_names,
    data={
        **{f"u_{i}": u[:, i] for i in range(u.shape[1])},
        **{f"z_{i}": z[:, i] for i in range(z.shape[1])},
        **{f"h_{i}": h[:, i] for i in range(h.shape[1])}
    }
)
representation_df = pd.concat([representation_df, adata.obs], axis=1)
representation_df.to_csv("results/latent_representations.csv")
del u, z, h
gc.collect()

# Compute UMAP embedding
print("Computing UMAP embedding...")
sc.pp.neighbors(adata, use_rep="u")
sc.tl.umap(adata)
gc.collect()

# Format group names for plotting - remove underscores
print("Formatting group names for display...")
adata.obs['display_group'] = adata.obs['group'].astype(str).str.replace('_', ' ')

# Create UMAP visualizations for all latent representations
print("Creating UMAP visualizations for latent spaces...")
for rep_name, rep_key in [("u", "u"), ("z", "z"), ("h", "h")]:
    # Create neighbors and UMAP for this representation
    sc.pp.neighbors(adata, use_rep=rep_key, key_added=f"neighbors_{rep_name}")
    sc.tl.umap(adata, neighbors_key=f"neighbors_{rep_name}", min_dist=0.3, copy=False, key_added=f"X_umap_{rep_name}")
    
    # Store current UMAP temporarily in the main slot for plotting
    original_umap = None
    if "X_umap" in adata.obsm:
        original_umap = adata.obsm["X_umap"].copy()
    adata.obsm["X_umap"] = adata.obsm[f"X_umap_{rep_name}"]
    
    # Define variables to plot
    plot_vars = ['display_group', 'batch']
    
    # Save UMAP plots for each representation
    sc.pl.umap(
        adata,
        color=plot_vars,
        frameon=False,
        ncols=1,
        title=[f"{rep_name.upper()} Latent Space: {var.capitalize()}" for var in plot_vars],
        legend_loc='upper left',
        save=f"_{rep_name}_latent_space.pdf"
    )
    
    # Restore original UMAP if it existed
    if original_umap is not None:
        adata.obsm["X_umap"] = original_umap
    
    gc.collect()

# Create specific comparisons: F_Disease vs F_Control and M_Disease vs M_Control
print("\n=== Running Targeted Differential Expression and Abundance Analysis ===")

# Make sure group is categorical for ordering
if not pd.api.types.is_categorical_dtype(model.sample_info['group']):
    model.sample_info['group'] = model.sample_info['group'].astype('category')

# Define the specific comparisons we want
comparisons = [
    ("Female_Control", "Female_Disease"),
    ("Male_Control", "Male_Disease")
]

# Create variables to store results for the final plot
de_results = {}
da_results = {}

# Process DE for each comparison
for ref, target in comparisons:
    print(f"\nAnalyzing {target} vs {ref}")
    
    # Reorder categories to make 'ref' the reference group
    ordered_groups = [ref] + [g for g in adata.obs['group'].unique() if g != ref]
    model.sample_info['group'] = model.sample_info['group'].cat.reorder_categories(ordered_groups)
    
    # Run DE analysis
    print(f"Computing differential expression: {target} vs {ref}")
    de_res = model.differential_expression(
        sample_cov_keys=["group"],
        store_lfc=True,
    )
    
    # Save the full DE results
    de_res.to_netcdf(f"results/de_{target}_vs_{ref}.nc")
    
    # Extract the covariate name for this comparison
    cov_name = None
    for cov in de_res["effect_size"].coords["covariate"].values:
        if target in cov:
            cov_name = cov
            break
    
    if cov_name:
        # Store effect size in adata
        effect_col = f"{target}_vs_{ref}_DE"
        adata.obs[effect_col] = de_res["effect_size"].sel(covariate=cov_name).values
        
        # Save for final plot
        de_results[(ref, target)] = effect_col
        
        # Save top DE genes
        if "lfc" in de_res:
            try:
                mean_lfc = np.abs(de_res["lfc"].sel(covariate=cov_name).mean(dim="cell_name"))
                top_genes = mean_lfc.to_pandas().nlargest(10).index.tolist()
                pd.DataFrame({f"{target} vs {ref}": top_genes}).to_csv(
                    f"results/top_genes_{target}_vs_{ref}.csv", index=False
                )
            except Exception as e:
                print(f"Error getting top genes: {str(e)}")

# Compute differential abundance
print("\nComputing differential abundance")
da_res = model.differential_abundance(sample_cov_keys=["group"])
da_res.to_netcdf("results/da_results.nc")

# Process DA for each comparison
logp = da_res["group_log_probs"]

for ref, target in comparisons:
    # Calculate ratio
    print(f"Computing differential abundance: {target} vs {ref}")
    ratio = logp.sel({"group": target}) - logp.sel({"group": ref})
    
    # Store in adata
    ratio_col = f"{target}_vs_{ref}_DA"
    adata.obs[ratio_col] = ratio.values
    
    # Save for final plot
    da_results[(ref, target)] = ratio_col

# Create the combined plot (Group UMAP + 2x2 grid of comparisons)
print("\nCreating comparison grid plot")
fig = plt.figure(figsize=(17, 12))

# Define grid layout
gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1, 1])

# Left panel - Group UMAP
ax_groups = fig.add_subplot(gs[:, 0])
sc.pl.umap(adata, color="display_group", ax=ax_groups, show=False, title="Cell Groups", legend_loc='upper left')

# Right panel - 2x2 grid of DE and DA
axes = [
    fig.add_subplot(gs[0, 1]),  # Male DE
    fig.add_subplot(gs[0, 2]),  # Female DE
    fig.add_subplot(gs[1, 1]),  # Male DA
    fig.add_subplot(gs[1, 2]),  # Female DA
]

# Plot DE and DA
comparison_titles = [
    "Male disease vs control (DE)",
    "Female disease vs control (DE)",
    "Male disease vs control (DA)",
    "Female disease vs control (DA)"
]

i = 0
for c_type, results in [("DE", de_results), ("DA", da_results)]:
    for (ref, target) in comparisons:
        if (ref, target) in results:
            col = results[(ref, target)]
            
            # Get proper axis based on comparison
            if "Male" in target:
                ax_idx = 0 if c_type == "DE" else 2
            else:
                ax_idx = 1 if c_type == "DE" else 3
            
            # Set proper cmap and limits
            if c_type == "DE":
                cmap = "viridis"
                vmin = None
                vmax = np.quantile(adata.obs[col].values, 0.95)
            else:
                cmap = "coolwarm"
                vmin = -1
                vmax = 1
            
            # Plot
            sc.pl.umap(
                adata, 
                color=col,
                ax=axes[ax_idx],
                show=False,
                title=comparison_titles[ax_idx],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax
            )
            i += 1

# Adjust layout and save
plt.tight_layout()
plt.savefig("figs/comparison_grid.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figs/comparison_grid.png", dpi=300, bbox_inches="tight")
plt.close()

print("Analysis complete!")
