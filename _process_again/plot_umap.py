#!/usr/bin/env python3
import os, glob
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# 1) Paths & files (unchanged) …
ann_dirs = [
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE147528/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE160936/ann_data",
    "/l/users/darya.taratynova/scbb_project/_process_again/GSE174367/ann_data",
]
h5ad_files = []
for d in ann_dirs:
    h5ad_files += glob.glob(os.path.join(d, "*.h5ad"))
if not h5ad_files:
    raise RuntimeError("No .h5ad files!")
out_umap = "combined_umap.png"
out_bar  = "cell_type_distribution.png"

# 2) Load & concatenate (unchanged) …
adatas = []
for fn in h5ad_files:
    ad = sc.read_h5ad(fn)
    ad.var_names_make_unique()
    series = os.path.basename(os.path.dirname(os.path.dirname(fn)))
    ad.obs["dataset"] = series
    adatas.append(ad)
combined = sc.concat(adatas, join="outer")

# 3) Preprocess & UMAP (unchanged) …
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="seurat", subset=True)
sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, svd_solver="arpack")
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=20)
sc.tl.umap(combined)

rename_dict = {
    "astocytes":                       "Astrocytes",
    "endothelial":                      "Endothelial",
    "excitatory_neuron":                "Excitatory Neurons",
    "inhibitory_neuron":                "Inhibitory Neurons",
    "microglia":                        "Microglia",
    "nkt":                              "NK/T Cells",
    "oligodendrocyte_precursor_cells":  "OPCs",
    "oligodendrocytes":                 "Oligodendrocytes",
    "pericytes":                        "Pericytes",
}
combined.obs["cell_type"] = (
    combined.obs["cell_type_marker"]
    .map(rename_dict)
    .fillna(combined.obs["cell_type_marker"])
)
combined = combined[combined.obs["cell_type"] != ""].copy()
combined.obs["cell_type"] = combined.obs["cell_type"].astype("category")
combined.obs["cell_type"] = combined.obs["cell_type"].cat.remove_unused_categories()

# define your custom palette again…
custom_colors = {
    "Astrocytes":        "#FF6B35",
    "Endothelial":       "#03A791",
    "Excitatory Neurons":"#16C47F",
    "Inhibitory Neurons":"#C1121F",
    "Microglia":         "#FFEB00",
    "NK/T Cells":        "#780000",
    "OPCs":              "#01BAEF",
    "Oligodendrocytes":  "#A4907C",
    "Pericytes":         "#D67BFF",
}

# — 1) ONLY the cell_type UMAP —
fig, ax = plt.subplots(figsize=(6,6))
sc.pl.umap(
    combined,
    color="cell_type",
    palette=custom_colors,
    ax=ax,
    show=False,
    title="Annotated Cell Types UMAP"
)
# remove axes
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
plt.tight_layout()
plt.savefig("celltype_umap.png", dpi=300)


# — 2) Stacked-bar of proportions per dataset —
# raw counts:
counts = combined.obs.groupby(["dataset","cell_type"]).size()

# convert to proportions within each dataset:
props = (
    counts
    .groupby(level=0)
    .apply(lambda x: x / x.sum())
    .unstack(fill_value=0)
)

fig2, ax2 = plt.subplots(figsize=(8,5))
props.plot(kind="bar", stacked=True, ax=ax2)
ax2.set_ylabel("Proportion of cells")
ax2.set_xlabel("Dataset")
ax2.set_title("Cell-type proportions per dataset")
ax2.legend(title="Cell type", bbox_to_anchor=(1.02,1), loc="upper left")
plt.tight_layout()
plt.savefig("celltype_proportions.png", dpi=300)

fig, (ax_umap, ax_bar) = plt.subplots(
    1, 2,
    figsize=(12, 6),
    gridspec_kw={'width_ratios': [3, 0.5]}
)
# 1) UMAP of cell types
sc.pl.umap(
    combined,
    color="cell_type",
    palette=custom_colors,
    ax=ax_umap,
    show=False
)
# Clean up UMAP axes
ax_umap.set_xticks([]); ax_umap.set_yticks([])
for spine in ax_umap.spines.values():
    spine.set_visible(False)
category_order = list(custom_colors.keys())
props = props.reindex(columns=category_order, fill_value=0)
bar_colors = [custom_colors[ct] for ct in props.columns]

fig, ax = plt.subplots(figsize=(6,6))
sc.pl.umap(
    combined,
    color="cell_type",
    palette=custom_colors,
    ax=ax,
    show=False,
    title="UMAP: cell types"
)
# remove axes
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
plt.tight_layout()
plt.savefig("celltype_umap_only.png", dpi=300)
plt.show()