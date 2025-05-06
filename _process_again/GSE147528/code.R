# #===============================================================================
# # Cell-type annotation pipeline in R (Seurat + SingleR)
# #===============================================================================

# # 0) Install missing packages (run once)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")

# # CRAN packages
# install.packages(c("Seurat", "ggplot2", "readr", "dplyr"), dependencies=TRUE)

# # Bioconductor packages for SingleR
# BiocManager::install(c("SingleCellExperiment", "scater", "SingleR", "celldex"), ask=FALSE)

# #===============================================================================
# # 1) Load libraries
# #===============================================================================
# library(Seurat)
# library(SingleCellExperiment)
# library(scater)
# library(hdf5r)

# #===============================================================================
# # 2) Load 10x raw counts and create Seurat object
# #===============================================================================
# h5_path <- "/l/users/darya.taratynova/scbb_project/GSM4120422/GSM4432635_SFG2_raw_gene_bc_matrices_h5.h5"
# counts  <- Read10X_h5(filename = h5_path)

# seurat_obj <- CreateSeuratObject(
#   counts = counts,
#   project = "AD_SingleNuclei",
#   min.cells = 3,
#   min.features = 200
# )
# percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# #===============================================================================
# # 3) Basic QC and filtering (optional)
# #===============================================================================
# seurat_obj[["percent.mt"]] <- percent.mt
# seurat_obj <- subset(
#   seurat_obj,
#   subset = nFeature_RNA > 200 &
#            nFeature_RNA < 6000 &
#            percent.mt < 5
# )
# #===============================================================================
# # 4) Preprocessing: Normalize, HVG, Scale, PCA, Neighbors, UMAP, Clustering
# #===============================================================================
# seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
# seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
# seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("percent.mt", "nCount_RNA"))
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30)

# # Neighbors & clustering
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, k.param = 20)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# # UMAP
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# # Save cluster UMAP
# pdf("figures/umap_clusters.pdf", width=6, height=5)
# DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")
# dev.off()

# #===============================================================================
# # 5) Load marker lists and score cells (marker-based annotation)
# #===============================================================================
# marker_dir <- "/l/users/darya.taratynova/scbb_project/_process_again/cell_types"
# marker_files <- list.files(marker_dir, pattern = "\\.csv$", full.names = TRUE)

# marker_list <- list()
# for (f in marker_files) {
#   ct <- tools::file_path_sans_ext(basename(f))
#   df <- read.csv(f, stringsAsFactors = FALSE)
#   marker_list[[ct]] <- df$marker
# }

# # Add module scores
# for (ct in names(marker_list)) {
#   score_name <- paste0(ct, "_score")
#   seurat_obj <- AddModuleScore(
#     object = seurat_obj,
#     features = list(marker_list[[ct]]),
#     name = score_name
#   )
# }

# # Determine top score per cell
# score_cols <- paste0(names(marker_list), "_score1")
# seurat_obj$cell_type_marker <- apply(
#   seurat_obj@meta.data[, score_cols, drop = FALSE],
#   1,
#   function(x) names(marker_list)[which.max(x)]
# )


# # Save marker-based UMAP
# pdf("figures/umap_marker_cell_types.pdf", width=6, height=5)
# DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type_marker")
# dev.off()
# #===============================================================================
# saveRDS(seurat_obj, file = "results/seurat_obj_celltyped.rds")
# write_csv(seurat_obj@meta.data, path = "results/obs_cell_types.csv")

# # Optionally export AnnData for Python
# # library(zellkonverter)
# # writeH5AD(as.SingleCellExperiment(seurat_obj), "results/adata_celltyped.h5ad")


#===============================================================================
# Batch cell‑typing pipeline (Seurat → SingleCellExperiment → AnnData)
#===============================================================================

#===============================================================================
# Batch cell‑typing pipeline (Seurat → AnnData via sceasy)
#===============================================================================

# 0) One‑off installs
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# CRAN packages
install.packages(c("Seurat", "ggplot2", "readr", "dplyr"), dependencies = TRUE)

# Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "scater", "hdf5r"), ask = FALSE)

# sceasy for AnnData export
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("cellgeni/sceasy")

#===============================================================================
# Batch cell‑typing pipeline (Seurat → flat files for Python → AnnData)
#===============================================================================

# 0) one‑off installs
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
install.packages(c("Seurat","ggplot2","readr","dplyr","Matrix"), dependencies=TRUE)
BiocManager::install(c("SingleCellExperiment","scater","hdf5r"), ask=FALSE)

#===============================================================================
# 1) libraries
#===============================================================================
library(Seurat)
library(scater)
library(hdf5r)
library(Matrix)    # for writeMM()
library(tools)     # for file_path_sans_ext()

#===============================================================================
# 2) load marker lists
#===============================================================================
marker_dir   <- "/l/users/darya.taratynova/scbb_project/_process_again/cell_types"
marker_files <- list.files(marker_dir, "\\.csv$", full.names=TRUE)
marker_list  <- lapply(marker_files, function(f) read.csv(f, stringsAsFactors=FALSE)$marker)
names(marker_list) <- file_path_sans_ext(basename(marker_files))

#===============================================================================
# 3) find your 10x `.h5` files
#===============================================================================
h5_files <- list.files(
  "/l/users/darya.taratynova/scbb_project/GSM4120422",
  "\\.h5$", full.names=TRUE
)

#===============================================================================
# 4) loop
#===============================================================================
for (h5 in h5_files) {
  sample_name <- file_path_sans_ext(basename(h5))
  message("Processing ", sample_name, " …")
  
  # a) Read + QC + preprocess
  counts <- Read10X_h5(h5)
  seu    <- CreateSeuratObject(counts, project=sample_name,
                               min.cells=3, min.features=200)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  seu <- subset(seu,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 6000 &
                         percent.mt < 5)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=2000)
  seu <- ScaleData(seu, vars.to.regress=c("percent.mt","nCount_RNA"))
  seu <- RunPCA(seu, features=VariableFeatures(seu), npcs=30)
  seu <- FindNeighbors(seu, dims=1:20, k.param=20)
  seu <- FindClusters(seu, resolution=0.8)
  seu <- RunUMAP(seu, dims=1:20)
  
  # b) Module‑score & winner‑takes‑all
  for (ct in names(marker_list)) {
    seu <- AddModuleScore(
      object   = seu,
      features = list(marker_list[[ct]]),
      name     = paste0(ct, "_score")
    )
  }
  score_cols <- paste0(names(marker_list), "_score1")
  sc_mat     <- seu@meta.data[, score_cols, drop=FALSE]
  winner     <- apply(sc_mat, 1, which.max)
  seu$cell_type_marker <- names(marker_list)[winner]
  
  # c) One‑hot encode
  for (ct in names(marker_list)) {
    seu[[ct]] <- as.integer(seu$cell_type_marker == ct)
  }
  
  # d) Save UMAP
  pdf(file.path("figures", paste0(sample_name, "_umap_marker.pdf")),
      width=6, height=5)
    DimPlot(seu, reduction="umap", group.by="cell_type_marker") +
      ggtitle(sample_name)
  dev.off()
  
  # e) Export flat files for Python assembly
  dir.create("results", showWarnings=FALSE)
  
  # 1) counts matrix (genes × cells) in MatrixMarket
  mat <- GetAssayData(seu, assay="RNA", slot="counts")
  Matrix::writeMM(mat, file = file.path("results", paste0(sample_name, "_counts.mtx")))
  
  # 2) gene list and barcode list
  write.table(
    rownames(mat),
    file      = file.path("results", paste0(sample_name, "_genes.tsv")),
    sep       = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  write.table(
    colnames(mat),
    file      = file.path("results", paste0(sample_name, "_barcodes.tsv")),
    sep       = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  # 3) cell metadata (including cell_type_marker & one‑hots)
  meta <- seu@meta.data
  write.csv(
    meta,
    file      = file.path("results", paste0(sample_name, "_metadata.csv")),
    quote     = FALSE, row.names = TRUE
  )
  
  # 4) UMAP coordinates
  umap <- Embeddings(seu, "umap")
  write.csv(
    umap,
    file      = file.path("results", paste0(sample_name, "_umap.csv")),
    quote     = FALSE, row.names = TRUE
  )
  
  message("  → done, wrote flat files for ", sample_name, "\n")
}
