library(Seurat)
library(hdf5r)
library(Matrix)
library(tools)    # for file_path_sans_ext()

# 0) Paths to your snRNA‐seq files
h5f    <- "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
meta_f <- "/l/users/darya.taratynova/scbb_project/GSE174367/sn_data/GSE174367_snRNA-seq_cell_meta.csv"

# 1) Load the counts matrix from the HDF5
counts <- Read10X_h5(h5f)
seu    <- CreateSeuratObject(counts, project = "GSE174367_snRNA",
                             min.cells = 3, min.features = 200)

# 2) Bring in the author’s metadata
meta <- read.csv(meta_f, stringsAsFactors = FALSE)
# make sure rownames match barcodes in the matrix
rownames(meta) <- meta$Barcode
seu <- AddMetaData(seu, meta)

# 3) (Optional) QC filter by nFeature and percent.mt
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 5)

# 4) Pre‐process (normalization, PCA, UMAP—so you can visualize)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, vars.to.regress = c("percent.mt","nCount_RNA"))
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.8)
seu <- RunUMAP(seu, dims = 1:20)

# 5) Load your marker lists (same cell_types folder as before)
marker_dir   <- "/l/users/darya.taratynova/scbb_project/_process_again/cell_types"
marker_files <- list.files(marker_dir, "\\.csv$", full.names = TRUE)
marker_list  <- setNames(
  lapply(marker_files, \(f) read.csv(f, stringsAsFactors=FALSE)$marker),
  file_path_sans_ext(basename(marker_files))
)

# 6) Prepare & filter marker lists against your data’s genes
# uppercase everything to avoid mismatches
rownames(seu)   <- toupper(rownames(seu))
marker_list     <- lapply(marker_list, toupper)
marker_list     <- lapply(marker_list, intersect, rownames(seu))

# 7) AddModuleScore for each marker set (≥3 genes), else zero
for (ct in names(marker_list)) {
  genes_present <- marker_list[[ct]]
  score_name    <- paste0(ct, "_score")
  if (length(genes_present) >= 3) {
    seu <- AddModuleScore(seu,
                          features = list(genes_present),
                          name     = score_name)
  } else {
    seu[[paste0(score_name, "1")]] <- 0
  }
}

# 8) Winner‐takes‐all
score_cols <- paste0(names(marker_list), "_score1")
winner     <- apply(seu@meta.data[, score_cols, drop = FALSE], 1, which.max)
seu$cell_type_marker <- names(marker_list)[winner]

# 9) Save out your cell‐type calls
if (!dir.exists("results")) dir.create("results")
celltypes <- data.frame(
  barcode   = colnames(seu),
  cell_type = seu$cell_type_marker,
  stringsAsFactors = FALSE
)
write.csv(celltypes,
          file      = file.path("results", "GSE174367_snRNA_cell_types.csv"),
          quote     = FALSE,
          row.names = FALSE)

dev.off()
