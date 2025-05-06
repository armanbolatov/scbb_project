library(Seurat)
library(Matrix)
library(tools)

base_dir    <- "/l/users/darya.taratynova/scbb_project/GSE160936/filtered_matrices"
all_dirs    <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

# keep only those with an rds/barcodes.tsv.gz inside
sample_dirs <- all_dirs[
  file.exists(file.path(all_dirs, "rds", "barcodes.tsv.gz"))
]

marker_dir   <- "/l/users/darya.taratynova/scbb_project/_process_again/cell_types"
marker_files <- list.files(marker_dir, "\\.csv$", full.names=TRUE)
marker_list  <- setNames(
  lapply(marker_files, function(f) read.csv(f, stringsAsFactors=FALSE)$marker),
  file_path_sans_ext(basename(marker_files))
)

if (!dir.exists("figures")) dir.create("figures")

for (sd in sample_dirs) {
  sample_name <- basename(sd)
  message("→ Processing sample ", sample_name)
  
  mat_dir <- file.path(sd, "rds")
  counts <- Read10X(data.dir = mat_dir, gene.column = 2)
  
  seu <- CreateSeuratObject(counts,
                            project = sample_name,
                            min.cells = 3,
                            min.features = 200)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, vars.to.regress = c("percent.mt","nCount_RNA"))
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.8)
  seu <- RunUMAP(seu, dims = 1:20)

  for (ct in names(marker_list)) {
    seu <- AddModuleScore(
      object   = seu,
      features = list(marker_list[[ct]]),
      name     = paste0(ct, "_score")
    )
  }
  score_cols <- paste0(names(marker_list), "_score1")
  winner     <- apply(seu@meta.data[, score_cols, drop=FALSE], 1, which.max)
  seu$cell_type_marker <- names(marker_list)[winner]
  if (!dir.exists("results")) dir.create("results")

  # 1) Save just the cell_type_marker per barcode
  celltypes <- data.frame(
    barcode    = colnames(seu),
    cell_type  = seu$cell_type_marker,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  write.csv(
    celltypes,
    file = file.path("results", paste0(sample_name, "_cell_types.csv")),
    quote = FALSE,
    row.names = FALSE
  )
  
  # 2) (Optionally) save the full metadata, including your one-hot scores
  write.csv(
    seu@meta.data,
    file = file.path("results", paste0(sample_name, "_metadata_full.csv")),
    quote = FALSE
  )

  message("✅ Done ", sample_name, "\n")
}