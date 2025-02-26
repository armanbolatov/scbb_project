if (!requireNamespace("qs", quietly = TRUE)) install.packages("qs")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)
library(qs)

options(timeout = 600)

# create directories
dirs <- c("raw_data", "processed_data")
sapply(dirs, function(d) if (!dir.exists(d)) dir.create(d))

# project ids
project_ids <- sprintf("AD001%02d", 1:10)

# process each dataset
for (proj_id in project_ids) {
  qsave_file <- file.path("raw_data", paste0(proj_id, ".qsave"))
  h5seurat_file <- file.path("processed_data", paste0(proj_id, ".h5Seurat"))
  h5ad_file <- file.path("processed_data", paste0(proj_id, ".h5ad"))
  
  if (!file.exists(qsave_file)) {
    cat(sprintf("File %s does not exist. Skipping.\n", qsave_file))
    next
  }

  tryCatch({
    cat(sprintf("Processing %s...\n", proj_id))
    data_obj <- qs::qread(qsave_file)
    
    # convert to h5ad format
    SaveH5Seurat(data_obj, filename = h5seurat_file, overwrite = TRUE)
    Convert(h5seurat_file, dest = "h5ad", overwrite = TRUE)
    
    # remove intermediate files to save space
    file.remove(h5seurat_file)
    
    cat(sprintf("Successfully processed %s\n", proj_id))
  }, error = function(e) {
    cat(sprintf("Error processing %s: %s\n", proj_id, e$message))
  })
}