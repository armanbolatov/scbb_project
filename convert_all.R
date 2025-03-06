library(qs)
library(Seurat)
library(SeuratDisk)

input_dir <- "./downloaded"   
output_dir <- "./ann_data" 

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

qsave_files <- list.files(input_dir, pattern = "\\.qsave$", full.names = TRUE)

convert_qsave_to_h5ad <- function(qsave_path) {
  file_name <- basename(qsave_path)  
  file_base <- sub("\\.qsave$", "", file_name)  

  seurat_h5_path <- file.path(output_dir, paste0(file_base, ".h5Seurat"))
  h5ad_path <- file.path(output_dir, paste0(file_base, ".h5ad"))

  data_obj <- qs::qread(qsave_path)

  SaveH5Seurat(data_obj, filename = seurat_h5_path, overwrite = TRUE)

  Convert(seurat_h5_path, dest = "h5ad", overwrite = TRUE)

  final_h5ad_path <- file.path(output_dir, paste0(file_base, ".h5ad"))
  file.rename(h5ad_path, final_h5ad_path)

  unlink(seurat_h5_path)

  cat(paste("Converted:", qsave_path, "->", final_h5ad_path, "and deleted", seurat_h5_path, "\n"))
}

lapply(qsave_files, convert_qsave_to_h5ad)

cat("All conversions completed!\n")
