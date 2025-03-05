import rpy2.robjects as ro
from rpy2.robjects.packages import importr

r_code = """
library(qs)
library(Seurat)
library(SeuratDisk)

# Load Seurat object
data_obj <- qs::qread("./unprocessed_data/AD00202.qsave")

# Save as H5Seurat format
SaveH5Seurat(data_obj, filename = "seurat_data.h5Seurat", overwrite = TRUE)

# Convert to H5AD (Scanpy format)
Convert("seurat_data.h5Seurat", dest = "h5ad", overwrite = TRUE)
"""

ro.r(r_code)

print("Conversion completed! seurat_data.h5ad is ready.")