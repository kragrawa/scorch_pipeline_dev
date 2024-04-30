#LoadDataSampleSheet - R script

library(Seurat) 
library(scales)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(viridis)
library(grid)
library(dplyr)
library(Matrix)
library(R.utils)

# Parse the command-line options
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
metadata <- strsplit(args[2], ",")[[1]]
output <- args[3]

print(input)

# Load the Cell Ranger data into a Seurat object
data <- Read10X(input)
seurat_object <- CreateSeuratObject(counts = data$`Gene Expression`)

# Add the metadata to the Seurat object
for (i in seq_along(metadata)) {
    key_value <- strsplit(metadata[i], "=")[[1]]
    seurat_object@meta.data[, key_value[1]] <- key_value[2]
}

# Save the Seurat object to an RDS file
saveRDS(seurat_object, file = output)