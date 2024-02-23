# MergeData

# DoubletDetection

library(Seurat)

# source Junchen's helper functions
source("/data/kriti/pipeline_dev/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
output_file <- args[2]

data_S_list <- readRDS(loaded_data)

message("Merging seurat objects")
data_S <- MergeData(data_S_list)
rm(data_S_list)

data_S <- subset(data_S, cells = colnames(data_S)[data_S@meta.data$scDblFinder.class == "singlet"])

# 4. Data normalization and dimensionality reduction
DefaultAssay(data_S) <- "RNA"
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 30, verbose = F)

data_S_list
# Process Sample for doublet detection
data_S_list <- ProcessSingleSample(data_S_list)

# Doublet_Detection using DoubletFinder --> this step can take some time
data_S_list <- scDoubletDetection(data_S_list)

# save individual objects with doublets
saveRDS(data_S_list, output_file)

sessionInfo()