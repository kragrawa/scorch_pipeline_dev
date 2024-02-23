# MergeData

library(Seurat)
library(dplyr)
library(ggplot2)

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

# Remove the Doublets
data_S <- subset(data_S, cells = colnames(data_S)[data_S@meta.data$scDblFinder.class == "singlet"]) # nolint: line_length_linter.

# Split the data into RNA for integration in the future
data_S[["RNA"]] <- split(data_S[["RNA"]], f = data_S$sample_ids)

# 4. Data normalization and dimensionality reduction
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 30, verbose = F)

# Unintegrated clustering
data_S <- FindNeighbors(data_S, dims = 1:30, reduction = "pca")
data_S <- FindClusters(data_S, resolution = 2, cluster.name = "unintegrated_clusters")
data_S <- RunUMAP(data_S, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

p1 <- DimPlot(data_S, reduction = "umap.unintegrated", group.by = c("condition"))
p2 <- DimPlot(data_S, reduction = "umap.unintegrated", group.by = c("sample_ids"))
p3 <- DimPlot(data_S, reduction = "umap.unintegrated", group.by = c("seurat_clusters"))

ggsave(filename = "umap_unintegrated_condition.png", plot = p1, width = 10, height = 10, dpi = 700)
ggsave(filename = "umap_unintegrated_sample_ids.png", plot = p2, width = 10, height = 10, dpi = 700)
ggsave(filename = "umap_unintegrated_seurat_clusters.png", plot = p3, width = 10, height = 10, dpi = 700)

# save individual objects with doublets
saveRDS(data_S, output_file)

sessionInfo()