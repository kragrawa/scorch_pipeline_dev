# Integrate Data

library(Seurat)
library(dplyr)
library(ggplot2)

# source Junchen's helper functions
source("/data/kriti/pipeline_dev/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
output_file <- args[2]

options(future.globals.maxSize = 2000 * 1024^2)

merged_data <- readRDS(loaded_data)
#merged_data[["RNA"]] <- as(object = merged_data[["RNA"]], Class = "Assay5")

integrated_data <- IntegrateLayers(
  object = merged_data, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "LogNormalize")

integrated_data <- FindNeighbors(integrated_data, reduction = "integrated.rpca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 2, cluster.name = "rpca_clusters")

integrated_data <- RunUMAP(integrated_data, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

p1 <- DimPlot(integrated_data, reduction = "umap.rpca", group.by = c("condition"))
p2 <- DimPlot(integrated_data, reduction = "umap.rpca", group.by = c("sample_ids"))
p3 <- DimPlot(integrated_data, reduction = "umap.rpca", group.by = c("rpca_clusters"))

ggsave(filename = "umap_integrated_rpca_condition.png", plot = p1, width = 10, height = 10, dpi = 700)
ggsave(filename = "umap_integrated_rpca_sample_ids.png", plot = p2, width = 10, height = 10, dpi = 700)
ggsave(filename = "umap_integrated_rpca_seurat_clusters.png", plot = p3, width = 10, height = 10, dpi = 700)

#write the embeddings
write.csv(Embeddings(integrated_data[["integrated.rpca"]]), "embeddings_umap_integrated_rpca.csv")

# save individual objects with doublets
saveRDS(integrated_data, output_file)

sessionInfo()