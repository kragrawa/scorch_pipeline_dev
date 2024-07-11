# Label Transfer
library(Seurat)
library(SeuratObject)
library(scales)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(viridis)
library(grid)

# source Junchen's helper functions
source("/data/kriti/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
reference <- args[2]
ref_map_col <- args[3]
output_file <- args[4]

integrated_data <- readRDS(loaded_data)
reference_data <- readRDS(reference)

anchors <- FindTransferAnchors(reference_data,
    integrated_data, reference.reduction = "integrated.rpca", dims = 1:30)

predictions <- TransferData(
  anchorset = anchors,
  refdata = reference_data[[]][, ref_map_col]
)

integrated_data <- AddMetaData(object = integrated_data, metadata = predictions)

p1 <- DimPlot(integrated_data, reduction = "umap.rpca", group.by = c("predicted.id"))

ggsave(filename = "umap_integrated_rpca_labeled.png", plot = p1, width = 10, height = 10, dpi = 700)

saveRDS(integrated_data, output_file)

sessionInfo()