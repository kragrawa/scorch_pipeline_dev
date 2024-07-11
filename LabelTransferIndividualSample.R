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

data_S_list <- readRDS(loaded_data)
reference_data <- readRDS(reference)

for (i in names(data_S_list)) {
  #Run UMAP for each sample for visualization, PCA was run in the doublet detection step (using 50 dims for UMAP)
  data_S_list[[i]] <- RunUMAP(data_S_list[[i]], dims = 1:50)

  #get the reference anchors
  anchors <- FindTransferAnchors(reference_data,
    data_S_list[[i]], reference.reduction = "integrated.rpca", dims = 1:50)
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_data[[]][, ref_map_col]
  )

  #add the metadata to the object
  data_S_list[[i]] <- AddMetaData(object = data_S_list[[i]], metadata = predictions)

  #plot the predicted labels
  p1 <- DimPlot(data_S_list[[i]], reduction = "umap", group.by = c("predicted.id"))
  ggsave(filename = paste0(i,"_umap_labeled.png"), plot = p1, width = 10, height = 10, dpi = 700)
}

saveRDS(data_S_list, output_file)

sessionInfo()