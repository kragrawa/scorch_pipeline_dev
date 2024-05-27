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
source("/data/kriti/pipeline_dev/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
biccn_ref <- args[2]
biccn_map_col <- args[3]
ma_ref <- args[4]
ma_map_col <- args[5]
sample_name <- args[6]
output_file <- args[7]

sample <- readRDS(loaded_data)
biccn_reference <- readRDS(biccn_ref)
ma_reference <- readRDS(ma_ref)

#both the ma et al reference are normalized using SCTransform
DefaultAssay(sample) <- "RNA"
sample <- SCTransform(sample)
sample <- RunPCA(sample)

#get the reference anchors for biccn
transfer_anchors_BICCN <- FindTransferAnchors(reference = biccn_reference,
                          query = sample, dims = 1:30, reduction = "pcaproject")
predictions_BICCN <- TransferData(
  anchorset = transfer_anchors_BICCN,
  refdata = biccn_reference[[]][, biccn_map_col], dims = 1:30
)

#add the biccn prefix to the column names
colnames(predictions_BICCN) <- paste("biccn", colnames(predictions_BICCN), sep = "_")
write.csv(predictions_BICCN, paste0("biccn_pred_", sample_name, ".csv"))

#get the reference anchors for ma et al
transfer_anchors_MA <- FindTransferAnchors(reference = ma_reference,
                          query = sample, dims = 1:30, reduction = "pcaproject")
predictions_MA <- TransferData(
  anchorset = transfer_anchors_MA,
  refdata = ma_reference[[]][, ma_map_col], dims = 1:30
)
colnames(predictions_MA) <- paste("ma", colnames(predictions_MA), sep = "_")
write.csv(predictions_MA, paste0("ma_pred_", sample_name, ".csv"))

#add the metadata to the object
sample <- AddMetaData(object = sample, metadata = predictions_BICCN)
sample <- AddMetaData(object = sample, metadata = predictions_MA)

#plot the predicted labels
sample <- RunUMAP(sample, dims = 1:30)
p1 <- DimPlot(sample, reduction = "umap", group.by = c("biccn_predicted.id"))
ggsave(filename = paste0("umap_labeled_biccn_", sample_name, ".png"), plot = p1, width = 10, height = 10, dpi = 700)

p2 <- DimPlot(sample, reduction = "umap", group.by = c("ma_predicted.id"))
ggsave(filename = paste0("umap_labeled_ma_", sample_name, ".png"), plot = p2, width = 10, height = 10, dpi = 700)

saveRDS(sample, output_file)

sessionInfo()