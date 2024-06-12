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

#consensus annotations
high_level_map <- list(
  "MSNs" = "MSNs",
  "Oligos" = "Glial",
  "Microglia" = "Glial",
  "Oligos_Pre" = "Glial",
  "Astrocytes" = "Glial",
  "Mural/Fibroblast" = "Nonneuronal",
  "Endothelial" = "Nonneuronal",
  "Interneurons" = "InterN",
  "L2/3 IT" = "ExN",
  "L6 IT" = "ExN",
  "L4 IT" = "ExN",
  "L5 IT" = "ExN",
  "L6 IT Car3" = "ExN",
  "L5 ET" = "ExN",
  "L6 CT" = "ExN",
  "L5/6 NP" = "ExN",
  "L6b" = "ExN",
  "Lamp5" = "InN",
  "Vip" = "InN",
  "Pax6" = "InN",
  "Sncg" = "InN",
  "Lamp5 Lhx6" = "InN",
  "Chandelier" = "InN",
  "Pvalb" = "InN",
  "Sst" = "InN",
  "Sst Chodl" = "InN",
  "Astro" = "Glial",
  "Endo" = "Nonneuronal",
  "VLMC" = "Nonneuronal",
  "Oligo" = "Glial",
  "Micro/PVM" = "Glial",
  "OPC" = "Glial",
  "L3-5 IT-1" = "ExN",
  "L2-3 IT" = "ExN",
  "SST HGF" = "InN",
  "VIP" = "InN",
  "L3-5 IT-3" = "ExN",
  "L6 IT-1" = "ExN",
  "Micro" = "Glial",
  "L5-6 NP" = "ExN",
  "L6 IT-2" = "ExN",
  "PC" = "Nonneuronal",
  "L3-5 IT-2" = "ExN",
  "LAMP5 LHX6" = "InN",
  "Immune" = "Nonneuronal",
  "PVALB" = "InN",
  "ADARB2 KCNG1" = "InN",
  "SST" = "InN",
  "L6B" = "ExN",
  "LAMP5 RELN" = "InN",
  "SMC" = "Nonneuronal",
  "SST NPY" = "InN",
  "PVALB ChC" = "InN",
  "RB" = "Nonneuronal"
)

high_level_vector <- unlist(high_level_map)

metadata <- sample@meta.data

# add high_level_mapping and consensus_annotation to metadata
metadata <- metadata %>%
  mutate(
    high_level_ct_ma = high_level_vector[ma_predicted.id],
    high_level_ct_biccn = high_level_vector[biccn_predicted.id],
    high_level_ct = coalesce(high_level_ct_biccn, high_level_ct_ma),
    consensus_annotation = if_else(high_level_ct_biccn %in% c("ExN", "InN"), biccn_predicted.id, ma_predicted.id) # nolint
  )

#add metadata back into the object
sample@meta.data <- metadata

#plot the predicted labels
sample <- RunUMAP(sample, dims = 1:30)
p1 <- DimPlot(sample, reduction = "umap", group.by = c("biccn_predicted.id"))
ggsave(filename = paste0("umap_labeled_biccn_", sample_name, ".png"), plot = p1, width = 10, height = 10, dpi = 700)

p2 <- DimPlot(sample, reduction = "umap", group.by = c("ma_predicted.id"))
ggsave(filename = paste0("umap_labeled_ma_", sample_name, ".png"), plot = p2, width = 10, height = 10, dpi = 700)

p3 <- DimPlot(sample, reduction = "umap", group.by = c("consensus_annotation"), label = TRUE)
ggsave(filename = paste0("umap_consensus_", sample_name, ".png"), plot = p3, width = 10, height = 10, dpi = 700)

p4 <- DimPlot(sample, reduction = "umap", group.by = c("high_level_ct"), label = TRUE)
ggsave(filename = paste0("umap_high_level_ct_", sample_name, ".png"), plot = p4, width = 10, height = 10, dpi = 700)

write.csv(metadata, paste0("labeled_metadata_", sample_name, ".csv"))
final_counts <- sample[["RNA"]]$counts
saveRDS(final_counts, paste0("final_counts_", sample_name, ".rds"))
saveRDS(sample, output_file)

sessionInfo()