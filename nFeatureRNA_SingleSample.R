# nFeatureRNA Filter

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
output_file <- args[2]
metadata_filename <- args[3]

sample <- readRDS(loaded_data)

metadata <- data.frame()

sample$nFeatureRNA_passed <- ifelse(sample$nFeature_RNA > 500 & sample$nFeature_RNA <= 7500, "yes", "no")

metadata <- sample[[]]

sample <- subset(
  x = sample,
  subset = nFeature_RNA > 500 & nFeature_RNA <= 7500
)

# save the objects
write.csv(metadata, metadata_filename)
saveRDS(sample, output_file)

sessionInfo()