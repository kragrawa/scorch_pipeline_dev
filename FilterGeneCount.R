#Allen Processing

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

individual_data <- readRDS(loaded_data)

metadata <- data.frame()

for(i in names(individual_data)){
  individual_data[[i]] <- RenameCells(individual_data[[i]], add.cell.id = i)

  #add into the metadata
  individual_data[[i]]$nFeatureRNA_passed <- ifelse(individual_data[[i]]$nFeature_RNA > 500 & individual_data[[i]]$nFeature_RNA <= 7500, "yes", "no")

  #add the metadata
  metadata <- rbind(metadata, individual_data[[i]][[]])

  #subset the object
  individual_data[[i]] <- subset(
    x = individual_data[[i]],
    subset = nFeature_RNA > 500 & nFeature_RNA <= 7500
  )
}

#save the objects
write.csv(metadata, "metadata_nFeatureRNA.csv")
saveRDS(individual_data, output_file)

sessionInfo()