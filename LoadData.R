#LoadData - R script

library(Seurat) 
library(scales)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(viridis)
library(grid)
library(dplyr)

#parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

#assuming the first argument is a file path and the second is an output path
sample_data_dir <- args[1]
metadata <- args[2]
output_file <- args[3]
print(sample_data_dir)
print(output_file)

#source Junchen's helper functions
source("/data/kriti/pipeline_dev/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

sample_metadata <- readRDS(metadata)

sample_names <- list.dirs(sample_data_dir, full.names = F, recursive = F)
sample_names
data_S_list_v0 <- Load10xData(sample_data_dir, sample_names, with.rna.only = T, sub_rna_dir = "outs/filtered_feature_bc_matrix")

#add in the metadata
for(sample in sample_names){
  metadata_row <- sample_metadata[sample_metadata$sample_names %in% sample, ]
  metadata_list <- as.list(metadata_row[-which(names(metadata_row) == "sample_names")])
  data_S_list_v0[[sample]] <- AddMetaData(data_S_list_v0[[sample]], metadata = metadata_list)
}

saveRDS(data_S_list_v0, output_file)