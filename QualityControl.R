#QualityControl - R script

library(Seurat)
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
output_file <- args[2]
mito_percent <- as.numeric(args[3])

individual_data <- readRDS(loaded_data)
data_S_list <- FilterCells(individual_data, ngene_lth = 500,
                           ngene_hth = 7500, mt_hth = mito_percent)
saveRDS(data_S_list, output_file)

sessionInfo()