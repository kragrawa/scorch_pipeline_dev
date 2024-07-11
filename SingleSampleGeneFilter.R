#FilterMitoAndSexGenesSingleSample

library(Seurat)
library(SeuratObject)
library(scales)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(viridis)
library(grid)
library(Matrix)
library(R.utils)

# source Junchen's helper functions
source("/data/kriti/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
output_file <- args[2]
chrX_txt <- args[3]
chrY_txt <- args[4]
sample_name_counts <- args[5]

individual_data <- readRDS(loaded_data)
individual_data
chrX_genes <- read.table(chrX_txt, header = FALSE)
chrY_genes <- read.table(chrY_txt, header = FALSE)


mt_genes <- grep("^MT-", rownames(individual_data), value = TRUE)

genes_to_remove <- c(chrX_genes$V2, chrY_genes$V2, mt_genes)

options(Seurat.object.assay.version = "v5")

data <- individual_data[["RNA"]]$counts
sample_metadata <- individual_data[[]]
filteredData <- data[!rownames(data) %in% genes_to_remove, ]

filteredSample <- CreateSeuratObject(counts = filteredData)
filteredSample <- AddMetaData(filteredSample, sample_metadata[Cells(filteredSample), ])
# save the counts matrix for each file
filtered_counts <- filteredSample[["RNA"]]$counts
writeMM(filtered_counts, sample_name_counts)
gzip(sample_name_counts)

saveRDS(filteredSample, output_file)

sessionInfo()