#FilterMitoAndSexGenes

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

individual_data <- readRDS(loaded_data)
individual_data
chrX_genes <- read.table(chrX_txt, header = FALSE)
chrY_genes <- read.table(chrY_txt, header = FALSE)

mt_genes <- grep("^MT-", rownames(individual_data[[1]]), value = TRUE)

genes_to_remove <- c(chrX_genes$V2, chrY_genes$V2, mt_genes)

filteredSamples <- list()

options(Seurat.object.assay.version = "v5")

for(i in names(individual_data)){
  data <- individual_data[[i]][["RNA"]]$counts
  sample_metadata <- individual_data[[i]][[]]
  filteredData <- data[!rownames(data) %in% genes_to_remove, ]

  filteredSamples[[i]] <- CreateSeuratObject(counts = filteredData, project = i)
  filteredSamples[[i]] <- AddMetaData(filteredSamples[[i]], sample_metadata[Cells(filteredSamples[[i]]), ])
  # save the counts matrix for each file
  filtered_counts <- filteredSamples[[i]][["RNA"]]$counts
  writeMM(filtered_counts, paste0(i, "_no_mito_sex_counts.mtx"))
  gzip(paste0(i, "_no_mito_sex_counts.mtx"))
}

saveRDS(filteredSamples, output_file)

sessionInfo()