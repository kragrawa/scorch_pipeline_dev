# DoubletDetection

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

# source Junchen's helper functions
source("/data/kriti/scorch_pipeline_dev/r_helper_functions/seurat_wrapper_funs.R")

# parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
loaded_data <- args[1]
output_file <- args[2]

data_S_list <- readRDS(loaded_data)

data_S_list
# Process Sample for doublet detection
data_S_list <- ProcessSingleSample(data_S_list)

# Doublet_Detection using DoubletFinder --> this step can take some time
data_S_list <- scDoubletDetection(data_S_list)

metadata <- data.frame()
embeddings <- data.frame()
for(i in names(data_S_list)){
  metadata <- rbind(metadata, data_S_list[[i]][[]])
  embeddings <- rbind(embeddings, Embeddings(object = data_S_list[[i]][["pca"]]))
}

write.csv(metadata, "metadata_doublet_dectection.csv")
write.csv(embeddings, "embeddings_doublet_dectection_pca.csv")

# save individual objects with doublets
saveRDS(data_S_list, output_file)

sessionInfo()