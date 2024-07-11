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
clean_file <- args[3]
metadata_file <- args[4]
pca_file <- args[5]

sample <- readRDS(loaded_data)

# process the single sample for doublet detection
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
sample <- ScaleData(sample)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))

# Doublet_Detection using scDblFinder --> this step can take some time
# convert Seurat object to v3 and then to SingleCellExperiment object
sample[["RNA"]] <- as(object = sample[["RNA"]], Class = "Assay")
data_sce <- as.SingleCellExperiment(sample)

# run scDblFinder
data_sce <- scDblFinder(data_sce)

# save the results to seurat object
sample$scDblFinder.class <- data_sce$scDblFinder.class
sample$scDblFinder.score <- data_sce$scDblFinder.score

# save the metadata and embeddings
metadata <- data.frame()
embeddings <- data.frame()

metadata <- rbind(metadata, sample[[]])
embeddings <- rbind(embeddings, Embeddings(object = sample[["pca"]]))

write.csv(metadata, metadata_file)
write.csv(embeddings, pca_file)

# save individual objects with doublets
saveRDS(sample, output_file)

# remove the doublets and save the cleaned object
Idents(sample) <- "scDblFinder.class"
sample <- subset(sample, idents = "singlet")
saveRDS(sample, clean_file)

sessionInfo()