# CreateMetaData
# Function to parse sample names and create conditions
library(dplyr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_file <- args[2]

counter_env <- new.env()
counter_env$CTR <- 0
counter_env$HIV <- 0
counter_env$HIVOUD <- 0
counter_env$OUD <- 0

parse_conditions_with_counter <- function(name, env) {
  # Determine the condition based on the pattern in the name
  if (grepl("_CTR_", name)) {
    condition <- "CTR"
  } else if (grepl("_HIV_", name)) {
    condition <- "HIV"
  } else if(grepl("_HIVOUD_", name)){
    condition <- "HIVOUD"
  } else if(grepl("_OUD_", name)){
    condition <- "OUD"
  } else {
    condition <- "Unknown"
  }

  # Increment the counter for the identified condition
  if (condition != "Unknown") {
    env[[condition]] <- env[[condition]] + 1
    # Combine the condition with the counter value to create the identifier
    condition_identifier <- paste(condition, env[[condition]], sep = "_")
  } else {
    condition_identifier <- condition
  }
  return(condition_identifier)
}

# Initialize the sample_metadata data frame
sample_metadata <- data.frame("sample_names" = list.dirs(data_dir, recursive = F, full.names = F))

# Apply the function to each sample name using a loop
condition_identifiers <- character(length(sample_metadata$sample_names))
for (i in seq_along(sample_metadata$sample_names)) {
  condition_identifiers[i] <- parse_conditions_with_counter(sample_metadata$sample_names[i], counter_env)
}

sample_metadata$sample_ids <- condition_identifiers
sample_metadata$condition <- sapply(sample_metadata$sample_ids, function(id) strsplit(id, split = "_")[[1]][1])
saveRDS(sample_metadata, output_file)
