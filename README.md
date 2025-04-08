# Y-Scorch Pipeline Development

We used nextflow for initial sample-level and region level processing. This pipeline has not been written to be used by others without modification at the moment. Some of the file paths are hard coded and would need to altered for use from others, but all files to reproduce the analysis can be found here. 

Description of support files:
1. seurat_v5.dockerfile -> Docker file with packages needed to run the pipeline
2. nextflow.config -> some processes are data intensive and require more memory and cpus

Pipeline files:
scorch_R_test1.nf -> Nexflow file to run the pipeline (equivalent to main.nf)
The data for the paper using testing as the workflow, so please focus on the following formats. Other versions will be removed and cleaned before publication.
Testing Workflow:
  1. LoadDataFromSampleSheet.R - Loads data from a samplesheet (example can be found in data/test_sample_sheet.csv)
  2. FilterMitoAndSexGenesSingleSample.R - Filters out mitochondrial and x and y chromosome genes
  3. nFeatureRNA_SingleSample.R - Filter cells with 500 < nFeatureRNA < 7500 
  4. DoubletDetectionSingleSample.R - Performs doublet detection using scDblFinder
  5. LabelTransferSingleSample.R - Performs label transfer and consensus annotations using Ma et al and BICCN
  6. LabelTransferVST.R - Performs label transfer and consensus annotations using Ma et al, BICCN, and non-human primate reference for VST

Other Workflows:
Originally we did not use a sample sheet so there are some alterative files that run essentially the same pipeline, but actually takes in a data directory as the input. Those files are not relevant for the paper as the output from those workflows was not used for the paper. 
