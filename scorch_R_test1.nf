nextflow.enable.dsl=2

params.output_dir = "/data/kriti/pipeline_dev/test_output/"
params.input_dir = "/banach2/SCORCH/data/raw/10xMultiome-PFC-HIVOUD_OUD-2pairs-12152022/cellranger_arc"
params.chrX_genes = "/data/kriti/pipeline_dev/scorch_pipeline_dev/data/X_chromosome_gene_names.txt"
params.chrY_genes = "/data/kriti/pipeline_dev/scorch_pipeline_dev/data/Y_chromosome_gene_names.txt"

process LoadData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path data_dir
    path sample_metadata

    output:
    path "individual_samples.rds", emit: raw_samples
    path "*.mtx.gz", emit: raw_counts
    path "raw_metadata.csv", emit: raw_metadata

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/LoadData.R "${data_dir}" $sample_metadata individual_samples.rds
    """
}


process FilterMitoAndSexGenes {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}" , mode: 'copy'
    
    input:
    path raw_samples

    output:
    path "sample_no_mito_no_sex.rds", emit: no_mito_no_sex_samples
    path "*.mtx.gz", emit: no_mito_no_sex_counts

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/FilterMitoAndSexGenes.R $raw_samples sample_no_mito_no_sex.rds ${params.chrX_genes} ${params.chrY_genes}
    """
}

process QualityControl {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}" , mode: 'copy'
    
    input:
    path raw_samples

    output:
    path "seperate_mito_gene_filtered.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/QualityControl.R $raw_samples seperate_mito_gene_filtered.rds 2
    """
}

process QualityControlGenes {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path raw_samples

    output:
    path "allen_seperate_mito_gene_filtered.rds", emit: filtered_samples
    path "metadata_nFeatureRNA.csv", emit: metadata_nFeatureRNA

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/FilterGeneCount.R $raw_samples allen_seperate_mito_gene_filtered.rds
    """
}


process DoubletDetection {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}" , mode: 'copy'
    
    input:
    path filtered_samples

    output:
    path "seperate_doublets_filtered.rds", emit: doublet_samples
    path "metadata_doublet_dectection.csv", emit: metadata_doublet_detection
    path "embeddings_doublet_dectection_pca.csv", emit: embeddings_doublet_detection_pca

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/DoubleDetection.R $filtered_samples seperate_doublets_filtered.rds
    """
}


process CreateMetaData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}" , mode: 'copy'

    input:
    path data_dir

    output:
    path "sample_metadata.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/CreateMetaData.R $data_dir sample_metadata.rds
    """
}

process MergeData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}" , mode: 'copy'
    
    input:
    path clean_samples

    output:
    path "merged_data.rds", emit: merged_data
    path "*.png"
    path "embeddings_umap_unintegrated.csv", emit: embeddings_umap_unintegrated

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/MergeData.R $clean_samples merged_data.rds
    """
}

process IntegrateData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path merged_data
   
    output:
    path "integrated_data.rds", emit: integrated_data
    path "*.png"
    path "embeddings_umap_integrated_rpca.csv", emit: embeddings_umap_integrated_rpca
 
    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/IntegrateData.R $merged_data integrated_data.rds
    """
}


process LabelTransfer{
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path integrated_data
    path reference_data

    output:
    path "labeled_data.rds", emit: labeled_data
    path "*.png"
 
    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/LabelTransfer.R $integrated_data $reference_data cell_types_level1_predicted labeled_data.rds
    """
}

// Workflow definition
workflow mito_filter_processing {
    def data_dir = Channel.fromPath("${params.input_dir}")
    def reference_data = Channel.fromPath('/banach2/SCORCH/analysis/231118-combinedAnalysisPFC/seurat_integrated_v2.RDS')
    sample_metadata = CreateMetaData(data_dir) 
    LoadData(data_dir, sample_metadata) | QualityControl | DoubletDetection | MergeData
    IntegrateData(MergeData.out.merged_data)
    LabelTransfer(IntegrateData.out.integrated_data, reference_data)
}

workflow allen_processing{
    def data_dir = Channel.fromPath("${params.input_dir}")
    def reference_data = Channel.fromPath('/banach2/SCORCH/analysis/231118-combinedAnalysisPFC/seurat_integrated_v2.RDS')
    sample_metadata = CreateMetaData(data_dir) 
    LoadData(data_dir, sample_metadata) | FilterMitoAndSexGenes | QualityControlGenes | DoubletDetection | MergeData
    IntegrateData(MergeData.out.merged_data)
    LabelTransfer(IntegrateData.out.integrated_data, reference_data)
}

workflow {
    allen_processing()
    /*allen_processing()*/
}

workflow testing{
    def data_dir = Channel.fromPath("${params.input_dir}")
    def reference_data = Channel.fromPath('/banach2/SCORCH/analysis/231118-combinedAnalysisPFC/seurat_integrated_v2.RDS')
    sample_metadata = CreateMetaData(data_dir) 
    LoadData(data_dir, sample_metadata)
    FilterMitoAndSexGenes(LoadData.out.raw_samples)
    QualityControlGenes(FilterMitoAndSexGenes.out.no_mito_no_sex_samples)
    DoubletDetection(QualityControlGenes.out.filtered_samples)
    MergeData(DoubletDetection.out.doublet_samples)
    IntegrateData(MergeData.out.merged_data)
}
