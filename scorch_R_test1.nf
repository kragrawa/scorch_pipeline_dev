nextflow.enable.dsl=2

process LoadData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
    input:
    path data_dir
    path sample_metadata

    output:
    path "individual_samples.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/LoadData.R "${data_dir}" $sample_metadata individual_samples.rds
    """
}

process QualityControl {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
    input:
    path raw_samples

    output:
    path "seperate_mito_gene_filtered.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/QualityControl.R $raw_samples seperate_mito_gene_filtered.rds 2
    """
}

process DoubletDetection {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
    input:
    path filtered_samples

    output:
    path "seperate_doublets_filtered.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/DoubleDetection.R $filtered_samples seperate_doublets_filtered.rds
    """
}


process CreateMetaData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
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
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
    input:
    path clean_samples

    output:
    path "merged_data.rds", emit: merged_data
    path "*.png"


    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/MergeData.R $clean_samples merged_data.rds
    """
}

process IntegrateData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    
    input:
    path merged_data

    output:
    path "integrated_data.rds", emit: integrated_data
    path "*.png"
 
    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/IntegrateData.R $merged_data integrated_data.rds
    """
}

// Workflow definition
workflow {
    def data_dir = Channel.fromPath('/banach2/SCORCH/data/raw/10xMultiome-PFC-HIVOUD_OUD-2pairs-12152022/cellranger_v7_RNA/')
    sample_metadata = CreateMetaData(data_dir) 
    LoadData(data_dir, sample_metadata) | QualityControl | DoubletDetection | MergeData
    IntegrateData(MergeData.out.merged_data)
}
