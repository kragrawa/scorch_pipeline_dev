nextflow.enable.dsl=2

process loadData {
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir '/data/kriti/pipeline_dev/test_output/', mode: 'copy'
    

    input:
    path data_dir

    output:
    path "individual_samples.rds"

    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/LoadData.R "${data_dir}" individual_samples.rds
    """
}

// Workflow definition
workflow {
    def data_dir = Channel.fromPath('/banach2/SCORCH/data/raw/10xMultiome-PFC-HIVOUD_OUD-2pairs-12152022/cellranger_v7_RNA/')
    loadData(data_dir)
}