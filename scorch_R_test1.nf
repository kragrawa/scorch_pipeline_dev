nextflow.enable.dsl=2

params.output_dir = "/data/kriti/pipeline_dev/samplesheet_test/"
params.input_dir = "/banach2/SCORCH/data/raw/10xMultiome-PFC-HIVOUD_OUD-2pairs-12152022/cellranger_arc"
params.chrX_genes = "/data/kriti/scorch_pipeline_dev/data/X_chromosome_gene_names.txt"
params.chrY_genes = "/data/kriti/scorch_pipeline_dev/data/Y_chromosome_gene_names.txt"
params.label_transfer_colname = "cell_types_level1_predicted"
params.samplesheet = "/data/kriti/scorch_pipeline_dev/data/test_sample_sheet.csv"
params.biccn_reference = "/banach2/SCORCH/data/analysis/resources/BICCN_withMetadata_sct.RDS"
params.biccn_map_col = "within_area_subclass"
params.ma_reference = "/banach2/SCORCH/data/analysis/resources/Ma_Sestan_seuratV5_sct.rds"
params.ma_map_col = "subclass"
params.nhp_vst_reference = "/banach2/SCORCH/data/analysis/resources/primate_nacc_vst_monkeyP_sct.RDS"
params.nhp_vst_map_col = "cell_type_2"
params.vst = false


/*
Sample Sheet and Single Sample Processing
*/
process LoadDataFromSampleSheet {
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}", mode: 'copy'
    
    input:
    tuple val(sample_name), val(cellranger_path), val(metadata)

    output:
    tuple val(sample_name), path("raw_${sample_name}.rds"), emit: raw_data

    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/LoadDataSampleSheet.R "${cellranger_path}" '${metadata}' raw_${sample_name}.rds
    """
}

process FilterMitoAndSexGenesSingleSample {
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}/mito_sex_filter/" , mode: 'copy'
    
    input:
    tuple val(sample_name), val(raw_data)

    output:
    tuple val(sample_name), path ("filtered_${sample_name}.rds"), emit: filtered_data
    path "*.mtx.gz", emit: no_mito_no_sex_counts

    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/SingleSampleGeneFilter.R ${raw_data} filtered_${sample_name}.rds ${params.chrX_genes} ${params.chrY_genes} counts_mitosex_${sample_name}.mtx
    """
}

process FilterNFeatureRNASingleSample {
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}/nRNAFilter/", mode: 'copy'
    
    input:
    tuple val(sample_name), val(filtered_data)

    output:
    tuple val(sample_name), path ("filter_complete_${sample_name}.rds"), emit: rna_filtered_data
    path "*.csv", emit: metadata_nFeatureRNA

    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/nFeatureRNA_SingleSample.R ${filtered_data} filter_complete_${sample_name}.rds metadata_nFeatureRNA_${sample_name}.csv
    """
}

process DoubletDetectionSingleSample {
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}/doublet_detection/" , mode: 'copy'
    
    input:
    tuple val(sample_name), val(rna_filtered_data)

    output:
    tuple val(sample_name), path ("clean_${sample_name}.rds"), emit: no_doublets
    tuple val(sample_name), path ("doublets_${sample_name}.rds"), emit: doublets
    path "metadata_doublet_dectection_${sample_name}.csv", emit: metadata_doublet_detection
    path "embeddings_doublet_dectection_pca_${sample_name}.csv", emit: embeddings_doublet_detection_pca

    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/DoubletDetectionSingleSample.R ${rna_filtered_data} doublets_${sample_name}.rds clean_${sample_name}.rds metadata_doublet_dectection_${sample_name}.csv embeddings_doublet_dectection_pca_${sample_name}.csv
    """
}

process LabelTransferSingleSample{
    maxForks 4
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}/label_transfer/", mode: 'copy'
    
    input:
    tuple val(sample_name), val(no_doublets)

    output:
    tuple val(sample_name), path ("labeled_${sample_name}.rds"), emit: labeled
    path "final_counts_${sample_name}.rds", emit: final_counts
    path "*.png", emit: plots
    path "biccn_pred_${sample_name}.csv", emit: biccn_predictions
    path "ma_pred_${sample_name}.csv", emit: ma_predictions
    path "labeled_metadata_${sample_name}.csv", emit: metadata

 
    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/LabelTransferSingleSample.R ${no_doublets} "${params.biccn_reference}" "${params.biccn_map_col}" "${params.ma_reference}" "${params.ma_map_col}" ${sample_name} labeled_${sample_name}.rds
    """
}

process LabelTransferSingleSampleVST{
    maxForks 1
    tag {sample_name}
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/:/data/kriti/ -v /banach2/SCORCH/data:/banach2/SCORCH/data'
    publishDir "${params.output_dir}/${sample_name}/label_transfer/", mode: 'copy'
    
    input:
    tuple val(sample_name), val(no_doublets)

    output:
    tuple val(sample_name), path ("labeled_${sample_name}.rds"), emit: labeled
    path "*.png", emit: plots
    path "final_counts_${sample_name}.rds", emit: final_counts
    path "biccn_pred_${sample_name}.csv", emit: biccn_predictions
    path "ma_pred_${sample_name}.csv", emit: ma_predictions
    path "nhp_pred_${sample_name}.csv", emit: nhp_predictions
    path "labeled_metadata_${sample_name}.csv", emit: metadata

 
    script:
    """
    Rscript /data/kriti/scorch_pipeline_dev/LabelTransferVST.R ${no_doublets} \
        "${params.biccn_reference}" \
        "${params.biccn_map_col}" \
        "${params.ma_reference}" \
        "${params.ma_map_col}" \
        "${params.nhp_vst_reference}" \
        "${params.nhp_vst_map_col}" \
        ${sample_name} \
        labeled_${sample_name}.rds
    """
}


/* Older Directory Based Processing */
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
    maxForks 4
    memory '32GB'
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


process LabelTransferIndividualSample{
    container 'seurat_v5_image'
    containerOptions = '-v /data/kriti/pipeline_dev:/data/kriti/pipeline_dev'
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path sample_data_list
    path reference_data

    output:
    path "labeled_data_list.rds", emit: labeled_data_list
    path "*.png"
 
    script:
    """
    Rscript /data/kriti/pipeline_dev/scorch_pipeline_dev/LabelTransferSingleSample.R $sample_data_list $reference_data "${params.label_transfer_colname}" labeled_data_list.rds
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

workflow quick_test{
    def data_dir = Channel.fromPath("${params.input_dir}")
    def reference_data = Channel.fromPath('/banach2/SCORCH/analysis/231118-combinedAnalysisPFC/seurat_integrated_v2.RDS')
    sample_metadata = CreateMetaData(data_dir) 
    LoadData(data_dir, sample_metadata)
    FilterMitoAndSexGenes(LoadData.out.raw_samples)
    QualityControlGenes(FilterMitoAndSexGenes.out.no_mito_no_sex_samples)
    DoubletDetection(QualityControlGenes.out.filtered_samples)
    LabelTransferSingleSample(DoubletDetection.out.doublet_samples, reference_data)
    /*MergeData(LabelTransferSingleSample.out.labeled_data_list)
    IntegrateData(MergeData.out.merged_data)*/
}

/* read in the sample sheet*/
samplesheet = Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> tuple(row.sample_name, row.cellranger_path, row.collect { k, v -> "$k=$v" }.join(',')) }

workflow testing{
    LoadDataFromSampleSheet(samplesheet)
    FilterMitoAndSexGenesSingleSample(LoadDataFromSampleSheet.out.raw_data)
    FilterNFeatureRNASingleSample(FilterMitoAndSexGenesSingleSample.out.filtered_data)
    DoubletDetectionSingleSample(FilterNFeatureRNASingleSample.out.rna_filtered_data)
    if(params.vst){
        LabelTransferSingleSampleVST(DoubletDetectionSingleSample.out.no_doublets)
    }
    else{
        LabelTransferSingleSample(DoubletDetectionSingleSample.out.no_doublets)
    }
}