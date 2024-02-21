nextflow.enable.dsl=2

process loadData {
    container 'seurat_v5_image'

    input:
    file(input_file) from input_channel
    output:
    file("individual_samples.rds") into output_channel

    // Specify the script to run
    script:
    """
    Rscript your_r_script.R $input_file individual_samples.rds
    """
}

process runRScript {
    container 'my_r_seurat_image'

    output:
        path("output.txt"), emit: result

    script:
        """
        Rscript /data/kriti/pipeline_dev/seurat_test_script.r > output.txt
        """
}

workflow {
    runRScript()
}