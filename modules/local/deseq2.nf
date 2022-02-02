process DESEQ2 {
    //TODO change container?
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'qbicpipelines/rnadeseq:1.3.2' :
        'qbicpipelines/rnadeseq:1.3.2' }"

    input:
    path gene_counts
    path metadata
    path model
    path contrast_matrix
    path relevel
    path contrast_list
    path contrast_pairs
    path genelist

    output:
    path '*.zip'             , emit: deseq2
    path "contrast_names.txt", emit: contrnames

    script:
    def gene_list_opt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def contrast_mat_opt = contrast_matrix.name != 'DEFAULT' ? "--contrasts_matrix $contrast_matrix" : ''
    def contrast_list_opt = contrast_list.name != 'DEFAULT1' ? "--contrasts_list $contrast_list" : ''
    def contrast_pairs_opt = contrast_pairs.name != 'DEFAULT2' ? "--contrasts_pairs $contrast_pairs" : ''
    def relevel_opt = relevel.name != 'NO_FILE2' ? "--relevel $relevel" : ''
    def batch_effect_opt = params.batch_effect ? "--batchEffect" : ''
    """
    DESeq2.R --counts $gene_counts --metadata $metadata --design $model \
    --logFCthreshold $params.logFCthreshold $relevel_opt $contrast_mat_opt \
    $contrast_list_opt $contrast_pairs_opt $gene_list_opt $batch_effect_opt
    zip -r differential_gene_expression.zip differential_gene_expression
    """
}
