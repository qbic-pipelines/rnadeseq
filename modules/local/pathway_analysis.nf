process PATHWAY_ANALYSIS {
    //TODO: Is publishdir still used? Or do I have to adapt that somehow? rnaseq does not use it
    //TODO change container

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'qbicpipelines/rnadeseq:1.3.2' :
        'qbicpipelines/rnadeseq:1.3.2' }"

    input:
    path deseq_output
    path metadata
    path model
    path genelist
    path keggblacklist

    output:         //TODO: remove _for_report? YES
    path '*.zip', emit: pathway_analysis

    script:
    def genelistopt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def keggblacklistopt = keggblacklist.name != 'NO_FILE3' ? "--kegg_blacklist $keggblacklist" : ''
    """
    unzip $deseq_output
    pathway_analysis.R --dirContrasts 'differential_gene_expression/DE_genes_tables/' --metadata $metadata \
    --model $model --normCounts 'differential_gene_expression/gene_counts_tables/rlog_transformed_gene_counts.tsv' \
    --species $params.species $genelistopt $keggblacklistopt --min_DEG_pathway $params.min_DEG_pathway
    zip -r pathway_analysis.zip pathway_analysis/
    """
}
