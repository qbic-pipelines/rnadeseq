process PATHWAY_ANALYSIS {

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    deseq_output
    metadata
    model
    genelist
    keggblacklist

    output:         //TODO: remove _for_report?
    path '*.zip'       , emit: pathway_analysis_for_report

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
