process PATHWAY_ANALYSIS {
    //TODO change container?

    input:
    path deseq_output
    path metadata
    path model
    path genelist
    path keggblacklist

    output:
    path '*.zip', emit: pathway_analysis

    script:
    def genelistopt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def keggblacklistopt = keggblacklist.name != 'NO_FILE3' ? "--kegg_blacklist $keggblacklist" : ''
    def basepath = 'differential_gene_expression/gene_counts_tables/'
    def normInput = params.skip_rlog ? basepath + 'vst_transformed_gene_counts.tsv' : basepath + 'rlog_transformed_gene_counts.tsv'
    """
    unzip $deseq_output
    pathway_analysis.R --dirContrasts 'differential_gene_expression/DE_genes_tables/!(*allgenes*)' --metadata $metadata \
    --model $model --normCounts $normInput \
    --species $params.species $genelistopt $keggblacklistopt --min_DEG_pathway $params.min_DEG_pathway
    zip -r pathway_analysis.zip pathway_analysis/
    """
}
