process REPORT {

    input:
    path gene_counts
    path metadata
    path model
    path contrast_matrix
    path contrast_list
    path contrast_pairs
    path relevel
    path genelist
    path gtf

    path proj_summary
    path softwareversions
    path multiqc
    path quote

    output:
    path "*.zip"
    path "RNAseq_report.html", emit: rnaseq_report

    script:

    def genelist_opt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def contrast_matrix_opt = contrast_matrix.name != 'DEFAULT' ? "--contrast_matrix $contrast_matrix" : ''
    def contrast_list_opt = contrast_list.name != 'DEFAULT1' ? "--contrast_list $contrast_list" : ''
    def contrast_pairs_opt = contrast_pairs.name != 'DEFAULT2' ? "--contrast_pairs $contrast_pairs" : ''
    def relevel_opt = relevel.name != 'NO_FILE2' ? "--relevel $relevel" : ''
    def batch_effect_opt = params.batch_effect ? "--batchEffect TRUE" : ''
    def rlog_opt = params.use_vst ? '--rlog FALSE' : ''
    def quoteopt = quote.name != 'NO_FILE4' ? "$quote" : ''
    def pathwayopt = params.skip_pathway_analysis ? '' : "--pathway_analysis"

    """
    if [ "$multiqc" != "NO_FILE3" ]; then
        unzip $multiqc
        mkdir QC
        mv MultiQC/multiqc_plots/ MultiQC/multiqc_data/ MultiQC/multiqc_report.html QC/
    fi
    Execute_report.R --report '$baseDir/assets/RNAseq_report.Rmd' \
    --gene_counts $gene_counts \
    --metadata $metadata \
    $contrast_matrix_opt \
    $contrast_list_opt \
    $contrast_pairs_opt \
    $genelist_opt \
    $relevel_opt \
    --output 'RNAseq_report.html' \
    --proj_summary $proj_summary \
    --versions $softwareversions \
    --model $model \
    --revision $workflow.manifest.version \
    $genelist_opt \
    --gtf $gtf \
    --organism $params.organism \
    --log_FC_threshold $params.logFCthreshold \
    $batch_effect_opt \
    --min_DEG_pathway $params.min_DEG_pathway \
    --species_library $params.library \
    --nsub_genes $params.vst_genes_number \
    --keytype $params.keytype \
    --input_type $params.input_type \
    $pathwayopt \
    $rlog_opt
    # Remove allgenes dir as the contained files do not contain only DE genes
    rm -r differential_gene_expression/allgenes
    if [ "$pathwayopt" == "--pathway_analysis" ]; then
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ pathway_analysis/ $quoteopt
    else
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ $quoteopt
    fi
    """
}

