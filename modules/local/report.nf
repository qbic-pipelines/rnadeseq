process REPORT {

    container 'qbicpipelines/rnadeseq:dev'

    input:
    path gene_counts
    path metadata
    path model
    path gtf

    path contrast_matrix
    path contrast_list
    path contrast_pairs
    path genelist
    path relevel

    path proj_summary
    path software_versions
    path multiqc
    path custom_gmt
    path custom_background

    output:
    path "*.zip"
    path "RNAseq_report.html", emit: rnaseq_report

    script:

    def contrast_matrix_opt = contrast_matrix.name != 'DEFAULT' ? "--contrast_matrix $contrast_matrix" : ''
    def contrast_list_opt = contrast_list.name != 'DEFAULT1' ? "--contrast_list $contrast_list" : ''
    def contrast_pairs_opt = contrast_pairs.name != 'DEFAULT2' ? "--contrast_pairs $contrast_pairs" : ''
    def genelist_opt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def relevel_opt = relevel.name != 'NO_FILE2' ? "--relevel $relevel" : ''
    def batch_effect_opt = params.batch_effect ? "--batch_effect TRUE" : ''
    def rlog_opt = params.use_vst ? '--rlog FALSE' : ''
    def round_DE_opt = params.round_DE ? "--round_DE $params.round_DE" : ''

    def pathway_opt = params.run_pathway_analysis ? "--pathway_analysis" : ''
    def custom_gmt_opt = custom_gmt.name != 'NO_FILE3' ? "--custom_gmt $custom_gmt" : ''
    def custom_background_opt = custom_background.name != 'NO_FILE6' ? "--custom_background $custom_background" : ''

    def quote_opt = params.quote != 'NO_FILE5' ? "--path_quote $params.quote" : ''
    def software_versions_opt = params.software_versions != 'NO_FILE6' ? "--software_versions $params.software_versions" : ''

    def citest_opt = params.citest ? "--citest TRUE" : ''

    """
    if [ "$multiqc" != "NO_FILE4" ]; then
        unzip $multiqc
        mkdir QC
        mv MultiQC/multiqc_plots/ MultiQC/multiqc_data/ MultiQC/multiqc_report.html QC/
    fi
    Execute_report.R \
        --report '$baseDir/assets/RNAseq_report.Rmd' \
        --output 'RNAseq_report.html' \
        --input_type $params.input_type \
        --gene_counts $gene_counts \
        --metadata $metadata \
        --model $model \
        --gtf $gtf \
        $contrast_matrix_opt \
        $contrast_list_opt \
        $contrast_pairs_opt \
        $genelist_opt \
        $relevel_opt \
        $batch_effect_opt \
        --logFC_threshold $params.logFC_threshold \
        --pval_threshold $params.pval_threshold \
        $rlog_opt \
        --nsub_genes $params.vst_genes_number \
        $round_DE_opt \
        $pathway_opt \
        $custom_gmt_opt \
        $custom_background_opt \
        --organism $params.organism \
        --species_library $params.species_library \
        --keytype $params.keytype \
        --min_DEG_pathway $params.min_DEG_pathway \
        $quote_opt \
        $software_versions_opt \
        --proj_summary $proj_summary \
        --revision $workflow.manifest.version \
        $citest_opt

    # Remove allgenes dir as the contained files do not contain only DE genes
    rm -r differential_gene_expression/allgenes
    # If citest, copy results before zipping as unzip does not work properly in the container
    if [ "$params.citest" == true ]; then
        mkdir ../../../results_test
        cp -r RNAseq_report.html differential_gene_expression/ ../../../results_test
        if [ "$pathway_opt" == "--pathway_analysis" ]; then
            cp -r pathway_analysis/ ../../../results_test
        fi
    fi
    if [ "$pathway_opt" == "--pathway_analysis" ]; then
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ pathway_analysis/
    else
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/
    fi
    """

}

