process REPORT {

    container 'ghcr.io/qbic-pipelines/rnadeseq:dev'

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
    path report_file
    path references_file
    path css
    path logo

    output:
    path "*.zip"

    script:

    def contrast_matrix_opt = contrast_matrix.name != 'DEFAULT' ? "--contrast_matrix ${contrast_matrix}" : ''
    def contrast_list_opt = contrast_list.name != 'DEFAULT1' ? "--contrast_list ${contrast_list}" : ''
    def contrast_pairs_opt = contrast_pairs.name != 'DEFAULT2' ? "--contrast_pairs ${contrast_pairs}" : ''
    def genelist_opt = genelist.name != 'NO_FILE' ? "--genelist ${genelist}" : ''
    def relevel_opt = relevel.name != 'NO_FILE2' ? "--relevel ${relevel}" : ''
    def batch_effect_opt = params.batch_effect ? "--batch_effect TRUE" : ''

    def pathway_opt = params.run_pathway_analysis ? "--pathway_analysis" : ''
    def custom_gmt_opt = custom_gmt.name != 'NO_FILE3' ? "--custom_gmt ${custom_gmt}" : ''
    def set_background_opt = params.set_background ? "--set_background TRUE" : "--set_background FALSE"
    def custom_background_opt = custom_background.name != 'NO_FILE7' ? "--custom_background ${custom_background}" : ''
    def datasources_opt = params.datasources ? "--datasources ${params.datasources}" : ''
    def heatmaps_cluster_rows_opt = params.heatmaps_cluster_rows ? "--heatmaps_cluster_rows TRUE" : ''
    def heatmaps_cluster_cols_opt = params.heatmaps_cluster_cols ? "--heatmaps_cluster_cols TRUE" : ''
    def pathway_adj_pval_threshold_opt = params.pathway_adj_pval_threshold == -1 ? "--pathway_adj_pval_threshold ${params.adj_pval_threshold}" : "--pathway_adj_pval_threshold ${params.pathway_adj_pval_threshold}"


    def quote_opt = params.quote != 'NO_FILE5' ? "--path_quote ${params.quote}" : ''
    def software_versions_opt = params.software_versions != 'NO_FILE6' ? "--software_versions ${params.software_versions}" : ''

    def citest_opt = params.citest ? "--citest TRUE" : ''

    """
    if [ "${multiqc}" != "NO_FILE4" ]; then
        unzip ${multiqc}
        mkdir QC
        mv MultiQC/multiqc_plots/ MultiQC/multiqc_data/ MultiQC/multiqc_report.html QC/ || mv multiqc/*/multiqc_plots/ multiqc/*/multiqc_data/ multiqc/*/multiqc_report.html QC/ || mv multiqc_plots/ multiqc_data/ multiqc_report.html QC/
    fi
    Execute_report.R \
        --report '${report_file}' \
        --output 'rnadeseq_report.html' \
        --input_type ${params.input_type} \
        --gene_counts ${gene_counts} \
        --metadata ${metadata} \
        --model ${model} \
        --gtf ${gtf} \
        ${contrast_matrix_opt} \
        ${contrast_list_opt} \
        ${contrast_pairs_opt} \
        ${genelist_opt} \
        ${relevel_opt} \
        ${batch_effect_opt} \
        --logFC_threshold ${params.logFC_threshold} \
        --adj_pval_threshold ${params.adj_pval_threshold} \
        --norm_method ${params.norm_method} \
        --nsub_genes ${params.vst_genes_number} \
        --round_DE ${params.round_DE} \
        ${pathway_opt} \
        ${custom_gmt_opt} \
        ${set_background_opt} \
        ${custom_background_opt} \
        --organism ${params.organism} \
        --species_library ${params.species_library} \
        --keytype ${params.keytype} \
        --min_DEG_pathway ${params.min_DEG_pathway} \
        ${datasources_opt} \
        ${heatmaps_cluster_rows_opt} \
        ${heatmaps_cluster_cols_opt} \
        ${pathway_adj_pval_threshold_opt} \
        ${quote_opt} \
        ${software_versions_opt} \
        --proj_summary ${proj_summary} \
        --revision ${workflow.manifest.version} \
        --logo ${logo} \
        ${citest_opt}

    # Remove allgenes dir as the contained files do not contain only DE genes
    rm -r differential_gene_expression/allgenes
    # If citest, copy results before zipping as unzip does not work properly in the container
    if [ "${params.citest}" == true ]; then
        mkdir ../../../results_test
        cp -r rnadeseq_report.html differential_gene_expression/ ../../../results_test
        if [ "${pathway_opt}" == "--pathway_analysis" ]; then
            cp -r enrichment_analysis/ ../../../results_test
        fi
    fi
    if [ "${pathway_opt}" == "--pathway_analysis" ]; then
        zip -r report.zip rnadeseq_report.html differential_gene_expression/ QC/ enrichment_analysis/ rnadeseq_software_versions.yml
    else
        zip -r report.zip rnadeseq_report.html differential_gene_expression/ QC/ rnadeseq_software_versions.yml
    fi
    """
}
