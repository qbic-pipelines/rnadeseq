process REPORT {
    //TODO change container?

    input:
    path proj_summary
    path softwareversions
    path model
    path report_options
    path contrnames
    path deseq2
    path multiqc
    path genelist
    path gprofiler
    path quote

    output:
    path "*.zip"
    path "RNAseq_report.html", emit: rnaseq_report

    //TODO: remove multiqc
    script:
    def genelistopt = genelist.name != 'NO_FILE' ? "--genelist $genelist" : ''
    def batchopt = params.batch_effect ? "--batch_effect" : ''
    def quoteopt = quote.name != 'NO_FILE4' ? "$quote" : ''
    def pathwayopt = params.skip_pathway_analysis ? '' : "--pathway_analysis"
 //TODO:   def pa_done = gprofiler.name == "FALSE" ? false: "$gprofiler"
    """
    unzip $deseq2
    unzip $multiqc
    if [ "$gprofiler" != "FALSE" ]; then
        unzip $gprofiler
    fi
    mkdir QC
    mv MultiQC/multiqc_plots/ MultiQC/multiqc_data/ MultiQC/multiqc_report.html QC/
    Execute_report.R --report '$baseDir/assets/RNAseq_report.Rmd' \
    --output 'RNAseq_report.html' \
    --proj_summary $proj_summary \
    --versions $softwareversions \
    --model $model \
    --report_options $report_options \
    --revision $workflow.revision \
    --contrasts $contrnames \
    $genelistopt \
    --organism $params.species \
    --log_FC $params.logFCthreshold \
    $batchopt \
    --min_DEG_pathway $params.min_DEG_pathway \
    $pathwayopt
    if [ "$gprofiler" != "FALSE" ]; then
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ pathway_analysis/ $quoteopt
    else
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ $quoteopt
    fi
    """
}
