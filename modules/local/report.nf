process REPORT {
    //TODO change container?

    input:
    path proj_summary
    path softwareversions
    path model
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
    def rlogopt = params.use_vst ? '' : "--rlog"
    """
    unzip $deseq2
    unzip $multiqc
    if [ "$pathwayopt" == "--pathway_analysis" ]; then
        unzip $gprofiler
    fi
    mkdir QC
    mv MultiQC/multiqc_plots/ MultiQC/multiqc_data/ MultiQC/multiqc_report.html QC/
    Execute_report.R --report '$baseDir/assets/RNAseq_report.Rmd' \
    --output 'RNAseq_report.html' \
    --proj_summary $proj_summary \
    --versions $softwareversions \
    --model $model \
    --revision $workflow.revision \
    --contrasts $contrnames \
    $genelistopt \
    --organism $params.organism \
    --log_FC $params.logFCthreshold \
    $batchopt \
    --min_DEG_pathway $params.min_DEG_pathway \
    --species_library $params.library \
    $pathwayopt \
    $rlogopt
    if [ "$pathwayopt" == "--pathway_analysis" ]; then
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ pathway_analysis/ $quoteopt
    else
        zip -r report.zip RNAseq_report.html differential_gene_expression/ QC/ $quoteopt
    fi
    """
}
