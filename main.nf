#!/usr/bin/env nextflow
/*
========================================================================================
    qbic-pipelines/rnadeseq
========================================================================================
    Github : https://github.com/qbic-pipelines/rnadeseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
 * Create a channel for input  files
 */

// Required parameters
Channel.fromPath("${params.rawcounts}", checkIfExists: true)
           .ifEmpty{exit 1, "Please provide raw counts file!"}
           .set {ch_counts_file}
Channel.fromPath("${params.metadata}", checkIfExists: true)
           .ifEmpty{exit 1, "Please provide metadata file!"}
           .into { ch_metadata_file_for_deseq2; ch_metadata_file_for_pathway }
Channel.fromPath("${params.model}", checkIfExists: true)
            .ifEmpty{exit 1, "Please provide linear model file!"}
            .into { ch_model_for_deseq2_file; ch_model_for_report_file; ch_model_file_for_pathway}
Channel.fromPath("${params.project_summary}", checkIfExists: true)
            .ifEmpty{exit 1, "Please provide project summary file!"}
            .set { ch_proj_summary_file }
Channel.fromPath("${params.versions}", checkIfExists: true)
            .ifEmpty{exit 1, "Please provide sofware versions file!"}
            .set { ch_softwareversions_file }
Channel.fromPath("${params.multiqc}", checkIfExists: true)
            .ifEmpty{exit 1, "Please provide multiqc.zip folder!"}
            .set { ch_multiqc_file }

// Optional parameters
Channel.fromPath("${params.contrast_matrix}")
            .set { ch_contrast_matrix_for_deseq2 }
Channel.fromPath("${params.contrast_list}")
            .set { ch_contrast_list_for_deseq2 }
Channel.fromPath("${params.contrast_pairs}")
            .set { ch_contrast_pairs_for_deseq2 }
Channel.fromPath("${params.relevel}")
            .set { ch_relevel_for_deseq2 }
Channel.fromPath("${params.quote}")
           .set { ch_quote_file }
Channel.fromPath("${params.genelist}")
            .into { ch_genes_for_deseq2_file; ch_genes_for_report_file; ch_genes_for_pathway }
Channel.fromPath("${params.report_options}", checkIfExists: true)
            .set { ch_report_options_file }
Channel.fromPath("${params.kegg_blacklist}")
            .set { ch_kegg_blacklist_for_pathway }



ch_fastqc_file = file(params.fastqc)

/*
 * Check mandatory parameters
 */
if (!params.species) {
  exit 1, "No species has been specified!"
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Gene Counts'] = params.rawcounts
summary['Metadata'] = params.metadata
summary['Model'] = params.model
summary['Contrast matrix'] = params.contrast_matrix
summary['Contrast list'] = params.contrast_list
summary['Contrast pairs'] = params.contrast_pairs
summary['Relevel'] = params.relevel
summary['LogFCthreshold'] = params.logFCthreshold
summary['Gene list'] = params.genelist
summary['Project summary'] = params.project_summary
summary['Software versions'] = params.versions
summary['Report options'] = params.report_options
summary['Multiqc results'] = params.multiqc
summary['Species'] = params.species
summary['Min DEG gene'] = params.min_DEG_pathway
summary['Quote'] = params.quote
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'qbic-pipelines-rnadeseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'qbic-pipelines/rnadeseq Workflow Summary'
    section_href: 'https://github.com/qbic-pipelines/rnadeseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.tsv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo \$(R --version 2>&1) > v_R.txt
    Rscript -e "library(RColorBrewer); write(x=as.character(packageVersion('RColorBrewer')), file='v_rcolorbrewer.txt')"
    Rscript -e "library(reshape2); write(x=as.character(packageVersion('reshape2')), file='v_reshape2.txt')"
    Rscript -e "library(genefilter); write(x=as.character(packageVersion('genefilter')), file='v_genefilter.txt')"
    Rscript -e "library(DESeq2); write(x=as.character(packageVersion('DESeq2')), file='v_deseq2.txt')"
    Rscript -e "library(ggplot2); write(x=as.character(packageVersion('ggplot2')), file='v_ggplot2.txt')"
    Rscript -e "library(plyr); write(x=as.character(packageVersion('plyr')), file='v_plyr.txt')"
    Rscript -e "library(vsn); write(x=as.character(packageVersion('vsn')), file='v_vsn.txt')"
    Rscript -e "library(gplots); write(x=as.character(packageVersion('gplots')), file='v_gplots.txt')"
    Rscript -e "library(pheatmap); write(x=as.character(packageVersion('pheatmap')), file='v_pheatmap.txt')" 
    Rscript -e "library(optparse); write(x=as.character(packageVersion('optparse')), file='v_optparse.txt')"
    Rscript -e "library(svglite); write(x=as.character(packageVersion('svglite')), file='v_svglite.txt')"
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - DE analysis
 */
process DESeq2 {
    publishDir "${params.outdir}/differential_gene_expression", mode: 'copy'

    input:
    file(gene_counts) from ch_counts_file
    file(metadata) from ch_metadata_file_for_deseq2
    file(model) from ch_model_for_deseq2_file
    file(contrast_matrix) from ch_contrast_matrix_for_deseq2
    file(relevel) from ch_relevel_for_deseq2
    file(contrast_list) from ch_contrast_list_for_deseq2
    file(contrast_pairs) from ch_contrast_pairs_for_deseq2
    file(genelist) from ch_genes_for_deseq2_file

    output:
    file "*.zip" into ch_deseq2_for_report, ch_deseq2_for_pathway
    file "contrast_names.txt" into ch_contrnames_for_report

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

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
WorkflowMain.initialise(workflow, params, log)
/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { RNADESEQ } from './workflows/rnadeseq'
//
// WORKFLOW: Run main qbic-pipelines/rnadeseq analysis pipeline
//
workflow QBIC_RNADESEQ {
    RNADESEQ ()
}
/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    QBIC_RNADESEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
