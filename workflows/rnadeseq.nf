/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnadeseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.gene_counts, params.metadata, params.model,
    params.project_summary, params.versions,
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// Check mandatory parameters
if (params.rawcounts) { ch_counts_file = Channel.fromPath(params.rawcounts) } else { exit 1, 'Please provide raw counts file!' }
if (params.metadata) { ch_metadata_file = Channel.fromPath(params.metadata) } else { exit 1, 'Please provide metadata file!' }
if (params.model) { ch_model_file = Channel.fromPath(params.model) } else { exit 1, 'Please provide linear model file!' }
if (params.project_summary) { ch_proj_summary_file = Channel.fromPath(params.project_summary) } else { exit 1, 'Please provide project summary file!' }
if (params.versions) { ch_softwareversions_file = Channel.fromPath(params.versions) } else { exit 1, 'Please provide software versions file!' }


//TODO: Remove _for_blabla??
// Create channels for optional parameters
ch_contrast_matrix = Channel.fromPath(params.contrast_matrix)
ch_contrast_list = Channel.fromPath(params.contrast_list)
ch_contrast_pairs = Channel.fromPath(params.contrast_pairs)
ch_relevel = Channel.fromPath(params.relevel)
ch_quote_file = Channel.fromPath(params.quote)
ch_genes = Channel.fromPath(params.genelist)
ch_report_options_file = Channel.fromPath(params.report_options)
ch_kegg_blacklist = Channel.fromPath(params.kegg_blacklist)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/
//TODO: Are any configs needed?
//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
// TODO: Import here the three processes/modules
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { DESEQ2 } from '../modules/local/deseq2'
include { PATHWAY_ANALYSIS } from '../modules/local/pathway_analysis'
include { REPORT } from '../modules/local/report'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
//TODO: def multiqc_report = []

workflow RNADESEQ {


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    // TODO: Here, I have to call the modules
//
//  MODULE: DE analysis
//
    print "STARTED"
    DESEQ2 (
        ch_counts_file,
        ch_metadata_file,
        ch_model_file,
        ch_contrast_matrix,
        ch_relevel,
        ch_contrast_list,
        ch_contrast_pairs,
        ch_genes
    )
    //TODO: Do I need this part?
    print "alsikfhliahf"
    ch_deseq2 = DESEQ2.out.deseq2
    ch_contrnames = DESEQ2.out.contrnames
   // print ch_deseq2
    print "DE YAY!!!"



//
//  MODULE: Pathway analysis
//
    PATHWAY_ANALYSIS (
        ch_deseq2,
        ch_metadata_file,
        ch_model_file,
        ch_genes,
        ch_kegg_blacklist
    )
    //TODO: Do I need this part?
    ch_pathway_analysis = PATHWAY_ANALYSIS.out.pathway_analysis
    print "PATHWAY YAY!!!"

//
//  MODULE: RNAseq Report
//

    REPORT (
        ch_proj_summary_file,
        ch_softwareversions_file,
        ch_model_file,
        ch_report_options_file,
        ch_contrnames,
        ch_deseq2,
        ch_genes,
        ch_pathway_analysis,     //TODO: remove _for_report? YES
        ch_quote_file
    )
    //TODO: Do I need this part? Delete? Dont add report to email
    ch_rnaseq_report = REPORT.out.rnaseq_report
    print "REPORT YAY!!!"
    // This channel contains the versions of all tools of the current module
    //TODO:
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_softwareversions_file.unique().collectFile(name: 'collated_versions.yml')
    // )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/
//TODO: delete email_on_fail or add it somewhere else to the code with a definition?
//Here I have to change the params?

/*
//TODO: At a later point, remove multiqc and also edit /lib/nfcoretemplate to add deseq report instead of multiqc
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}
*/
/*
========================================================================================
    THE END
========================================================================================
*/
